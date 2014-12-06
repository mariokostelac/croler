#include "./overlap.h"
#include <zlib.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include "AMOS/src/foundation_AMOS.hh"
#include "./parser/parser.h"
#include "./timer/timer.h"
#include "./memory.cpp"
#include "./align/align.h"
#include "./parsero/parsero.h"
using std::vector;
using std::string;

// init fasta/fastq reader
KSEQ_INIT(gzFile, gzread)

int THREADS_NUM = 4;
int ALIGNMENT_BAND_RADIUS = 5;
int OFFSET_WIGGLE = 3;
int MERGE_RADIUS = 5 * ALIGNMENT_BAND_RADIUS;
double MAXIMUM_ERROR_RATE = 0.05;

char *INPUT_FILE = NULL;
FILE *OUTPUT_FD = stdout;

typedef unsigned int uint;

ThreadPool* pool;

bool sort_offsets(offset_t a, offset_t b) {
    if (a.index == b.index) {
        return a.lo_offset < b.lo_offset;
    }
    return a.index < b.index;
}

void read_from_bank(vector<char *>& string_list, AMOS::BankStream_t& read_bank) {
  AMOS::Read_t read;
  read_bank.seekg(0, AMOS::BankStream_t::BEGIN);
  while (read_bank >> read) {
    AMOS::Range_t clear = read.getClearRange();
    string seq = read.getSeqString(clear);
    //clr_list.push_back(clear);
    //s.push_back(tmp);
  }
  return;
}

void read_from_fasta(vector<char *>& string_list, const char *filename) {
    Timer* timer = new Timer("reading");

    fprintf(stderr, "* Reading from file %s...\n", filename);

    gzFile fp = gzopen(filename, "r");      // STEP 2: open the file handler
    kseq_t *seq = kseq_init(fp);            // STEP 3: initialize seq

    int len = 0;
    while ((len = kseq_read(seq)) > 0) {    // STEP 4: read sequence
        char *read_string = new char[len + 1];
        strcpy(read_string, seq->seq.s);
        string_list.push_back(read_string);
    }

    kseq_destroy(seq);      // STEP 5: destroy seq
    gzclose(fp);            // STEP 6: close the file handler

    fprintf(stderr, "* Read %lu strings...\n", string_list.size());

    timer->end(true);
    delete timer;
}

void output_overlap(int index1, int index2, int len1, int len2, int score, std::pair<int, int>& start, std::pair<int, int>& end,
        double errate, bool forward_overlap) {

    int a_hang = abs(len1 - end.first - start.first);
    int b_hang = abs(len2 - end.second - start.second);
    if (start.first > 0) {
        a_hang *= -1;
        b_hang *= -1;
    }
    fprintf(OUTPUT_FD, "{OVL adj:%c rds:%d,%d scr:%d ahg:%d bhg:%d }\n",
        forward_overlap ? 'N' : 'I',
        index1 + 1,
        index2 + 1,
        score,
        a_hang,
        b_hang
   );
}

char base_complement(char base) {
    if (base == 'A')        return 'T';
    else if (base == 'C')   return 'G';
    else if (base == 'G')   return 'C';
    else                    return 'A';
}

char* reversed_complement(const char *seq) {

    int len = strlen(seq);
    char* result = new char[len + 1];
    result[len] = 0;

    for (int i = 0; i < len; ++i) {
        result[i] = base_complement(seq[len - i - 1]);
    }

    return result;
}

void add_offset(vector<offset_t>& offsets, unsigned int str_index, int offset, int wiggle) {

    // extend existing offset if it is possible
    for (int i = 0, len = offsets.size(); i < len; ++i) {
        if (offsets[i].index != str_index) continue;

        // merge it
        if (offsets[i].lo_offset - wiggle <= offset && offsets[i].hi_offset + wiggle >= offset) {
            offsets[i].lo_offset = std::min(offsets[i].lo_offset, offset);
            offsets[i].hi_offset = std::max(offsets[i].hi_offset, offset);
            return;
        }
    }

    // not merged, add it as standalone
    offsets.push_back(offset_t(str_index, offset));
}

void merge_offsets(vector<offset_t>& offsets, uint radius) {

    if (offsets.size() == 0) return;

    // index of current offset
    uint curr = 0;
    for (uint i = 1, len = offsets.size(); i < len; ++i) {
        if (offsets[curr].index == offsets[i].index && offsets[curr].hi_offset + radius >= offsets[i].lo_offset - radius) {
            // extend offset to the right
            offsets[curr].hi_offset = std::max(offsets[curr].hi_offset, offsets[i].hi_offset);
        } else {
            // finish current offset
            curr++;
            if (curr < i) offsets[curr] = offsets[i];
        }
    }

    offsets.resize(curr + 1);
}

void find_overlaps_from_offsets(vector<char *>& string_list, int t, const char *target, vector<offset_t>& offsets,
        overlap_band_t* band = NULL, bool forward_overlaps = true) {

    if (offsets.size() == 0) return;

    int len_t = strlen(target);

    // variables used for tracking the best overlap for current pair
    int last_q = offsets[0].index;
    int best_q = -1, best_len_q = -1, best_score = -1;
    double best_errate = 0;
    std::pair<int, int> best_start, best_end;

    for (int j = 0, jlen = offsets.size(); j < jlen; ++j) {

        auto& offset = offsets[j];
        int q = offset.index;
        int len_q = strlen(string_list[q]);

        std::pair<int, int> start, end;
        int score = banded_overlap(target, len_t, string_list[q], len_q,
                offset.lo_offset - ALIGNMENT_BAND_RADIUS, offset.hi_offset + ALIGNMENT_BAND_RADIUS, &start, &end, band);

        double len = (end.first - start.first + end.second - start.second) / 2.;
        double errors = (score - len)/(INDEL_SCORE + GAP_SCORE + MISMATCH_SCORE);
        double error_rate = errors / len;

        // output best overlap for previous pair (t, q)
        if (q != last_q) {
            if (best_errate < MAXIMUM_ERROR_RATE) {
                output_overlap(t, best_q, len_t, best_len_q, best_score, best_start, best_end, best_errate, forward_overlaps);
            }
            best_score = -1;
            best_errate = 1000;
        }

        // set best score for current (t, q)
        if (best_q == -1 || score > best_score) {
            best_q = q;
            best_len_q = len_q;
            best_score = score;
            best_start = start;
            best_end = end;
            best_errate = error_rate;
        }

        last_q = q;
    }

    // output best overlap for last pair (t, q)
    if (best_errate < MAXIMUM_ERROR_RATE) {
        output_overlap(t, best_q, len_t, best_len_q, best_score, best_start, best_end, best_errate, forward_overlaps);
    }
}

void find_overlaps(vector<char *>& string_list, Minimizer *minimizer, int wiggle, int merge_radius, bool forward_overlaps = true) {

    // find longest string so we could create band
    unsigned long max_len = 0;
    for (int i = 0, len = string_list.size(); i < len; ++i) {
        max_len = std::max(max_len, strlen(string_list[i]));
    }

    const minimizers_t& minimizers = minimizer->get_minimizers();

    vector<std::future<void>> results;
    for (int t = 0, tlen = string_list.size(); t < tlen; ++t) {

        results.push_back(pool->enqueue([&string_list, forward_overlaps, max_len, wiggle, merge_radius, &minimizer, &minimizers, t]() {

            char *target = forward_overlaps ? string_list[t] : reversed_complement(string_list[t]);
            overlap_band_t band(max_len + 1, max_len + 1);
            vector<minimizer_t> curr_minimizers;
            vector<offset_t> curr_offsets;

            minimizer->calculate_and_get(curr_minimizers, target);

            for (uint m = 0, mlen = curr_minimizers.size(); m < mlen; ++m) {
                auto list = minimizers.get_list(curr_minimizers[m].str);

                for (auto kp = list.begin(); kp != list.end(); ++kp) {
                    int k = (*kp).first;
                    if (t >= k) continue;

                    add_offset(curr_offsets, k, (*kp).second - curr_minimizers[m].pos, wiggle);
                }
            }

            std::sort(curr_offsets.begin(), curr_offsets.end(), sort_offsets);
            merge_offsets(curr_offsets, merge_radius);

            find_overlaps_from_offsets(string_list, t, target, curr_offsets, &band, forward_overlaps);

            if (forward_overlaps == false) {
                delete[] target;
            }
        }));
    }

    // wait for the results
    for (int i = 0, len = results.size(); i < len; ++i) {
        results[i].get();
    }
}

void setup_cmd_interface(int argc, char **argv) {

    parsero::add_option("t:", "number of threads",
        [] (char *option) { THREADS_NUM = atoi(option);
    });

    parsero::add_option("a:", "alignment band radius",
        [] (char *option) { ALIGNMENT_BAND_RADIUS = atoi(option); }
    );

    parsero::add_option("w:", "offset wiggle",
        [] (char *option) { OFFSET_WIGGLE = atoi(option); }
    );

    parsero::add_option("r:", "merge radius",
        [] (char *option) { MERGE_RADIUS = atoi(option); }
    );

    parsero::add_option("e:", "maximum error rate",
        [] (char *option) { sscanf(option, "%lf", &MAXIMUM_ERROR_RATE); }
    );

    parsero::add_option("o:", "output file",
        [] (char *filename) { OUTPUT_FD = fopen(filename, "w"); }
    );

    parsero::add_argument("input file",
        [] (char *filename) { INPUT_FILE = filename; }
    );

    parsero::parse(argc, argv);
}

int main(int argc, char **argv) {

    setup_cmd_interface(argc, argv);

    // initialize a thread pool used for finding overlaps
    pool = new ThreadPool(THREADS_NUM);

    vector<char *> string_list;
    read_from_fasta(string_list, INPUT_FILE);

    Minimizer *m = new Minimizer(16, 20);

    Timer mtimer("calculating minimizers");
    for (int i = 0, len = string_list.size(); i < len; ++i) {
        m->calculate_and_store(i, string_list[i]);
    }
    mtimer.end();

    Timer ftimer("calculating forward overlaps");
    find_overlaps(string_list, m, OFFSET_WIGGLE, MERGE_RADIUS, true);
    ftimer.end();

    Timer btimer("calculating backward overlaps");
    find_overlaps(string_list, m, OFFSET_WIGGLE, MERGE_RADIUS, false);
    btimer.end();

    // cleaning up the mess
    for (int i = 0, len = string_list.size(); i < len; ++i) {
        delete[] string_list[i];
    }

    delete m;
    delete pool;

    fclose(OUTPUT_FD);

    return 0;
}
