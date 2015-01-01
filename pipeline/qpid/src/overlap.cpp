#include "AMOS/src/foundation_AMOS.hh"
#include "overlap.h"
#include "parser/parser.h"
#include "timer/timer.h"
#include "memory/memory.cpp"
#include "align/align.h"
#include "parsero/parsero.h"
#include "read.h"

#include <algorithm>
#include <mutex>
#include <vector>
#include <string>
#include <utility>
#include <zlib.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
using std::vector;
using std::string;
using std::pair;
using std::swap;

typedef unsigned int uint;
typedef void (*output_funptr)(const Overlap& overlap);

// init fasta/fastq reader
KSEQ_INIT(gzFile, gzread)

int THREADS_NUM = 4;
int ALIGNMENT_BAND_RADIUS = 5;
int OFFSET_WIGGLE = 3;
int MERGE_RADIUS = 5 * ALIGNMENT_BAND_RADIUS;
double MAXIMUM_ERROR_RATE = 0.03;

char *INPUT_FILE = NULL;
FILE *OUTPUT_FD = stdout;

char *BANK_INPUT = NULL;
char *BANK_OUTPUT = NULL;
AMOS::BankStream_t *bank_output = NULL;
std::mutex bank_mutex;

output_funptr output;

ThreadPool* pool;

void usage(int argc, char **argv) {
  fprintf(stderr, "usage:\n");
  fprintf(stderr, "\t%s [-f <fastq_reads>] [-O <output_overlaps_afg>] [-b <amos_input_bank>] [-B <amos_output_bank>]\n", argv[0]);
}

bool sort_offsets(offset_t a, offset_t b) {
    if (a.index == b.index) {
        return a.lo_offset < b.lo_offset;
    }
    return a.index < b.index;
}

void read_from_fasta(vector<Read>& reads, const char *filename) {
    Timer* timer = new Timer("reading");

    fprintf(stderr, "* Reading from file %s...\n", filename);

    gzFile fp = gzopen(filename, "r");      // STEP 2: open the file handler
    kseq_t *seq = kseq_init(fp);            // STEP 3: initialize seq

    int len = 0;
    int id = 0;
    while ((len = kseq_read(seq)) > 0) {    // STEP 4: read sequence
        char *read_string = new char[len + 1];
        strcpy(read_string, seq->seq.s);
        reads.push_back(Read(id++, read_string));
    }

    kseq_destroy(seq);      // STEP 5: destroy seq
    gzclose(fp);            // STEP 6: close the file handler

    timer->end(true);
    delete timer;
}

void read_from_bank(vector<Read>& reads, const char *bank_name) {

  AMOS::BankStream_t bank(AMOS::Read_t::NCODE);
  assert(bank.exists(bank_name));

  bank.open(bank_name, AMOS::B_READ);

  AMOS::Read_t read;
  bank.seekg(0, AMOS::BankStream_t::BEGIN);
  while (bank >> read) {
    AMOS::Range_t clear = read.getClearRange();
    string seq = read.getSeqString(clear);
    char* cpy = new char[seq.length() + 1];
    strcpy(cpy, seq.c_str());
    int id = read.getIID();
    reads.push_back(Read(id, cpy));
  }

  bank.close();
  return;
}

void output_overlap_file(const Overlap& overlap) {
    fprintf(OUTPUT_FD, "{OVL\nadj:%c\nrds:%d,%d\nscr:%d\nahg:%d\nbhg:%d\n}\n",
        overlap.normal_overlap ? 'N' : 'I',
        overlap.r1.id,
        overlap.r2.id,
        (int) overlap.score,
        overlap.a_hang,
        overlap.b_hang
   );
}

void output_overlap_bank(const Overlap& overlap) {
  AMOS::Overlap_t amos_overlap;
  std::pair<AMOS::ID_t, AMOS::ID_t> read_pair(overlap.r1.id, overlap.r2.id);

  if (overlap.normal_overlap) {
    amos_overlap.setAdjacency(AMOS::Overlap_t::NORMAL);
  } else {
    amos_overlap.setAdjacency(AMOS::Overlap_t::INNIE);
  }
  amos_overlap.setReads(read_pair);
  amos_overlap.setAhang(overlap.a_hang);
  amos_overlap.setBhang(overlap.b_hang);

  bank_mutex.lock();
  (*bank_output) << amos_overlap;
  bank_mutex.unlock();
}

char base_complement(char base) {
    if (base == 'A')        return 'T';
    else if (base == 'C')   return 'G';
    else if (base == 'G')   return 'C';
    else                    return 'A';
}

const Read reversed_complement(const Read& read) {

    int len = strlen(read.sequence);
    char* result = new char[len + 1];
    result[len] = 0;

    for (int i = 0; i < len; ++i) {
        result[i] = base_complement(read.sequence[len - i - 1]);
    }

    return Read(read.id, result);
}

void add_offset(vector<offset_t>& offsets, const int& str_index, const int& offset, const int& wiggle) {

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

void find_overlaps_from_offsets(vector<Read>& reads, int t, const Read &target, vector<offset_t>& offsets,
    bool target_forward_oriented) {

    int len_t = strlen(target.sequence);

    int i = 0, offsets_len = offsets.size();
    while (i < offsets_len) {
      Overlap best_overlap;
      std::pair<int, int> start, end;
      int q = offsets[i].index;
      int len_q = strlen(reads[q].sequence);

      int first_j = i;
      for (int j = i; j < offsets_len && offsets[j].index == q; ++j, ++i) {
        offset_t& offset = offsets[j];
        int score = banded_overlap(
            target.sequence,
            len_t,
            reads[q].sequence,
            len_q,
            offset.lo_offset - ALIGNMENT_BAND_RADIUS,
            offset.hi_offset + ALIGNMENT_BAND_RADIUS,
            &start,
            &end
        );

        Overlap curr_overlap(reads[t], reads[q], score, start, end, target_forward_oriented, true);

        if (j == first_j || score > best_overlap.score) {
          best_overlap = curr_overlap;
        }
      }

      if (best_overlap.error_rate < MAXIMUM_ERROR_RATE) {
        output(best_overlap);
      }
    }
}

void find_overlaps(vector<Read>& reads, Minimizer *minimizer, int wiggle, int merge_radius, bool forward_overlaps = true) {

    const minimizers_t& minimizers = minimizer->get_minimizers();

    vector<std::future<void>> results;
    for (int t = 0, tlen = reads.size(); t < tlen; ++t) {

        results.push_back(pool->enqueue([&reads, forward_overlaps, wiggle, merge_radius, &minimizer, &minimizers, t]() {

            const Read &target = forward_overlaps ? reads[t] : reversed_complement(reads[t]);
            vector<minimizer_t> curr_minimizers;
            vector<offset_t> curr_offsets;

            minimizer->calculate_and_get(curr_minimizers, target.sequence);

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

            find_overlaps_from_offsets(reads, t, target, curr_offsets, forward_overlaps);
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

    parsero::add_option("O:", "output file",
        [] (char *filename) { OUTPUT_FD = fopen(filename, "w"); }
    );

    parsero::add_option("B:", "output bank",
        [] (char *dirname) { BANK_OUTPUT = dirname; }
    );

    parsero::add_option("b:", "input bank",
        [] (char *dirname) { BANK_INPUT = dirname; }
    );

    parsero::add_option("f:", "input file",
        [] (char *filename) { INPUT_FILE = filename; }
    );

    parsero::parse(argc, argv);
}

int main(int argc, char **argv) {

    setup_cmd_interface(argc, argv);

    if (BANK_INPUT == nullptr && INPUT_FILE == nullptr) {
      usage(argc, argv);
      exit(1);
    }

    vector<Read> reads;

    // initialize a thread pool used for finding overlaps
    pool = new ThreadPool(THREADS_NUM);

    // read rads
    if (BANK_INPUT == NULL) {
      // read the data
      read_from_fasta(reads, INPUT_FILE);
    } else {
      // read the data
      read_from_bank(reads, BANK_INPUT);
    }

    fprintf(stderr, "* Read %lu strings...\n", reads.size());

    // setup output
    if (BANK_OUTPUT == NULL) {
      output = &output_overlap_file;
    } else {
      bank_output = new AMOS::BankStream_t(AMOS::Overlap_t::NCODE);
      if (!bank_output->exists(BANK_OUTPUT)) {
        fprintf(stderr, "Output bank did not exist. Creating bank %s\n", BANK_OUTPUT);
        bank_output->create(BANK_OUTPUT);
      }
      bank_output->open(BANK_OUTPUT, AMOS::B_WRITE);
      output = &output_overlap_bank;
    }

    fprintf(stderr, "* Maximum error rate: %lf\n", MAXIMUM_ERROR_RATE);

    Timer mtimer("calculating minimizers");
    // create a bank of all minimizers so finding appropriate read pairs could be efficient.
    Minimizer *m = new Minimizer(16, 20);
    for (int i = 0, len = reads.size(); i < len; ++i) {
      m->calculate_and_store(i, reads[i].sequence);
    }
    mtimer.end();

    try {
      Timer ftimer("calculating forward overlaps");
      find_overlaps(reads, m, OFFSET_WIGGLE, MERGE_RADIUS, true);
      ftimer.end();

      Timer btimer("calculating backward overlaps");
      find_overlaps(reads, m, OFFSET_WIGGLE, MERGE_RADIUS, false);
      btimer.end();
    } catch (AMOS::Exception_t &e) {
      fprintf(stderr, "%s:%d %s", e.file(), e.line(), e.what());
      exit(1);
    }

    // cleaning up the mess
    for (int i = 0, len = reads.size(); i < len; ++i) {
        delete[] reads[i].sequence;
    }

    delete m;
    delete pool;

    fclose(OUTPUT_FD);

    if (bank_output != NULL) {
      bank_output->close();
    }

    return 0;
}
