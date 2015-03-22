#include "overlap.h"
#include "timer/timer.h"
#include "memory/memory.cpp"
#include "align/align.h"
#include "read.h"
#include "lib/amos/reader.cpp"
#include "lib/parsero/parsero.h"

#include <algorithm>
#include <mutex>
#include <vector>
#include <string>
#include <utility>
#include <unistd.h>
#include <zlib.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
using std::vector;
using std::string;
using std::pair;
using std::swap;

typedef unsigned int uint;

int THREADS_NUM = sysconf(_SC_NPROCESSORS_ONLN);
int ALIGNMENT_BAND_RADIUS = 5;
int OFFSET_WIGGLE = 3;
int MERGE_RADIUS = 5 * ALIGNMENT_BAND_RADIUS;
double MAXIMUM_ERROR_RATE = 0.03;

char *INPUT_FILE = NULL;
FILE *OUTPUT_FD = stdout;

ThreadPool* pool;

bool sort_offsets(offset_t a, offset_t b) {
    if (a.index == b.index) {
        return a.lo_offset < b.lo_offset;
    }
    return a.index < b.index;
}

int read_from_afg(vector<Read>& reads, const char *filename) {
    Timer* timer = new Timer("reading");
    fprintf(stderr, "* Reading from file %s...\n", filename);

    vector<AMOS::Read*> tmp_reads;
    int reads_size = AMOS::get_reads(tmp_reads, filename);
    for (int i = 0; i < reads_size; ++i) {
      const auto& r = tmp_reads[i];
      int r_len = r->clr_hi - r->clr_lo;

      char* cpy = new char[r_len + 1];
      strncpy(cpy, r->seq + r->clr_lo, r_len);
      cpy[r_len] = 0; // terminate

      reads.push_back(Read(r->iid, cpy));
      delete tmp_reads[i];
    }

    timer->end(true);
    delete timer;

    return reads_size;
}

void output_overlap(const Overlap& overlap) {
    fprintf(OUTPUT_FD, "{OVL\nadj:%c\nrds:%d,%d\nscr:%d\nahg:%d\nbhg:%d\n}\n",
        overlap.normal_overlap ? 'N' : 'I',
        overlap.r1.id,
        overlap.r2.id,
        (int) overlap.score,
        overlap.a_hang,
        overlap.b_hang
   );
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
        output_overlap(best_overlap);
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

            // delete reversed complement from memory
            if (!forward_overlaps) {
              delete[] target.sequence;
            }
        }));
    }

    // wait for the results
    for (int i = 0, len = results.size(); i < len; ++i) {
        results[i].get();
    }
}

void setup_cmd_interface(int argc, char **argv) {

  parsero::set_header("qpid if read overlapper, often used as a part of croler genome assembler.");
  parsero::set_footer("Found a bug? Drop an email to mario.kostelac@gmail.com.");

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

  parsero::add_option("o:", "output file; if omitted, goes to stdout",
      [] (char *filename) { OUTPUT_FD = fopen(filename, "w"); }
      );

  parsero::add_argument("reads.afg",
      [] (char *filename) { INPUT_FILE = filename; }
      );

  parsero::parse(argc, argv);
}

int main(int argc, char **argv) {

    setup_cmd_interface(argc, argv);

    if (INPUT_FILE == nullptr) {
      parsero::help(argv[0]);
      exit(1);
    }

    vector<Read> reads;

    // initialize a thread pool used for finding overlaps
    pool = new ThreadPool(THREADS_NUM);

    int reads_size = read_from_afg(reads, INPUT_FILE);
    fprintf(stderr, "* Read %d strings...\n", reads_size);

    fprintf(stderr, "* Alignment band radius: %d\n", ALIGNMENT_BAND_RADIUS);
    fprintf(stderr, "* Number of threads: %d\n", THREADS_NUM);
    fprintf(stderr, "* Maximum error rate: %lf\n", MAXIMUM_ERROR_RATE);
    fprintf(stderr, "* Merge radius: %d\n", MERGE_RADIUS);
    fprintf(stderr, "* Offset wiggle: %d\n", OFFSET_WIGGLE);

    Timer mtimer("calculating minimizers");
    // create a bank of all minimizers so finding appropriate read pairs could be efficient.
    Minimizer *m = new Minimizer(16, 20);
    for (int i = 0, len = reads.size(); i < len; ++i) {
      m->calculate_and_store(i, reads[i].sequence);
    }
    mtimer.end();

    Timer ftimer("calculating forward overlaps");
    find_overlaps(reads, m, OFFSET_WIGGLE, MERGE_RADIUS, true);
    ftimer.end();

    Timer btimer("calculating backward overlaps");
    find_overlaps(reads, m, OFFSET_WIGGLE, MERGE_RADIUS, false);
    btimer.end();

    // cleaning up the mess
    for (int i = 0, len = reads.size(); i < len; ++i) {
        delete[] reads[i].sequence;
    }

    delete m;
    delete pool;

    fclose(OUTPUT_FD);

    return 0;
}
