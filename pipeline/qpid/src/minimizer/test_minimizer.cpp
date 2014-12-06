#include "./minimizer.h"
#include <zlib.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "../parser/parser.h"
using std::vector;

typedef HashList<nstring_t, std::pair<unsigned int, unsigned int>> minimizers_t;

// init fasta/fastq reader
KSEQ_INIT(gzFile, gzread)

void read_file(vector<const char *>* string_list, const char *filename) {
    gzFile fp = gzopen(filename, "r");      // STEP 2: open the file handler
    kseq_t *seq = kseq_init(fp);            // STEP 3: initialize seq

    int len = 0;
    while ((len = kseq_read(seq)) > 0) {    // STEP 4: read sequence
        char *read_string = (char *) malloc(len * sizeof(char));
        memcpy(read_string, seq->seq.s, len);
        string_list->push_back(read_string);
    }

    kseq_destroy(seq);      // STEP 5: destroy seq
    gzclose(fp);            // STEP 6: close the file handler
}

int main(int argc, char** argv) {

    if (argc < 2) {
        printf("Usage: %s filename\n", argv[0]);
        exit(1);
    }

    vector<const char *> string_list;
    read_file(&string_list, argv[1]);

    for (int i = 0, len = string_list.size(); i < len; ++i) {
        Minimizer m(16, 20);
        m.calculate_and_store(i, string_list[i]);

        const minimizers_t& minimizers = m.get_minimizers();
        const std::unordered_set<nstring_t>& keys = minimizers.keys();

        vector<int> mins;
        for (auto iter = keys.begin(), eend = keys.end(); iter != eend; ++iter) {
            auto list = minimizers.get_list(*iter);

            for (auto min1 = list.begin(), end = list.end(); min1 != end; ++min1) {
                int pos = (*min1).second;
                mins.push_back(pos);
            }

        }

        std::sort(mins.begin(), mins.end());
        for (int i = 0, len =  mins.size(); i < len; ++i) {
            printf("%d ", mins[i]);
        }

        printf("\n");
    }

    return 0;
}
