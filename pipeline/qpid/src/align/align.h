#include <cstdlib>
#include <utility>

const int MATCH_SCORE = 1;
const int MISMATCH_SCORE = -3;
const int INDEL_SCORE = -3;
const int GAP_SCORE = -1;

char** create_char_matrix(int r, int c);
int** create_int_matrix(int r, int c);

struct overlap_band_t {

    int width;
    int height;

    int *gapa_curr;
    int *middle_prev;
    int *middle_curr;
    int *gapb_curr;

    std::pair<int, int>* gapa_start;
    std::pair<int, int>* middle_prev_start;
    std::pair<int, int>* middle_curr_start;
    std::pair<int, int>* gapb_start;

    int best_score;
    std::pair<int, int> best_start;
    std::pair<int, int> best_end;

    overlap_band_t(int height, int width);
    ~overlap_band_t();
    void initialize(int start_row, int height);
    void update_band(int lo, int hi, int row, const char *a, const char *b);
};

int local_alignment(const char* a, int alen, const char* b, int blen, std::pair<int, int>* start, std::pair<int, int>* end);

int overlap_alignment(const char* a, int alen, const char* b, int blen,
        std::pair<int, int>* start, std::pair<int, int>* end, int** score = NULL, char** backtrack = NULL);

int banded_overlap2(const char* a, int alen, const char* b, int blen, int d_min, int d_max,
        std::pair<int, int>* start, std::pair<int, int>* end, int** score = NULL, char** backtrack = NULL);

int banded_overlap(const char* a, int alen, const char* b, int blen, int d_min, int d_max,
        std::pair<int, int>* start = NULL, std::pair<int, int>* end = NULL);
