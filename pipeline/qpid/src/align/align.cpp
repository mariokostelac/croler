#include "./align.h"
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <climits>
#define NEG_INF (INT_MIN + 100)

char match_score[4][4] = {
    { MATCH_SCORE,    MISMATCH_SCORE, MISMATCH_SCORE, MISMATCH_SCORE },
    { MISMATCH_SCORE, MATCH_SCORE,    MISMATCH_SCORE, MISMATCH_SCORE },
    { MISMATCH_SCORE, MISMATCH_SCORE, MATCH_SCORE,    MISMATCH_SCORE },
    { MISMATCH_SCORE, MISMATCH_SCORE, MISMATCH_SCORE, MATCH_SCORE }
};

int translate_letter(const char c) {

    if (c == 'A')   return 0;
    if (c == 'C')   return 1;
    if (c == 'G')   return 2;
    if (c == 'T')   return 3;
    return -1;
}

void set_pair(std::pair<int, int>* pair, int first, int second) {

    pair->first = first;
    pair->second = second;
}

int max(int a, int b, int c, int* index) {

    if (a >= b) {
        if (a >= c) {
            *index = 0;
            return a;
        } else {
            *index = 2;
            return c;
        }
    } else {
        if (b >= c) {
            *index = 1;
            return b;
        } else {
            *index = 2;
            return c;
        }
    }
}

int max(int a, int b, int* index) {

    if (a >= b) {
        *index = 0;
        return a;
    } else {
        *index = 1;
        return b;
    }
}

overlap_band_t::overlap_band_t(int height, int width) : width(width), height(height), best_score(NEG_INF) {

    gapa_curr = new int[width];
    middle_prev = new int[width];
    middle_curr = new int[width];
    gapb_curr = new int[width];

    gapa_start = new std::pair<int, int>[width];
    middle_prev_start = new std::pair<int, int>[width];
    middle_curr_start = new std::pair<int, int>[width];
    gapb_start = new std::pair<int, int>[width];
}

overlap_band_t::~overlap_band_t() {

    delete[] gapa_curr;
    delete[] middle_prev;
    delete[] middle_curr;
    delete[] gapb_curr;

    delete[] gapa_start;
    delete[] middle_prev_start;
    delete[] middle_curr_start;
    delete[] gapb_start;
}

void overlap_band_t::initialize(int start_row, int height) {

    this->height = height;

    int init_value = start_row == 1 ? 0 : NEG_INF;

    // forbids using upper row for calculating the current one,
    // or allows skipping arbitrary number of characters if it is the first row
    for (int i = 0; i < width; ++i) {

        gapb_curr[i] = init_value;
        set_pair(&gapb_start[i], start_row - 1, i);

        middle_prev[i] = init_value;
        set_pair(&middle_prev_start[i], start_row - 1, i);

        gapa_curr[i] = init_value;
        set_pair(&gapa_start[i], start_row - 1, i);
    }

    best_score = NEG_INF;
    set_pair(&best_start, -1, -1);
    set_pair(&best_end, -1, -1);
}

void overlap_band_t::update_band(int lo, int hi, int row, const char *a, const char *b) {

    int ca, cb, max_index;

    // allows skipping arbitrary number of characters in first string
    middle_prev[0] = 0;
    set_pair(&middle_prev_start[0], row - 1, 0);

    // ensures that gapb won't use any field from last step that wasn't involved in band
    gapb_curr[lo - 1] = NEG_INF;
    middle_curr[lo - 1] = NEG_INF;

    // ensures that gapa won't use any field from last step that wasn't involved in band
    gapa_curr[hi] = NEG_INF;
    middle_prev[hi] = NEG_INF;

    for (int col = lo; col <= hi; ++col) {

        ca = translate_letter(a[row - 1]);
        cb = translate_letter(b[col - 1]);

        gapb_curr[col] = max(gapb_curr[col - 1] + INDEL_SCORE, middle_curr[col - 1] + INDEL_SCORE + GAP_SCORE, &max_index);
        if (max_index == 0) gapb_start[col] = gapb_start[col - 1];
        else                gapb_start[col] = middle_curr_start[col - 1];

        gapa_curr[col] = max(gapa_curr[col] + INDEL_SCORE, middle_prev[col] + INDEL_SCORE + GAP_SCORE, &max_index);
        if (max_index == 1) gapa_start[col] = middle_prev_start[col];

        middle_curr[col] = max(middle_prev[col - 1] + match_score[ca][cb], gapa_curr[col], gapb_curr[col], &max_index);
        if (max_index == 0)         middle_curr_start[col] = middle_prev_start[col - 1];
        else if (max_index == 1)    middle_curr_start[col] = gapa_start[col];
        else                        middle_curr_start[col] = gapb_start[col];
    }

    if (hi == width - 1) {
        if (middle_curr[hi] > best_score) {
            best_score = middle_curr[hi];
            best_start = middle_curr_start[hi];
            set_pair(&best_end, row, hi);
        }
    }

    if (row == height - 1) {
        for (int col = lo; col <= hi; ++col) {
            if (middle_curr[col] > best_score) {
                best_score = middle_curr[col];
                best_start = middle_curr_start[col];
                set_pair(&best_end, row, col);
            }
        }
    }

    // swaps current and prev buffer
    int *tmp = middle_curr;
    middle_curr = middle_prev;
    middle_prev = tmp;

    std::pair<int, int>* tmp_pair = middle_curr_start;
    middle_curr_start = middle_prev_start;
    middle_prev_start = tmp_pair;
}

char** create_char_matrix(int r, int c) {
    char **matrix = new char*[r];

    for (int i = 0; i < r; ++i) {
        matrix[i] = new char[c];
    }

    return matrix;
}

int** create_int_matrix(int r, int c) {
    int **matrix = new int*[r];

    for (int i = 0; i < r; ++i) {
        matrix[i] = new int[c];
    }

    return matrix;
}

int local_alignment(const char* a, int alen, const char* b, int blen, std::pair<int, int>* start, std::pair<int, int>* end) {

    int** score = create_int_matrix(alen + 1, blen + 1);
    char** backtrack = create_char_matrix(alen + 1, blen + 1);

    // set border penalties
    for (int i = 0; i < alen + 1; ++i) {
        score[i][0] = 0;
    }
    for (int i = 0; i < blen + 1; ++i) {
        score[0][i] = 0;
    }

    // calculate score
    for (int i = 1; i < alen + 1; ++i) {
        for (int j = 1; j < blen + 1; ++j) {
            int ca = translate_letter(a[i - 1]), cb = translate_letter(b[j - 1]);
            int match = score[i-1][j-1] + match_score[ca][cb];
            int gapa = score[i][j-1] - 3;
            int gapb = score[i-1][j] - 3;

            score[i][j] = match;
            backtrack[i][j] = 'M';

            if (gapa > score[i][j]) {
                score[i][j] = gapa;
                backtrack[i][j] = 'A';
            }
            if (gapb > score[i][j]) {
                score[i][j] = gapb;
                backtrack[i][j] = 'B';
            }
            if (0 > score[i][j]) {
                score[i][j] = 0;
                backtrack[i][j] = 'S';
            }
        }
    }

    int max_i = 0, max_j = 0;
    for (int i = 0; i < alen + 1; ++i) {
        for (int j = 0; j < blen + 1; ++j) {
            if (score[i][j] > score[max_i][max_j]) {
                max_i = i;
                max_j = j;
            }
        }
    }

    int best_score = score[max_i][max_j];

    // store end point
    if (end != NULL) {
        end->first = max_i;
        end->second = max_j;
    }

    // backtrack to start
    if (start != NULL) {
        int i = max_i, j = max_j;
        while (i >  0 && j > 0 && backtrack[i][j] != 'S') {
            if (backtrack[i][j] == 'M') {
                --i; --j;
            } else if (backtrack[i][j] == 'A') {
                --j;
            } else if (backtrack[i][j] == 'B') {
                --i;
            }
        }
        start->first = i;
        start->second = j;
    }

    for (int i = 0; i < alen + 1; ++i) {
        delete[] score[i];
        delete[] backtrack[i];
    }
    delete[] score;
    delete[] backtrack;

    return best_score;
}

int overlap_alignment(const char* a, int alen, const char* b, int blen,
        std::pair<int, int>* start, std::pair<int, int>* end, int** score_out, char** backtrack_out) {

    int** score = score_out == NULL ? create_int_matrix(alen + 1, blen + 1) : score_out;
    char** backtrack = backtrack_out == NULL ? create_char_matrix(alen + 1, blen + 1) : backtrack_out;

    // set border penalties
    for (int i = 0; i < alen + 1; ++i) {
        score[i][0] = 0;
    }
    for (int i = 0; i < blen + 1; ++i) {
        score[0][i] = 0;
    }

    // calculate score
    for (int i = 1; i < alen + 1; ++i) {
        for (int j = 1; j < blen + 1; ++j) {
            int ca = translate_letter(a[i - 1]), cb = translate_letter(b[j - 1]);
            int match = score[i-1][j-1] + match_score[ca][cb];
            int gapa = score[i][j-1] - 3;
            int gapb = score[i-1][j] - 3;

            score[i][j] = match;
            backtrack[i][j] = 'M';

            if (gapa > score[i][j]) {
                score[i][j] = gapa;
                backtrack[i][j] = 'A';
            }

            if (gapb > score[i][j]) {
                score[i][j] = gapb;
                backtrack[i][j] = 'B';
            }
        }
    }

    int max_i = alen, max_j = blen;

    // check right edge
    for (int i = 0; i < alen; ++i) {
        if (score[i][blen] > score[max_i][max_j]) {
            max_i = i;
            max_j = blen;
        }
    }

    // check bottom edge
    for (int i = 0; i < blen; ++i) {
        if (score[alen][i] > score[max_i][max_j]) {
            max_i = alen;
            max_j = i;
        }
    }

    int best_score = score[max_i][max_j];

    // store end point
    if (end != NULL) {
        end->first = max_i;
        end->second = max_j;
    }

    // backtrack to start
    if (start != NULL) {
        int i = max_i, j = max_j;
        while (i >  0 && j > 0) {
            if (backtrack[i][j] == 'M') {
                --i; --j;
            } else if (backtrack[i][j] == 'A') {
                --j;
            } else if (backtrack[i][j] == 'B') {
                --i;
            }
        }
        start->first = i;
        start->second = j;
    }

    if (score_out == NULL) {
        for (int i = 0; i < alen + 1; ++i) {
            delete[] score[i];
        }
        delete[] score;
    }

    if (backtrack_out == NULL) {
        for (int i = 0; i < alen + 1; ++i) {
            delete[] backtrack[i];
        }
        delete[] backtrack;
    }

    return best_score;
}

int banded_overlap2(const char* a, int alen, const char* b, int blen, int d_min, int d_max,
        std::pair<int, int>* start, std::pair<int, int>* end, int** score_out, char** backtrack_out) {

    int** score = score_out == NULL ? create_int_matrix(alen + 1, blen + 1) : score_out;
    char** backtrack = backtrack_out == NULL ? create_char_matrix(alen + 1, blen + 1) : backtrack_out;

    // set border penalties
    for (int i = 0; i < alen + 1; ++i) {
        score[i][0] = 0;
    }
    for (int i = 0; i < blen + 1; ++i) {
        score[0][i] = 0;
    }

    // calculate score
    for (int i = 1; i < alen + 1; ++i) {

        // calculate band (lo inclusive, hi inclusive)
        int lo = std::max(d_min + i, 1);
        int hi = std::min(d_max + i, blen);
        int ca, cb, match, gapa, gapb;

        // if second string haven't started yet (both diagonals negative)
        if (lo > hi) continue;

        // if both diagonals are positive, we don't need full height of matrix
        if (lo > blen) break;

        // if lo == hi, we have case where we have to calculate band of width 1
        // so we have to be careful to use just already calculated values
        if (lo == hi) {
            if (lo == 1) {
                // band is on the left edge of matrix, we just use upper left field
                ca = translate_letter(a[i - 1]);
                cb = translate_letter(b[lo - 1]);
                match = score[i-1][lo-1] + match_score[ca][cb];
                score[i][lo] = match;
                backtrack[i][lo] = 'M';
            } else {
                // band is on the right edge of matrix, we just use upper left and upper field
                ca = translate_letter(a[i - 1]);
                cb = translate_letter(b[hi - 1]);
                match = score[i-1][hi-1] + match_score[ca][cb];
                gapb = score[i-1][hi] + MISMATCH_SCORE;

                if (match > gapb) {
                    score[i][hi] = match;
                    backtrack[i][hi] = 'M';
                } else {
                    score[i][hi] = gapb;
                    backtrack[i][hi] = 'B';
                }
            }
            continue;
        }

        // calculate first field in the band
        ca = translate_letter(a[i - 1]);
        cb = translate_letter(b[lo - 1]);
        match = score[i-1][lo-1] + match_score[ca][cb];
        gapb = score[i-1][lo] + MISMATCH_SCORE;

        if (match > gapb) {
            score[i][lo] = match;
            backtrack[i][lo] = 'M';
        } else {
            score[i][lo] = gapb;
            backtrack[i][lo] = 'B';
        }

        // first and last in band are calculated on different way
        for (int j = lo + 1; j < hi; ++j) {
            ca = translate_letter(a[i - 1]), cb = translate_letter(b[j - 1]);
            match = score[i-1][j-1] + match_score[ca][cb];
            gapa = score[i][j-1] + MISMATCH_SCORE;
            gapb = score[i-1][j] + MISMATCH_SCORE;

            score[i][j] = match;
            backtrack[i][j] = 'M';

            if (gapa > score[i][j]) {
                score[i][j] = gapa;
                backtrack[i][j] = 'A';
            }

            if (gapb > score[i][j]) {
                score[i][j] = gapb;
                backtrack[i][j] = 'B';
            }
        }

        // calculate last field in the band
        ca = translate_letter(a[i - 1]);
        cb = translate_letter(b[hi - 1]);
        match = score[i-1][hi-1] + match_score[ca][cb];
        gapa = score[i][hi-1] + MISMATCH_SCORE;

        if (match > gapa) {
            score[i][hi] = match;
            backtrack[i][hi] = 'M';
        } else {
            score[i][hi] = gapa;
            backtrack[i][hi] = 'B';
        }
    }

    int max_i = -1, max_j = -1;

    // check bottom edge
    int b_start = d_min + alen, b_end = std::min(d_max + alen, blen);
    if (b_start <= blen) {
        max_i = alen;
        max_j = b_start;
        for (int j = b_start; j <= b_end; ++j) {
            if (score[alen][j] > score[alen][max_j]) {
                max_j = j;
            }
        }
    }

    // check right edge
    for (int i = 0; i <= alen; ++i) {
        int d = blen - i;
        if (d < d_min || d > d_max) continue;
        if (max_i == -1 || score[i][blen] > score[max_i][max_j]) {
            max_i = i;
            max_j = blen;
        }
    }

    int best_score = score[max_i][max_j];

    // store end point
    if (end != NULL) {
        end->first = max_i;
        end->second = max_j;
    }

    // backtrack to start
    if (start != NULL) {
        int i = max_i, j = max_j;
        while (j - i >= d_min && j - i <= d_max && i > 0 && j > 0) {
            if (backtrack[i][j] == 'M') {
                --i; --j;
            } else if (backtrack[i][j] == 'A') {
                --j;
            } else if (backtrack[i][j] == 'B') {
                --i;
            }
        }
        start->first = i;
        start->second = j;
    }

    if (score_out == NULL) {
        for (int i = 0; i < alen + 1; ++i) {
            delete[] score[i];
        }
        delete[] score;
    }

    if (backtrack_out == NULL) {
        for (int i = 0; i < alen + 1; ++i) {
            delete[] backtrack[i];
        }
        delete[] backtrack;
    }

    return best_score;
}


int banded_overlap(const char* a, int alen, const char* b, int blen, int d_min, int d_max,
        std::pair<int, int>* start, std::pair<int, int>* end, overlap_band_t* band_out) {

    overlap_band_t* band = band_out;

    if (band == NULL) band = new overlap_band_t(alen + 1, blen + 1);

    int start_row = 1;
    if (d_max < 0) start_row = -d_max + 1;

    band->initialize(start_row, alen + 1);

    // calculate score
    for (int i = start_row; i < alen + 1; ++i) {

        // calculate band (lo inclusive, hi inclusive)
        int lo = std::max(d_min + i, 1);
        int hi = std::min(d_max + i, blen);

        // if both diagonals are positive, we don't need full height of matrix
        if (lo > blen) break;

        band->update_band(lo, hi, i, a, b);
    }

    if (start != NULL)  *start = band->best_start;
    if (end != NULL)    *end = band->best_end;

    int best_score = band->best_score;

    if (band_out == NULL) delete band;

    return best_score;
}
