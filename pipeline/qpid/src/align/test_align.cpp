#include "./align.h"
#include <cstring>
#include <cassert>
#include <cstdio>
#include <utility>

template <class T>
void print_matrix(T** matrix, int r, int c, const char *format = "%3d ") {

    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            printf(format, (int) matrix[i][j]);
        }
        printf("\n");
    }
}

void test1() {

    char a[] = "ACGTACGT";
    char b[] = "ACGTACGT";
    std::pair<int, int> start, end;

    int score = banded_overlap(a, strlen(a), b, strlen(b), -10, 10, &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);

    assert(score == (int) strlen(a));
    assert(0 == start.first);
    assert(0 == start.second);
    assert(8 == end.first);
    assert(8 == end.second);
}

void test2() {

    char a[] = "ACGTAAAA";
    char b[] = "AAAAACGT";
    std::pair<int, int> start, end;

    int score = banded_overlap(a, strlen(a), b, strlen(b), -10, 10, &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);

    assert(score == 4);
    assert(0 == start.first);
    assert(4 == start.second);
    assert(4 == end.first);
    assert(8 == end.second);
}

void test3() {

    char a[] = "ACGTACGTACGTAAAA";
    char b[] = "AAAAACGTAACGTACGT";
    std::pair<int, int> start, end;

    int score = banded_overlap(a, strlen(a), b, strlen(b), -strlen(a), strlen(a), &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);

    // with current scoring, we have two overlaps that worth the same
    assert(score == 12 + INDEL_SCORE + GAP_SCORE);
    assert(0 == start.first);
    assert(4 == start.second || 9 == start.second);
    assert(8 == end.first || 13 == end.first);
    assert((int) strlen(b) == end.second);
}

void test4() {

    char a[] = "ACGTACGTT";
    char b[] = "ACGTACCTT";
    std::pair<int, int> start, end;

    int score = banded_overlap(a, strlen(a), b, strlen(b), -strlen(a), strlen(a), &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);

    assert(score == (int) (strlen(a) - 1 + MISMATCH_SCORE));
    assert(0 == start.first);
    assert(0 == start.second);
    assert((int) strlen(a) == end.first);
    assert((int) strlen(b) == end.second);
}

void test5() {

    char a[] = "ACGTACGTACGTACGTACGTAAAA";
    char b[] = "AAAAACGTACGTACGTAAACGTACGT";
    std::pair<int, int> start, end;

    // skip first four characters of second string + allow two inserts in first string
    int score = banded_overlap(a, strlen(a), b, strlen(b), 4, 6, &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);

    assert(score == 20 + 2 * INDEL_SCORE + GAP_SCORE);
    assert(0 == start.first);
    assert(4 == start.second);
    assert((int) strlen(a) - 4 == end.first);
    assert((int) strlen(b) == end.second);
}

void test6() {

    char a[] = "AAAAACGTACGTACGTAAACGTACGT";
    char b[] = "ACGTACGTACGTACGTACGTAAAA";
    std::pair<int, int> start, end;

    // skip first four characters of first string + allow two inserts in first string
    int score = banded_overlap(a, strlen(a), b, strlen(b), -6, -4, &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);

    assert(score == 20 + 2 * INDEL_SCORE + GAP_SCORE);
    assert(4 == start.first);
    assert(0 == start.second);
    assert((int) strlen(a) == end.first);
    assert((int) strlen(b) - 4 == end.second);
}

void test7() {

    char a[] = "AGTGTGGCGTATTGGGGGTATGGTACGAAAATTGCTCGGAATATCTACGAGGTCTTTAAAAGTTCGCCGACCTAGTACATCCCAGCCAAAAACCCTGATACAATATATTTCGGGGAGAATACTCA";
    char b[] = "GACCTAGTACATCCCAGCCAAAAACCCTGATACAATATATTTCGGGGAGATTACGCTAGATCAAATAACAAGCTCCCCGCCGCCTGGAATCACAGATCAATAGGCAAGACGACATGAAACCGAAG";

    assert(banded_overlap(a, strlen(a), b, strlen(b), -500, 500)
            == banded_overlap(b, strlen(b), a, strlen(a), -500, 500));
}

void test8() {

    char a[] = "ACGTACGTACGTAAACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    char b[] = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACGTACGTACGTAAACGTACGT";
    std::pair<int, int> start, end;

    int score = banded_overlap(a, strlen(a), b, strlen(b), -500, 500, &start, &end);

    printf("%s %s %d, s: %d %d, e: %d %d\n", a, b, score, start.first, start.second, end.first, end.second);
}

int main() {

    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();

    return 0;
}
