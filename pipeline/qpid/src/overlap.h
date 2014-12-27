
#ifndef _OVERLAP_H
#define _OVERLAP_H

#include <cstring>
#include <cassert>
#include <vector>
#include <utility>
#include <algorithm>
#include "minimizer/minimizer.h"
#include "../lib/ThreadPool/ThreadPool.h"
#include "align/align.h"
#include "read.h"
using std::pair;
using std::swap;
using std::vector;

struct offset_t {
  int index;
  int lo_offset;
  int hi_offset;
  offset_t() {}
  offset_t(unsigned int index, int offset) : index(index), lo_offset(offset), hi_offset(offset) {}
};

typedef HashList<nstring_t, std::pair<unsigned int, unsigned int>> minimizers_t;

typedef struct Overlap {
  Read r1;
  Read r2;
  double score;
  int a_hang;
  int b_hang;
  int errors;
  double error_rate;
  bool normal_overlap;

  Overlap() {}

  // from http://sourceforge.net/p/amos/mailman/message/19965222/.
  //
  // read a      -------------------|-------------->     bhang
  // read b            ahang     ---------------|--------------->
  //
  // read a           -ahang     ---------------|--------------->
  // read b      -------------------|-------------->     -bhang
  Overlap(Read& r1, Read& r2, const double& score, const pair<int, int>& start, const pair<int, int>& end,
    const bool& first_forward, const bool& second_forward)
    : r1(r1), r2(r2), score(score), normal_overlap(first_forward) {

      assert(second_forward == true);

      int len1 = strlen(r1.sequence);
      int len2 = strlen(r2.sequence);
      int overlap_len_a = end.first - start.first;
      int overlap_len_b = end.second - start.second;

      int a_not_matching = len1 - overlap_len_a;
      int b_not_matching = len2 - overlap_len_b;

      if (start.first == 0 && end.first == len1) {
        // first contained
        a_hang = -start.second;
        b_hang = len2 - end.second;
      } else if (start.second == 0 && end.second == len2) {
        // second contained
        a_hang = start.first;
        b_hang = -(len1 - end.first);
      } else if (end.first == len1) {
        // first case from the comment
        a_hang = a_not_matching;
        b_hang = b_not_matching;
      } else if (end.second == len2) {
        // second case from the comment
        a_hang = -b_not_matching;
        b_hang = -a_not_matching;
      } else {
        assert(false);
      }

      // if all the values are calculated for reversed_complement(r1), r2, we store overlap for
      // r1, reverse_complement(r2)
      if (!first_forward) {
        swap(a_hang, b_hang);
        a_hang *= -1;
        b_hang *= -1;
      }

      int len = abs(end.first - start.first + end.second - start.second) / 2.;
      errors = (score - len)/(INDEL_SCORE + GAP_SCORE + MISMATCH_SCORE);
      error_rate = (double) errors / len;
    }
} Overlap;

#endif
