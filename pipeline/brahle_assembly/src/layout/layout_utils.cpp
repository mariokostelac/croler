// Copyright 2014 Bruno Rahle

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>

#include "layout/layout_utils.h"

namespace layout {

/**
 * Reads data for one read from the .afg file.
 */
uint8_t* _ReadDataForOneRead(FILE *fd) {
  static char buff[1 << 20];
  char *buff_free = buff;
  fscanf(fd, " seq: ");
  buff[0] = 0;
  do {
    // We remove one character for '\n'. However, in the first step, strlen is
    // 0, so we don't want to move in the negative direction.
    buff_free += std::max(0, static_cast<int>(strlen(buff_free) - 1));
    fgets(buff_free, sizeof(buff) - sizeof(char) * (buff_free - buff), fd);
  } while (buff_free[0] != '.');
  int length = buff_free - buff;
  char* data = new char[length+1];
  strncpy(data, buff, length);
  data[length] = 0;
  return reinterpret_cast<uint8_t*>(data);
}

overlap::Read* ReadOneReadAfg(FILE *fd) {
  static char buff[1 << 20];
  char *buff_free = buff;
  int id;
  int length;
  fscanf(fd, " iid:%d eid:%*s", &id);
  uint8_t *data = _ReadDataForOneRead(fd);
  do {
    fgets(buff, sizeof(buff), fd);
  } while (buff[0] != '.');
  fscanf(fd, " frg:%*d clr:%*d,%d", &length);
  return new overlap::Read(data, length, (id-1), id);
}

overlap::ReadSet* ReadReadsAfg(FILE *fd) {
  clock_t start = clock();
  auto read_set = new overlap::ReadSet(10000);
  int i = 0;
  char buff[1 << 20];
  while (fscanf(fd, " %s", buff) == 1) {
    if (!strcmp(buff, "{RED")) {
      read_set->Add(ReadOneReadAfg(fd));
    }
  }
  printf(
      "Reads read in %.2lfs\n",
      (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));
  return read_set;
}

overlap::OverlapSet* ReadOverlapsAfg(overlap::ReadSet* read_set, FILE *fd) {
  clock_t start = clock();
  overlap::OverlapSet* overlap_set = new overlap::OverlapSet(10000);
  char type;
  int read_one;
  int read_two;
  int score;
  int hang_one;
  int hang_two;
  int len_one;
  int len_two;
  while (fscanf(
             fd,
             " {OVL adj:%c rds:%d,%d scr:%d ahg:%d bhg:%d }",
             &type,
             &read_one,
             &read_two,
             &score,
             &hang_one,
             &hang_two) == 6) {
    --read_one;
    --read_two;
    if (hang_one >= 0) {
      len_one = (*read_set)[read_one]->size() - abs(hang_one);
    } else {
      len_one = std::min(
          (*read_set)[read_one]->size(),
          (*read_set)[read_two]->size() - abs(hang_one));
    }
    if (hang_two >= 0) {
      len_two = (*read_set)[read_two]->size() - abs(hang_two);
    } else {
      len_two = std::min(
          (*read_set)[read_two]->size(),
          (*read_set)[read_one]->size() - abs(hang_two));
    }
    overlap_set->Add(new overlap::Overlap(
        read_one,
        read_two,
        len_one,
        len_two,
        overlap::Overlap::Type::EB,
        0));
  }
  printf(
      "Overlaps read in %.2lfs\n",
      (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));
  return overlap_set;
}

int n50(Unitigging::ContigSetPtr contig_set) {
  std::vector< int > v;
  int cnt = 0;
  v.reserve(contig_set->size()*10);
  for (size_t i = 0; i < contig_set->size(); ++i) {
    if ((*contig_set)[i]->IsUsable()) {
      printf("Contig %d: size = %d\n", ++cnt, (*contig_set)[i]->size());
      for (size_t j = 0; j < (*contig_set)[i]->size(); ++j) {
        v.push_back((*contig_set)[i]->size());
      }
    }
  }
  if (v.size() == 0) {
    return 0;
  }
  std::sort(v.begin(), v.end());
  return v[v.size()/2];
}

};  // namespace layout

