// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_LAYOUT_UTILS_H_
#define LAYOUT_LAYOUT_UTILS_H_

#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/contig.h>
#include <layout/unitigging.h>

#include <memory>
#include <unordered_map>

namespace layout {

  /**
   * Reads the next read from an .afg file. File pointer must point to the start of the read.
   */
  overlap::Read* ReadOneReadAfg(FILE *fd);

  /**
   * Reads all reads from the .afg file.
   */
  overlap::ReadSet* ReadReadsAfg(FILE *fd);

  /**
   * Reads all reads from the .seq (fastq format) file.
   */
  overlap::ReadSet* ReadReadsSeq(char *filename);

  /**
   * Reads all reads from AMOS bank.
   */
  overlap::ReadSet* ReadReadsAmos(const char *bank_name);

  /**
   * Reads all overlaps form the .afg file.
   */
  overlap::OverlapSet* ReadOverlapsAfg(overlap::ReadSet* read_set, FILE *fd);

  /**
   * Reads all overlaps form the AMOS bank.
   */
  overlap::OverlapSet* ReadOverlapsAmos(overlap::ReadSet* read_set, const char* bank_name);

  /**
   * Finds the n50 of the given contig set.
   * http://en.wikipedia.org/wiki/N50_statistic
   */
  int n50(Unitigging::ContigSetPtr contig_set);

  /**
   * Calculates lengths of overlaps between read_one and read_two.
   */
  std::pair<int, int> getOverlapLengths(const overlap::ReadSet* read_set, const int read_one, const int read_two, const int hang_one, const int hang_two);

};  // namespace layout

#endif  // LAYOUT_LAYOUT_UTILS_H_
