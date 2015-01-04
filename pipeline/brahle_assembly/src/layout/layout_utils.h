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
   * Reads all reads from the .afg file.
   */
  uint32_t ReadReadsAfg(overlap::ReadSet& container, const char *filename);

  /**
   * Reads all overlaps from the .afg file.
   */
  overlap::OverlapSet* ReadOverlapsAfg(overlap::ReadSet* read_set, FILE *fd);

  /**
   * Finds the n50 of the given contig set.
   * http://en.wikipedia.org/wiki/N50_statistic
   */
  int n50(Unitigging::ContigSetPtr contig_set);

  /**
   * Calculates lengths of overlaps between read_one and read_two.
   */
  std::pair<int, int> getOverlapLengths(const overlap::ReadSet* read_set, const int read_one, const int read_two, const int hang_one, const int hang_two);

  int ContigsToFile(std::shared_ptr<layout::ContigSet> contigs, const char *contigs_filename);

  std::string dot_graph(overlap::ReadSet* reads, overlap::OverlapSet* overlaps);

  std::string dot_graph(const BetterReadSet* reads, const BetterOverlapSet* overlaps);
};  // namespace layout

#endif  // LAYOUT_LAYOUT_UTILS_H_
