// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_LAYOUT_UTILS_H_
#define LAYOUT_LAYOUT_UTILS_H_

#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/contig.h>
#include <layout/unitigging.h>

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
 * Reads all overlaps form the .afg file.
 */
overlap::OverlapSet* ReadOverlapsAfg(overlap::ReadSet* read_set, FILE *fd);

/**
 * Finds the n50 of the given contig set.
 * http://en.wikipedia.org/wiki/N50_statistic
 */
int n50(Unitigging::ContigSetPtr contig_set);
};  // namespace layout

#endif  // LAYOUT_LAYOUT_UTILS_H_
