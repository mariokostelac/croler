// Copyright 2014 Bruno Rahle

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <utility>
#include <zlib.h>
#include <string>
#include <sstream>
#include <unordered_map>

#include "lib/amos/reader.cpp"

#include "layout/layout_utils.h"

using std::stringstream;
using std::string;
using std::swap;
using std::unordered_map;
using std::vector;

namespace layout {

  /**
   * Since ids in ReadSet start with 0 and real ids (reads read from afg file/AMOS bank) can start with arbitrary number,
   * we have to map real_id -> internal_id (sequence that starts with 0).
   * That's why introduced this type.
   */
  typedef std::unordered_map<int, int> ReadIdMap;

  ReadIdMap _MapIds(overlap::ReadSet* reads) {
    ReadIdMap mapped;

    int reads_len = reads->size();
    for (int i = 0; i < reads_len; ++i) {
      auto& read = (*reads)[i];
      if (mapped.count(read->orig_id())) {
        fprintf(stderr, "Read with orig_id '%d' already seen\n", read->orig_id());
        exit(2);
      }
      mapped[read->orig_id()] = read->id();
    }

    return mapped;
  }

  uint32_t ReadReadsAfg(overlap::ReadSet& container, const char* filename) {
    int records = 0;

    clock_t start = clock();
    auto read_set = new overlap::ReadSet(10000);
    vector< AMOS::Read* > tmp_reads;
    uint32_t reads_size = AMOS::get_reads(tmp_reads, filename);
    if (reads_size < 0) {
      return reads_size;
    }

    for (int i = 0; i < reads_size; ++i) {
      const auto& read = tmp_reads[i];
      container.Add(new overlap::Read(
            (uint8_t*) read->seq,
            read->clr_lo,
            read->clr_hi,
            i,
            read->iid
      ));
    }

    printf(
        "Reads read in %.2lfs\n",
        (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));

    return records;
  }

  overlap::OverlapSet* ReadOverlapsAfg(overlap::ReadSet* read_set, FILE *fd) {
    clock_t start = clock();

    ReadIdMap internal_id = _MapIds(read_set);

    overlap::OverlapSet* overlap_set = new overlap::OverlapSet(10000);
    char type;
    int read_one;
    int read_two;
    int score;
    int hang_one;
    int hang_two;
    while (fscanf(
          fd,
          " {OVL adj:%c rds:%d,%d scr:%d ahg:%d bhg:%d }",
          &type,
          &read_one,
          &read_two,
          &score,
          &hang_one,
          &hang_two) == 6) {

      if (!internal_id.count(read_one)) {
        fprintf(stderr, "Read with orig_id '%d' has not been found\n", read_one);
        exit(3);
      }

      if (!internal_id.count(read_two)) {
        fprintf(stderr, "Read with orig_id '%d' has not been found\n", read_two);
        exit(3);
      }

      if (type != 'N' && type != 'I') {
        fprintf(stderr, "Unkown overlap type '%c'\n", type);
        exit(3);
      }

      read_one = internal_id[read_one];
      read_two = internal_id[read_two];

      std::pair<int, int> lenghts = getOverlapLengths(read_set, read_one, read_two, hang_one, hang_two);
      overlap_set->Add(new overlap::Overlap(
            read_one,
            read_two,
            lenghts.first,
            lenghts.second,
            hang_one,
            hang_two,
            type == 'N' ? overlap::Overlap::Type::EB : overlap::Overlap::Type::EE,
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

  /**
   * Calculates lengths of overlaps between read_one and read_two.
   */
  std::pair<int, int> getOverlapLengths(const overlap::ReadSet* read_set, const int a, const int b, const int a_hang, const int b_hang) {
    int len_one, len_two;
    int a_size = (*read_set)[a]->size();
    int b_size = (*read_set)[b]->size();
    //
    // -------|--------------> b_hang
    // a_hang ---------------|------>
    //
    // -a_hang -------------|------->
    // --------|------------> -b_hang
    //
    // -------|-------------|------->
    // a_hang --------------> -b_hang
    //
    // -a_hang --------------> b_hang
    // --------|-------------|------>
    //
    if (a_hang >= 0 && b_hang >= 0) {
      len_one = a_size - a_hang;
      len_two = b_size - b_hang;
    } else if (a_hang <= 0 && b_hang <= 0) {
      len_one = a_size + b_hang;
      len_two = b_size + a_hang;
    } else if (a_hang >= 0 && b_hang <= 0) {
      len_one = a_size + b_hang - a_hang;
      len_two = b_size;
    } else if (a_hang <= 0 && b_hang >= 0) {
      len_one = a_size;
      len_two = b_size + a_hang - b_hang;
    } else {
      // case not covered
      assert(false);
    }

    assert(len_one > 0 && len_two > 0);
    return std::make_pair(len_one, len_two);
  }

  /**
   * Writes all usable contigs from contig set to a file with given filename.
   * If file does not exist, it will be created.
   * Returns the number of written contigs.
   */
  int ContigsToFile(std::shared_ptr<ContigSet> contigs, const char *contigs_filename) {
    FILE *contigs_file = fopen(contigs_filename, "w");
    if (contigs_file == nullptr) {
      return -1;
    }

    int written = 0;
    int contigs_size = contigs->size();
    for (int i = 0; i < contigs_size; ++i) {
      // skip non-usable contigs
      if (!((*contigs)[i]->IsUsable())) continue;


      const std::deque<BetterRead*> &reads = (*contigs)[i]->getReads();
      const std::deque<BetterOverlap*> &overlaps = (*contigs)[i]->getOverlaps();

      bool forward = true;
      uint32_t offset = 0;

      const auto& f_read = reads[0];
      const auto& f_overlap = overlaps[0]->overlap();

      if (f_read->id() == f_overlap->read_two && f_overlap->type == overlap::Overlap::Type::EE) {
        forward = false;
      }

      auto process_read =
        [&contigs_file, &forward, &offset] (const BetterRead* r, const BetterOverlap* o) {
          const auto& read = r->read();
          const auto& overlap = o->overlap();
          bool eb = overlap->type == overlap::Overlap::Type::EB;
          uint32_t lo = read->lo();
          uint32_t hi = read->hi();
          if (!forward) {
            swap(lo, hi);
          }

          fprintf(contigs_file, "{TLE\n");
          fprintf(contigs_file, "clr:%u,%u\n", lo, hi);
          fprintf(contigs_file, "off:%u\n", offset);
          fprintf(contigs_file, "src:%d\n}\n", read->orig_id());

          if (read->id() == overlap->read_one) {
            if (overlap->a_hang > 0) {
              offset += overlap->a_hang;
            } else {
              offset += abs(overlap->b_hang);
            }
          } else if (read->id() == overlap->read_two) {
            if (overlap->b_hang > 0) {
              offset += overlap->b_hang;
            } else {
              offset += abs(overlap->a_hang);
            }
          }

          if (overlap->type == overlap::Overlap::Type::EE) {
            forward = !forward;
          }
        };

      fprintf(contigs_file, "{LAY\n");

      process_read(reads[0], overlaps[0]);
      int num_reads = reads.size();
      for (int j = 1; j < num_reads - 1; ++j) {
        process_read(reads[j], overlaps[j]);
      }
      process_read(reads[num_reads-1], overlaps[num_reads-2]);
      fprintf(contigs_file, "}\n");

      written++;
    }

    fclose(contigs_file);
    return written;
  }

  std::string dot_graph(overlap::ReadSet* reads, overlap::OverlapSet* overlaps) {
    BetterReadSet brs(reads, false);
    BetterOverlapSet bos(reads, overlaps);
    return dot_graph(&brs, &bos);
  }

  std::string dot_graph(const BetterReadSet* reads, const BetterOverlapSet* overlaps) {
    stringstream graph;
    graph << "digraph overlaps {\n";
    int overlaps_size = overlaps->size();
    for (int i = 0; i < overlaps_size; ++i) {
      const auto& overlap = (*overlaps)[i];
      int read1 = (*reads)[overlap->overlap()->read_one]->read()->orig_id();
      int read2 = (*reads)[overlap->overlap()->read_two]->read()->orig_id();
      if (overlap->GoesFrom(overlap->overlap()->read_one)) {
        graph << read1 << " -> " << read2;
      } else {
        graph << read2 << " -> " << read1;
      }
      if (overlap->overlap()->type == overlap::Overlap::Type::EB) {
        graph << " [color=green] ";
      } else {
        graph << " [color=pink] ";
      }
      graph << ";\n";
    }
    graph << "}\n";
    return graph.str();
  }
};  // namespace layout
