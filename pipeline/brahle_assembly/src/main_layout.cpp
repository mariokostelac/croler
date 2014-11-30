// Copyright 2014 Bruno Rahle
#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/unitigging.h>
#include <layout/layout_utils.h>
#include <layout/string_graph.h>
#include <layout/contig.h>
#include <layout/better_read.h>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <deque>
#include <vector>
#include <map>
#include <utility>

void usage(char *argv[]) {
  fprintf(
      stderr,
      "Usage: %s <minimus_folder> | (<reads_file.seq> <overlaps_file.afg>)\n",
      argv[0]);
  exit(0);
}

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    usage(argv);
  }
  char reads_file_name[1024];
  char overlaps_file_name[1024];
  if (argc == 2) {
    snprintf(
        reads_file_name,
        sizeof(reads_file_name),
        "%s/minimus_results/reads.afg",
        argv[1]);
    snprintf(
        overlaps_file_name,
        sizeof(overlaps_file_name),
        "%s/minimus_results/overlaps.afg",
        argv[1]);
  } else {
    snprintf(reads_file_name, sizeof(reads_file_name), "%s", argv[1]);
    snprintf(overlaps_file_name, sizeof(overlaps_file_name), "%s", argv[2]);
  }

  FILE *reads_file = fopen(reads_file_name, "r");
  if (reads_file == nullptr) {
    fprintf(
        stderr,
        "ERROR: reads file ('%s') cannot be opened!\n",
        reads_file_name);
    exit(1);
  }
  fclose(reads_file);
  FILE *overlaps_file = fopen(overlaps_file_name, "r");
  if (overlaps_file == nullptr) {
    fprintf(
        stderr,
        "ERROR: overlaps file ('%s') cannot be opened!\n",
        overlaps_file_name);
    exit(1);
  }
  FILE *graphviz_file = fopen("graph.dot", "w");
  if (reads_file == nullptr) {
    fprintf(
        stderr,
        "ERROR: graphviz file ('%s') cannot be opened!\n",
        "graph.dot");
    exit(1);
  }

  std::shared_ptr< overlap::ReadSet > reads(layout::ReadReadsSeq(reads_file_name));
  fprintf(stderr, "Number of reads = %u\n", reads->size());

  std::shared_ptr< overlap::OverlapSet > overlaps(
      layout::ReadOverlapsAfg(reads.get(), overlaps_file));
  fprintf(stderr, "Number of overlaps = %u\n", overlaps->size());

  // create map for overlaps
  // @mculinovic
  std::map< std::pair<uint32_t, uint32_t>, uint32_t> overlaps_map; 
  size_t overlaps_size = (*overlaps).size();
  for (size_t i = 0; i < overlaps_size; ++i) {
      overlap::Overlap* overlap = (*overlaps)[i];
      std::pair <uint32_t, uint32_t> key = std::make_pair(overlap->read_one, overlap->read_two);
      overlaps_map[key] = (overlap->len_one + overlap->len_two) / 2;
  }
  // end creating map  

  clock_t start = clock();
  std::shared_ptr< layout::Unitigging > u(
      new layout::Unitigging(reads.get(), overlaps.get()));
  u->start();
  fprintf(
      stderr,
      "Unitigging finished in %.2lfs\n",
      (clock() - start) / static_cast<double>(CLOCKS_PER_SEC));

  int n50_value = layout::n50(u->contigs());
  fprintf(stderr, "n50 = %d\n", n50_value);

  start = clock();
  layout::Graph g = layout::Graph::create(u->readSet(), u->noTransitives());
  fprintf(
      stderr,
      "String graph constructed in %.2lfs\n",
      (clock() - start) / static_cast<double>(CLOCKS_PER_SEC));

  g.printToGraphviz(graphviz_file);


  // output afg
  // @mculinovic
  // TODO(mculinovic) ugly output
  FILE *afg_file = fopen("layout.afg", "w");

  std::shared_ptr<layout::ContigSet> contigs = u->contigs();
  int contigs_size = contigs->size();
  for (int i = 0; i < contigs_size; ++i) {
      // skip non-usable contigs
      if (!((*contigs)[i]->IsUsable())) continue;

      fprintf(afg_file, "{LAY\n");
      uint32_t offset = 0;
      const std::deque< layout::BetterRead* > &reads = (*contigs)[i]->getReads();
      if (reads.size() == 1) {
          fprintf(afg_file, "{TLE\n");
          fprintf(afg_file, "clr:%u,%u\n", 0, reads[0]->read()->size());
          fprintf(afg_file, "off:0\n");
          fprintf(afg_file, "src:%d\n}\n}\n", reads[0]->read()->id());
          continue;
      }

      int num_reads = reads.size();
      // printf("Reads size: %d\n", num_reads);
      for (int j = 0; j < num_reads - 1; ++j) {
          layout::BetterRead* read1 = reads[j];
          layout::BetterRead* read2 = reads[j + 1];
          read1->Finalize();
          // doesn't work - overlaps.size is 0
          //const std::vector< std::pair< uint32_t, layout::BetterOverlap* >> &overlaps = read1->overlaps();

          // find overlap between first and second read
/*          for (const auto& overlap: overlaps) {
              if (overlap.first == read2->id() && overlap.second != nullptr) {
                  fprintf(afg_file, "{TLE\n");
                  fprintf(afg_file, "clr:%u,%u\n", 0, read1->read()->size());
                  fprintf(afg_file, "off:%u\n", offset);
                  fprintf(afg_file, "src:%d\n}\n", read1->read()->id()); 
                  offset += read1->read()->size() - overlap.second->Length();
                  break;
              }
          }*/

          // workaround - access overlap in overlaps_map
          std::map< std::pair<uint32_t, uint32_t>, uint32_t>::iterator it;
          // printf("%d %d\n", read1->read()->id(),read2->read()->id());
          // why opposite (read2, read1) instead of (read1, read2) ????
          it = overlaps_map.find(std::make_pair(read2->read()->id(), read1->read()->id()));
          if (it != overlaps_map.end()) {
              fprintf(afg_file, "{TLE\n");
              fprintf(afg_file, "clr:%u,%u\n", 0, read1->read()->size());
              fprintf(afg_file, "off:%u\n", offset);
              fprintf(afg_file, "src:%d\n}\n", read1->read()->id()); 
              offset += read1->read()->size() - it->second;
          }
      }
      layout::BetterRead* read = reads[num_reads - 1];
      fprintf(afg_file, "{TLE\n");
      fprintf(afg_file, "clr:%u,%u\n", 0, read->read()->size());
      fprintf(afg_file, "off:%u\n", offset);
      fprintf(afg_file, "src:%d\n}\n", read->read()->id());

      fprintf(afg_file, "}\n");
  }

  fclose(afg_file);
  // end output afg

  fclose(overlaps_file);
  fclose(graphviz_file);
  return 0;
}
