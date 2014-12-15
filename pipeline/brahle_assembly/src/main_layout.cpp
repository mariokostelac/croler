// Copyright 2014 Bruno Rahle
#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/unitigging.h>
#include <layout/layout_utils.h>
#include <layout/string_graph.h>
#include <layout/contig.h>
#include <layout/better_read.h>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <memory>
#include <deque>
#include <vector>
#include <map>
#include <utility>

char amos_bank_name[1024] = {0};

void usage(char *argv[]) {
  fprintf(
      stderr,
      "Usage: %s <amos_bank> | (<reads_file.seq> <overlaps_file.afg>)\n",
      argv[0]);
  exit(0);
}

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    usage(argv);
  }
  char reads_file_name[1024] = {0};
  char overlaps_file_name[1024] = {0};

  if (argc == 2) {
    snprintf(amos_bank_name, sizeof(amos_bank_name), "%s", argv[1]);
  } else {
    snprintf(reads_file_name, sizeof(reads_file_name), "%s", argv[1]);
    snprintf(overlaps_file_name, sizeof(overlaps_file_name), "%s", argv[2]);
  }

  FILE *reads_file = nullptr;
  FILE *overlaps_file = nullptr;
  FILE *graphviz_file = nullptr;

  if (strlen(reads_file_name)) {
    reads_file = fopen(reads_file_name, "r");
    if (reads_file == nullptr) {
      fprintf(
          stderr,
          "ERROR: reads file ('%s') cannot be opened!\n",
          reads_file_name);
      exit(1);
    }
    fclose(reads_file);
  }

  if (strlen(overlaps_file_name)) {
    overlaps_file = fopen(overlaps_file_name, "r");
    if (overlaps_file == nullptr) {
      fprintf(
          stderr,
          "ERROR: overlaps file ('%s') cannot be opened!\n",
          overlaps_file_name);
      exit(1);
    }
  }

  graphviz_file = fopen("graph.dot", "w");
  if (graphviz_file == nullptr) {
    fprintf(
        stderr,
        "ERROR: graphviz file ('%s') cannot be opened!\n",
        "graph.dot");
    exit(1);
  }

  // getting reads
  std::shared_ptr< overlap::ReadSet > reads;
  if (strlen(amos_bank_name) > 0) {
    fprintf(stderr, "Reading reads from bank %s", amos_bank_name);
    reads.reset(layout::ReadReadsAmos(amos_bank_name));
  } else {
    fprintf(stderr, "Reading reads from fasta file %s", reads_file_name);
    reads.reset(layout::ReadReadsSeq(reads_file_name));
  }
  fprintf(stderr, "Number of reads = %u\n", reads->size());

  // getting overlaps
  std::shared_ptr< overlap::OverlapSet > overlaps;
  if (strlen(amos_bank_name) > 0) {
    fprintf(stderr, "Reading overlaps from bank %s", amos_bank_name);
    overlaps.reset(layout::ReadOverlapsAmos(reads.get(), amos_bank_name));
  } else {
    fprintf(stderr, "Reading from fasta file %s", reads_file_name);
    overlaps.reset(layout::ReadOverlapsAfg(reads.get(), overlaps_file));
  }
  fprintf(stderr, "Number of overlaps = %u\n", overlaps->size());

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
  FILE *afg_file = fopen("layout.afg", "w");

  std::shared_ptr<layout::ContigSet> contigs = u->contigs();
  int contigs_size = contigs->size();
  for (int i = 0; i < contigs_size; ++i) {
      // skip non-usable contigs
      if (!((*contigs)[i]->IsUsable())) continue;

      fprintf(afg_file, "{LAY\n");
      uint32_t offset = 0;
      const std::deque< layout::BetterRead* > &reads = (*contigs)[i]->getReads();

      int num_reads = reads.size();
      for (int j = 0; j < num_reads - 1; ++j) {
          layout::BetterRead* read1 = reads[j];
          layout::BetterRead* read2 = reads[j + 1];
          read1->Finalize();
          const std::vector< std::pair< uint32_t, layout::BetterOverlap* >> &overlaps = read1->overlaps();

          // find overlap between first and second read
          for (const auto& overlap: overlaps) {
              if (overlap.first == read2->id() && overlap.second != nullptr) {
                  fprintf(afg_file, "{TLE\n");
                  fprintf(afg_file, "clr:%u,%u\n", 0, read1->read()->size());
                  fprintf(afg_file, "off:%u\n", offset);
                  fprintf(afg_file, "src:%d\n}\n", read1->read()->id()); 
                  offset += read1->read()->size() - overlap.second->Length();
                  break;
              }
          }
      }

      // output last read
      layout::BetterRead* read = reads[num_reads - 1];
      fprintf(afg_file, "{TLE\n");
      fprintf(afg_file, "clr:%u,%u\n", 0, read->read()->size());
      fprintf(afg_file, "off:%u\n", offset);
      fprintf(afg_file, "src:%d\n}\n", read->read()->id());

      fprintf(afg_file, "}\n");
  }

  if (afg_file != nullptr) {
    fclose(afg_file);
  }

  if (overlaps_file != nullptr) {
    fclose(overlaps_file);
  }

  if (graphviz_file != nullptr) {
    fclose(graphviz_file);
  }

  return 0;
}
