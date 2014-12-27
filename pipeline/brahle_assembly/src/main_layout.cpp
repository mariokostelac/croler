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
#include <unistd.h>
#include <memory>
#include <deque>
#include <vector>
#include <map>
#include <utility>

char amos_bank_name[1024] = {0};

bool READS_TYPE_FASTQ = false;

void usage(char *argv[]) {
  fprintf(
      stderr,
      "Usage: %s <amos_bank> | ([-f] <reads_file> <overlaps_file.afg>)\n",
      argv[0]);
  fprintf(stderr, "\n");
  fprintf(stderr, "Flags\n");
  fprintf(stderr, "\t-f\t reads provided in fastq format, instead of default, AMOS format\n");
  fprintf(stderr, "\n");
  exit(1);
}

void parse_args(int argc, char **argv) {
  int opterr = 0;
  char c;
  while ((c = getopt (argc, argv, "f")) != -1) {
    switch (c) {
      case 'f':
        READS_TYPE_FASTQ = true;
        break;
      default:
        usage(argv);
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 4) {
    usage(argv);
  }
  parse_args(argc, argv);

  char reads_file_name[1024] = {0};
  char overlaps_file_name[1024] = {0};

  if (argc == 2) {
    snprintf(amos_bank_name, sizeof(amos_bank_name), "%s", argv[1]);
  } else if (argc == 3) {
    snprintf(reads_file_name, sizeof(reads_file_name), "%s", argv[1]);
    snprintf(overlaps_file_name, sizeof(overlaps_file_name), "%s", argv[2]);
  } else if (argc == 4) {
    snprintf(reads_file_name, sizeof(reads_file_name), "%s", argv[2]);
    snprintf(overlaps_file_name, sizeof(overlaps_file_name), "%s", argv[3]);
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
    fprintf(stderr, "Reading reads from bank '%s'\n", amos_bank_name);
    reads.reset(layout::ReadReadsAmos(amos_bank_name));
  } else if (strlen(reads_file_name) > 0) {
    if (READS_TYPE_FASTQ) {
      fprintf(stderr, "Reading reads from fasta file '%s'\n", reads_file_name);
      reads.reset(layout::ReadReadsSeq(reads_file_name));
    } else {
      fprintf(stderr, "Reading reads from afg file '%s'\n", reads_file_name);
      reads.reset(layout::ReadReadsAfg(reads_file));
    }
  } else {
    fprintf(stderr, "No input file with reads provided\n");
    usage(argv);
  }
  fprintf(stderr, "Number of reads = %u\n", reads->size());

  // getting overlaps
  std::shared_ptr< overlap::OverlapSet > overlaps;
  if (strlen(amos_bank_name) > 0) {
    fprintf(stderr, "Reading overlaps from bank '%s'\n", amos_bank_name);
    overlaps.reset(layout::ReadOverlapsAmos(reads.get(), amos_bank_name));
  } else if (strlen(overlaps_file_name) > 0) {
    fprintf(stderr, "Reading overlaps from afg file '%s'\n", overlaps_file_name);
    overlaps.reset(layout::ReadOverlapsAfg(reads.get(), overlaps_file));
  } else {
    fprintf(stderr, "No input file with overlaps provided\n");
    usage(argv);
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

  // trimming
  // @mculinovic
  g.trim();

  typedef std::shared_ptr< layout::BetterOverlapSet > BetterOverlapSetPtr;
  overlap::ReadSet* read_set = g.extractReads();
  BetterOverlapSetPtr overlap_set = g.extractOverlaps();

  fprintf(stderr, "Number of reads after trimmming: %d\n", read_set->size());
  fprintf(stderr, "Number of overlaps after trimming: %d\n", overlap_set->size());

  u->makeContigs(overlap_set, read_set);

  n50_value = layout::n50(u->contigs());
  fprintf(stderr, "After trimming n50 = %d\n", n50_value);

  std::shared_ptr<layout::ContigSet> contigs = u->contigs();
  int written = ContigsToFile(contigs, "layout.afg");
  fprintf(stderr, "Written %d contigs to a file 'layout.afg'\n", written);

  if (overlaps_file != nullptr) {
    fclose(overlaps_file);
  }

  if (graphviz_file != nullptr) {
    fclose(graphviz_file);
  }

  return 0;
}
