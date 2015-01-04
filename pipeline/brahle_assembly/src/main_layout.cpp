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
#include <string>
#include <iostream>
#include <fstream>
using std::string;
using std::cout;

void usage(char *argv[]) {
  fprintf(
      stderr,
      "Usage: %s <reads_file> <overlaps_file.afg>\n",
      argv[0]);
  fprintf(stderr, "\n");
  fprintf(stderr, "Flags\n");
  fprintf(stderr, "\t-f\t reads provided in fastq format, instead of default, AMOS format\n");
  fprintf(stderr, "\n");
  exit(1);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    usage(argv);
  }

  char reads_file_name[1024] = {0};
  char overlaps_file_name[1024] = {0};

  snprintf(reads_file_name, sizeof(reads_file_name), "%s", argv[1]);
  snprintf(overlaps_file_name, sizeof(overlaps_file_name), "%s", argv[2]);

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
  if (strlen(reads_file_name) > 0) {
    fprintf(stderr, "Reading reads from afg file '%s'\n", reads_file_name);
    reads.reset(layout::ReadReadsAfg(reads_file));
  } else {
    fprintf(stderr, "No input file with reads provided\n");
    usage(argv);
  }
  fprintf(stderr, "Number of reads = %u\n", reads->size());

  // getting overlaps
  std::shared_ptr< overlap::OverlapSet > overlaps;
  if (strlen(overlaps_file_name) > 0) {
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

  // create initial dotgraph
  {
    std::ofstream initial_graph_file;
    string initial_graph = layout::dot_graph(reads.get(), overlaps.get());
    initial_graph_file.open("all_overlaps.dot", std::fstream::out);
    initial_graph_file << initial_graph;
    initial_graph_file.close();
  }

  // create dotgraph after removing transitive edges and containment reads
  {
    std::ofstream no_transitives_graph;
    no_transitives_graph.open("no_transitives.dot", std::fstream::out);
    layout::BetterReadSet brs(reads.get(), false);
    no_transitives_graph << layout::dot_graph(&brs, u->noTransitives().get());
    no_transitives_graph.close();
  }

  start = clock();
  layout::Graph g = layout::Graph::create(u->readSet(), u->noTransitives());
  fprintf(
      stderr,
      "String graph constructed in %.2lfs\n",
      (clock() - start) / static_cast<double>(CLOCKS_PER_SEC));

  g.printToGraphviz(graphviz_file);

  // trimming
  // @mculinovic
  const uint32_t trimSeqLenThreshold = 300;
  g.trim(trimSeqLenThreshold);
  // g.removeBubbles();

  typedef std::shared_ptr< layout::BetterOverlapSet > BetterOverlapSetPtr;
  overlap::ReadSet* read_set = g.extractReads();
  BetterOverlapSetPtr overlap_set = g.extractOverlaps();

  fprintf(stderr, "Number of reads after trimmming: %d\n", read_set->size());
  fprintf(stderr, "Number of overlaps after trimming: %d\n", overlap_set->size());

  u->makeContigs(overlap_set, read_set);

  n50_value = layout::n50(u->contigs());
  fprintf(stderr, "After trimming n50 = %d\n", n50_value);

  // create dotgraph after trimming
  {
    std::ofstream after_trimming_graph;
    after_trimming_graph.open("after_trimming.dot", std::fstream::out);
    layout::BetterReadSet brs(reads.get(), false);
    after_trimming_graph << layout::dot_graph(&brs, g.extractOverlaps().get());
    after_trimming_graph.close();
  }

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
