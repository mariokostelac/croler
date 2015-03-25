// Copyright 2014 Bruno Rahle
#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/unitigging.h>
#include <layout/layout_utils.h>
#include <layout/string_graph.h>
#include <layout/contig.h>
#include <layout/better_read.h>
#include <lib/parsero/parsero.h>

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

const uint32_t EXPECT_READS = 1 << 10;
uint32_t TRIM_ROUNDS = 1;
// trimming read length threshold
uint32_t READ_LEN_THRESHOLD = 300;

uint32_t BUBBLE_ROUNDS = 1;
// maximum number of bfs nodes in bubble
uint32_t MAX_NODES = 500;
// maximum walk sequence length in bubble
uint64_t MAX_DISTANCE = 5000;
// maximum number of bubble walks
uint32_t MAX_WALKS = 10;
// maximum diff between walk sequences after alignment
double MAX_DIFF = 0.2;

char *reads_file_name = nullptr;
char *overlaps_file_name = nullptr;

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

void setup_cmd_interface(int argc, char **argv) {

  parsero::add_option("t:", "number of trimming rounds",
    [] (char *option) { TRIM_ROUNDS = atoi(option); });

  parsero::add_option("b:", "number of bubble popping rounds",
    [] (char *option) { BUBBLE_ROUNDS = atoi(option); });

  parsero::add_option("h:", "trimming read length threshold",
    [] (char *option) { READ_LEN_THRESHOLD = atoi(option); });

  parsero::add_option("n:", "maximum number of nodes during bfs in bubble popping",
    [] (char *option) { MAX_NODES = atoi(option); });

  parsero::add_option("d:", "maximum walk sequence length in bubble",
    [] (char *option) { MAX_DISTANCE = atoi(option); });

  parsero::add_option("w:", "maximum number of walks in bubble",
    [] (char *option) { MAX_WALKS = atoi(option); });

  parsero::add_option("a:", "maximum diff between aligned bubble walk sequences",
    [] (char *option) { MAX_DIFF = atof(option); });

  parsero::add_argument("reads.afg",
    [] (char *filename) { reads_file_name = filename; });

  parsero::add_argument("overlaps.afg",
    [] (char *filename) { overlaps_file_name = filename; });

  parsero::parse(argc, argv);
}

int main(int argc, char *argv[]) {

  setup_cmd_interface(argc, argv);

  if (reads_file_name == nullptr || overlaps_file_name == nullptr) {
    parsero::help(argv[0]);
    exit(1);
  }

  FILE *overlaps_file = nullptr;
  FILE *graphviz_file = nullptr;

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
  overlap::ReadSet reads(EXPECT_READS);
  if (strlen(reads_file_name) > 0) {
    fprintf(stderr, "Reading reads from afg file '%s'\n", reads_file_name);
    layout::ReadReadsAfg(reads, reads_file_name);
  } else {
    fprintf(stderr, "No input file with reads provided\n");
    usage(argv);
  }
  fprintf(stderr, "Number of reads = %u\n", reads.size());

  // getting overlaps
  std::shared_ptr< overlap::OverlapSet > overlaps;
  if (strlen(overlaps_file_name) > 0) {
    fprintf(stderr, "Reading overlaps from afg file '%s'\n", overlaps_file_name);
    overlaps.reset(layout::ReadOverlapsAfg(&reads, overlaps_file));
  } else {
    fprintf(stderr, "No input file with overlaps provided\n");
    usage(argv);
  }
  fprintf(stderr, "Number of overlaps = %u\n", overlaps->size());

  clock_t start = clock();
  std::shared_ptr< layout::Unitigging > u(
      new layout::Unitigging(&reads, overlaps.get()));
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
    string initial_graph = layout::dot_graph(&reads, overlaps.get());
    initial_graph_file.open("all_overlaps.dot", std::fstream::out);
    initial_graph_file << initial_graph;
    initial_graph_file.close();
  }

  // output .afg files after removing transitive edges and containment edges
  {
    layout::BetterReadSet brs(&reads, false);

    std::ofstream no_containment_graph;
    no_containment_graph.open("no_containment.afg", std::fstream::out);
    layout::write_overlaps(no_containment_graph, &brs, u->noContains().get());
    no_containment_graph.close();

    std::ofstream no_transitives_graph;
    no_transitives_graph.open("no_transitives.afg", std::fstream::out);
    layout::write_overlaps(no_transitives_graph, &brs, u->noTransitives().get());
    no_transitives_graph.close();
  }

  // create dotgraph after removing transitive edges and containment reads
  {
    std::ofstream no_transitives_graph;
    no_transitives_graph.open("no_transitives.dot", std::fstream::out);
    layout::BetterReadSet brs(&reads, false);
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

  // simplification
  // @mculinovic

  while (TRIM_ROUNDS-- > 0) {
    g.trim(READ_LEN_THRESHOLD);
  }

  while (BUBBLE_ROUNDS-- > 0) {
    g.removeBubbles(MAX_NODES, MAX_DISTANCE, MAX_WALKS, MAX_DIFF);
  }

  typedef std::shared_ptr< layout::BetterOverlapSet > BetterOverlapSetPtr;
  overlap::ReadSet* read_set = g.extractReads();
  BetterOverlapSetPtr overlap_set = g.extractOverlaps();

  fprintf(stderr, "Number of reads after graph simplification: %d\n", read_set->size());
  fprintf(stderr, "Number of overlaps after graph simplification: %d\n", overlap_set->size());

  u->makeContigs(overlap_set, read_set);

  n50_value = layout::n50(u->contigs());
  fprintf(stderr, "After simplification n50 = %d\n", n50_value);

  // create dotgraph after trimming
  {
    std::ofstream after_trimming_graph;
    after_trimming_graph.open("after_trimming.dot", std::fstream::out);
    layout::BetterReadSet brs(&reads, false);
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
