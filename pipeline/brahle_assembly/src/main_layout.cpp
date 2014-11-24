// Copyright 2014 Bruno Rahle
#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/unitigging.h>
#include <layout/layout_utils.h>
#include <layout/string_graph.h>
#include <cstdlib>
#include <ctime>
#include <memory>

void usage(char *argv[]) {
  fprintf(
      stderr,
      "Usage: %s <minimus_folder> | (<reads_file.afg> <overlaps_file.afg>)\n",
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

  std::shared_ptr< overlap::ReadSet > reads(layout::ReadReadsAfg(reads_file));
  fprintf(stderr, "Number of reads = %u\n", reads->size());

  std::shared_ptr< overlap::OverlapSet > overlaps(
      layout::ReadOverlapsAfg(reads.get(), overlaps_file));
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

  fclose(overlaps_file);
  fclose(reads_file);
  fclose(graphviz_file);
  return 0;
}
