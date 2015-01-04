#include "MultipleAligner.h"
#include "FastaFile.h"
#include "LayoutFile.h"

#include "lib/amos/reader.cpp"

#include <cassert>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <unistd.h>

using namespace acgt;
using namespace std;

void usage(char *path) {
  printf(
    "usage: %s [-p] [-e EPSILON] [-m MAX_RADIUS] <reads-afg> <layout-afg>\n"
    "\t-p is no-pack, disable packing non-interscted reads to same rows\n"
    "\t-e (=0.001) EPSILON, "
      "band size in edit-distance as percent of maximum offset\n"
    "\t-m (=300) MAX_RADIUS, maximum band size to use\n", path);
  exit(1);
}

void scan_args(int argc, char **argv,
    float *epsilon, int *max_radius, bool *pack) {
  int opt;
  while ((opt = getopt(argc, argv, "pe:m:")) != -1) {
    switch (opt) {
      case 'p':
        *pack = 0;
        break;
      case 'e':
        sscanf(optarg, "%f", epsilon);
        break;
      case 'm':
        sscanf(optarg, "%d", max_radius);
        break;
      default:
        usage(*argv);
    }
  }
}

int doit(
    float epsilon,
    int max_radius,
    bool pack,
    const char *afg_file,
    const char *layout_file,
    const char *cons_file = "consensus.fasta") {

  FastaFile output(cons_file);
  LayoutFile layout(layout_file);

  vector<AMOS::Read*> reads;
  if (AMOS::get_reads(reads, afg_file) < 0) {
    fprintf(stderr, "Error while reading '%s'\n", afg_file);
    return 1;
  }

  if (!layout.initialize()) {
    fprintf(stderr, "%s", layout.getError().c_str());
    return 1;
  }

  // erase contents of cons_file
  fclose(fopen(cons_file, "w"));

  // map reads so we can get them by real ids
  unordered_map<uint32_t, AMOS::Read*> read_by_id;
  int reads_size = reads.size();
  for (int i = 0; i < reads_size; ++i) {
    const auto& read = reads[i];
    assert(!read_by_id.count(read->iid));
    read_by_id[read->iid] = read;
  }

  auto contigs = layout.getContigs();
  for (const auto& contig : contigs) {
    printf("processing new contig with %ld reads..\n", contig.size());
    MultipleAligner MA(epsilon, max_radius, pack);
    for (auto c_part : contig) {
      AMOS::Read* read = read_by_id[c_part.read_id];
      int size = c_part.hi - c_part.lo; // can be neg!
      if (c_part.hi < c_part.lo) {
        swap(c_part.lo, c_part.hi);
      }
      MA.addSequence(
          read->seq + c_part.lo,
          size,
          c_part.offset);
    }
    output.appendRead(MA.getConsensus(NULL));
  }
  return 0;
}

int main(int argc, char **argv) {
  float epsilon = 0.001;
  int max_radius = 300;
  bool pack = true;
  char *path = *argv;
  {
    scan_args(argc, argv, &epsilon, &max_radius, &pack);
    argv += optind;
    argc -= optind;
  }
  if (argc != 2) {
    fprintf(stderr, "expected two arguments after options.\n\n");
    usage(path);
  }
  return doit(epsilon, max_radius, pack,
      argv[0], argv[1]);
}

