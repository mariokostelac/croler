// Copyright 2014 Bruno Rahle

#include "test/unitigging_test.h"
#include <overlap/read.h>
#include <overlap/overlap.h>
#include <layout/better_overlap.h>
#include <layout/contig.h>
#include <layout/layout_utils.h>

#include <cstdio>
#include <cstring>
#include <ctime>
#include <typeinfo>

namespace test {

UnitiggingTest::UnitiggingTest() {
}

UnitiggingTest::~UnitiggingTest() {
}

overlap::Read* UnitiggingTest::makeRead(const char* data) {
  size_t len = strlen(data);
  uint8_t* uint_data = new uint8_t[len];
  memcpy(uint_data, data, len);
  return new overlap::Read(uint_data, len-1, 0, 0);
}

UnitiggingIsTransitiveTest::UnitiggingIsTransitiveTest() {
}

UnitiggingIsTransitiveTest::~UnitiggingIsTransitiveTest() {
}

bool UnitiggingIsTransitiveTest::run() {
  overlap::ReadSet read_set(3);
  read_set.Add(makeRead("AAAAAAAAAABBBBBBBBBBBCCCCCCC"));
  read_set.Add(makeRead("BBBBBBBBBBCCCCCCCCCCCDDDDD"));
  read_set.Add(makeRead("CCCCCCCCCCDDDDDDDDDDDAAA"));
  overlap::OverlapSet overlap_set(4);
  overlap_set.Add(new overlap::Overlap(0, 1, 17, 17, overlap::Overlap::Type::EB, 17));
  overlap_set.Add(new overlap::Overlap(1, 2, 15, 15, overlap::Overlap::Type::EB, 15));
  overlap_set.Add(new overlap::Overlap(0, 2, 7, 7, overlap::Overlap::Type::EB, 7));
  overlap_set.Add(new overlap::Overlap(2, 0, 3, 3, overlap::Overlap::Type::EB, 3));
  layout::Unitigging* unitigging = new layout::Unitigging(&read_set, &overlap_set);
  if (!unitigging->isTransitive(
          unitigging->overlaps_[2],
          unitigging->overlaps_[0],
          unitigging->overlaps_[1])) {
    delete unitigging;
    return false;
  }
  if (unitigging->isTransitive(
          unitigging->overlaps_[3],
          unitigging->overlaps_[1],
          unitigging->overlaps_[0])) {
    delete unitigging;
    return false;
  }
  delete unitigging;
  return true;
}

UnitiggingContainmentTest::UnitiggingContainmentTest() {
}

UnitiggingContainmentTest::~UnitiggingContainmentTest() {
}

bool UnitiggingContainmentTest::run() {
  overlap::ReadSet read_set(3);
  read_set.Add(makeRead("AAAAAAAAAABBBBBBBBBBBCCCCCCC"));
  read_set.Add(makeRead("BBBBBBBBBBCCCCCC"));
  read_set.Add(makeRead("CCCCCCCCCCDDDDDDDDDDDAAA"));
  overlap::OverlapSet overlap_set(3);
  overlap_set.Add(new overlap::Overlap(0, 1, 15, 15, overlap::Overlap::Type::EB, 15));
  overlap_set.Add(new overlap::Overlap(1, 2, 5, 5, overlap::Overlap::Type::EB, 5));
  overlap_set.Add(new overlap::Overlap(0, 2, 7, 7, overlap::Overlap::Type::EB, 7));
  layout::Unitigging* unitigging = new layout::Unitigging(&read_set, &overlap_set);
  unitigging->removeContainmentEdges();
  if (unitigging->no_contains_->size() != 1) {
    delete unitigging;
    return false;
  }
  delete unitigging;
  return true;
}

UnitiggingTransitiveTest::UnitiggingTransitiveTest() {
}

UnitiggingTransitiveTest::~UnitiggingTransitiveTest() {
}

bool UnitiggingTransitiveTest::run() {
  overlap::ReadSet read_set(3);
  read_set.Add(makeRead("AAAAAAAAAABBBBBBBBBBBCCCCCCC"));
  read_set.Add(makeRead("BBBBBBBBBBCCCCCCCCCCCDDDDD"));
  read_set.Add(makeRead("CCCCCCCCCCDDDDDDDDDDDAAA"));
  overlap::OverlapSet overlap_set(3);
  overlap_set.Add(new overlap::Overlap(0, 1, 17, 17, overlap::Overlap::Type::EB, 17));
  overlap_set.Add(new overlap::Overlap(1, 2, 15, 15, overlap::Overlap::Type::EB, 15));
  overlap_set.Add(new overlap::Overlap(0, 2, 7, 7, overlap::Overlap::Type::EB, 7));
  layout::Unitigging* unitigging = new layout::Unitigging(&read_set, &overlap_set);
  unitigging->removeContainmentEdges();
  unitigging->removeTransitiveEdges();
  if (unitigging->no_transitives_->size() != 2) {
    delete unitigging;
    return false;
  }
  delete unitigging;
  return true;
}

UnitiggingContigTest::UnitiggingContigTest() {
}

UnitiggingContigTest::~UnitiggingContigTest() {
}

bool UnitiggingContigTest::run() {
  overlap::ReadSet read_set(2);
  read_set.Add(makeRead("AAAAAAAAAABBBBBBBBBBBCCCCCCC"));
  read_set.Add(makeRead("BBBBBBBBBBCCCCCCCCCCCDDDDD"));
  overlap::OverlapSet overlap_set(1);
  overlap_set.Add(new overlap::Overlap(0, 1, 17, 17, overlap::Overlap::Type::EB, 17));
  layout::Unitigging* unitigging = new layout::Unitigging(&read_set, &overlap_set);
  unitigging->removeContainmentEdges();
  unitigging->removeTransitiveEdges();
  unitigging->makeContigs();
  int contig_sizes[2] = {2, 0};
  for (int i = 0; i < unitigging->contigs_->size(); ++i) {
    auto contig = (*(unitigging->contigs_))[i];
    if (contig->size() != contig_sizes[i]) {
      return false;
    }
  }
  delete unitigging;
  return true;
}

UnitiggingCompleteTest::UnitiggingCompleteTest(
    const char *read_file,
    const char *overlap_file) {
  rfd_ = fopen(read_file, "r");
  ofd_ = fopen(overlap_file, "r");
}

UnitiggingCompleteTest::~UnitiggingCompleteTest() {
  fclose(rfd_);
  fclose(ofd_);
}

bool UnitiggingCompleteTest::run() {
  auto read_set = layout::ReadReadsAfg(rfd_);
  auto overlap_set = layout::ReadOverlapsAfg(read_set, ofd_);
  auto unitigging = new layout::Unitigging(read_set, overlap_set);
  double start = clock();
  unitigging->start();
  printf("Unitigging done in %.2lfs\n", (clock() - start)/CLOCKS_PER_SEC);
  auto contigs = unitigging->contigs();
  int number = 0;
  for (size_t i = 0; i < contigs->size(); ++i) {
    auto contig = (*contigs)[i];
    if (!contig->IsUsable()) {
      continue;
    }
    ++number;
    printf("%d: %d\n", i, contig->size());
  }
  printf("Number of contigs = %d\n", number);
  delete unitigging;
  delete read_set;
  delete overlap_set;
  return true;
}

UnitiggingTestRunner::UnitiggingTestRunner() {
}

UnitiggingTestRunner::~UnitiggingTestRunner() {
  for (auto test : tests_) {
    delete test;
  }
}

void UnitiggingTestRunner::run() {
  int total = 0;
  int successful = 0;

  for (auto test : tests_) {
    double test_start = clock();
    bool result = test->run();
    if (result == false) {
      fprintf(
          stderr,
          "\e[31mFAILED: %s (%.2lfs) \033[m\n",
          typeid(*test).name(),
          (clock()-test_start)/CLOCKS_PER_SEC);
    } else {
      fprintf(
          stderr,
          "\e[32mOK: %s (%.2lfs) \033[m\n",
          typeid(*test).name(),
          (clock()-test_start)/CLOCKS_PER_SEC);
    }
    total += 1;
    successful += result;
  }

  if (total == successful) {
    printf("\e[32mALL OK!\033[m (%d/%d)\n", total, successful);
  } else {
    printf(
        "\e[31mFAILURE: %d tests FAILED!\033[m (%d/%d = %.2lf%% OK)\n",
        total-successful,
        successful,
        total,
        ((double)successful)/total);
  }
}

void UnitiggingTestRunner::addTest(UnitiggingTest* test) {
  tests_.emplace_back(test);
}

};  // namespace test

int main() {
  int success = 0;
  int total = 0;
  test::UnitiggingTestRunner ut;
  ut.addTest(new test::UnitiggingIsTransitiveTest());
  ut.addTest(new test::UnitiggingContainmentTest());
  ut.addTest(new test::UnitiggingTransitiveTest());
  ut.addTest(new test::UnitiggingContigTest());
  //  ovo je jako "porculansko"
  ut.addTest(new test::UnitiggingCompleteTest(
      "sample/small/minimus_results/reads.afg",
      "sample/small/minimus_results/overlaps.afg"));
  ut.addTest(new test::UnitiggingCompleteTest(
      "sample/large/minimus_results/reads.afg",
      "sample/large/minimus_results/overlaps.afg"));
  ut.run();

  return 0;
}
