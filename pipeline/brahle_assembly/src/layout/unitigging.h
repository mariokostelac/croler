// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_UNITIGGING_H_
#define LAYOUT_UNITIGGING_H_

#include <overlap/read.h>
#include <overlap/overlap.h>

#include <layout/better_overlap.h>
#include <layout/better_read.h>
#include <layout/contig.h>

#include <memory>

namespace test {
class UnitiggingTest;
class UnitiggingIsTransitiveTest;
class UnitiggingContainmentTest;
class UnitiggingTransitiveTest;
class UnitiggingContigTest;
};  // namespace test

namespace layout {

/**
 * Class that solves the unitigging problem.
 *
 * Call Unitigging::start() to actually do the solving.
 */
class Unitigging {
 public:
  /**
   * Constructor - requires a read set and an overlap set.
   */
  Unitigging(
      overlap::ReadSet* reads,
      overlap::OverlapSet* overlaps);

  /**
   * Destructor.
   */
  virtual ~Unitigging();

  typedef std::shared_ptr< BetterOverlapSet > BetterOverlapSetPtr;
  typedef std::shared_ptr< BetterReadSet > BetterReadSetPtr;
  typedef std::shared_ptr< ContigSet > ContigSetPtr;

  /**
   * Call this to actaully do the computations.
   */
  void start();

  /**
   * Getter for the contig set.
   *
   * Available after the start method has been completed.
   */
  ContigSetPtr& contigs();

  /**
   * Getter for the set of overlaps without conatained reads.
   *
   * Available after the start method has been completed.
   */
  const BetterOverlapSetPtr& noContains() const;

  /**
   * Getter for the set of overlaps with transitive edges removed.
   *
   * Available after the start method has been completed.
   */
  const BetterOverlapSetPtr& noTransitives() const;

  /**
   * Getter for the read set.
   */
  const BetterReadSetPtr& readSet();

 private:
  overlap::ReadSet* reads_;
  overlap::OverlapSet* orig_overlaps_;
  BetterOverlapSet overlaps_;
  BetterOverlapSetPtr no_contains_;
  BetterOverlapSetPtr no_transitives_;
  ContigSetPtr contigs_;
  BetterReadSetPtr better_read_set_;

  void removeContainmentEdges();
  bool isTransitive(
      BetterOverlap* o1,
      BetterOverlap* o2,
      BetterOverlap* o3) const;
  void removeTransitiveEdges();
  void makeContigs();

  friend test::UnitiggingTest;
  friend test::UnitiggingIsTransitiveTest;
  friend test::UnitiggingContainmentTest;
  friend test::UnitiggingTransitiveTest;
  friend test::UnitiggingContigTest;
};

}  // namespace layout

#endif  // LAYOUT_UNITIGGING_H_
