// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_CONTIG_H_
#define LAYOUT_CONTIG_H_

#include <layout/better_overlap.h>
#include <layout/better_read.h>
#include <deque>
#include <vector>

namespace layout {

/**
 * A class used to describe contigs.
 */
class Contig {
  typedef BetterRead* BetterReadPtr;
  typedef BetterReadSet* BetterReadSetPtr;
  typedef BetterOverlap* BetterOverlapPtr;

 public:
  /**
   * Creates a contig from a single read.
   */
  Contig(BetterReadPtr starting, BetterReadSetPtr read_set);

  /**
   * Destructor.
   */
  virtual ~Contig();

  /**
   * Size of the contig (number of reads).
   */
  size_t size() const;

  /**
   * Joins another contig to this one. The other contig is then marked as
   * not usable.
   */
  void Join(BetterOverlapPtr better_overlap, Contig* contig);

  /**
   * Is the contig usable. True if it was not joined to another contig.
   */
  bool IsUsable() const;

  /**
   * Sets the valid state for the contig. Use this to mark the usability.
   */
  void SetValid(bool value);

 private:
  std::deque< BetterReadPtr > reads_;
  BetterReadSetPtr read_set_;
  bool valid_;
  bool alive_;
  bool IsLeftOverlap(BetterOverlapPtr better_overlap) const;
  void Kill();
};


/**
 * A set of contigs.
 */
class ContigSet {
  typedef Contig* ContigPtr;

 public:
  /**
   * Creates a contig for each read from the given set.
   */
  explicit ContigSet(BetterReadSet* read_set);

  /**
   * Destructor.
   */
  virtual ~ContigSet();

  /**
   * Returns the number of contigs in the set.
   */
  size_t size() const;

  /**
   * Returns the i-th contig.
   */
  ContigPtr& operator[](int i);

  /**
   * Returns the i-th contig (conts version).
   */
  const ContigPtr& operator[](int i) const;

 private:
  std::vector< ContigPtr > contigs_;
};

}  // namespace layout

#endif  // LAYOUT_CONTIG_H_
