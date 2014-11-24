// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_BETTER_OVERLAP_H_
#define LAYOUT_BETTER_OVERLAP_H_

#include <overlap/read.h>
#include <overlap/overlap.h>
#include <cassert>
#include <vector>

namespace layout {

const double EPSILON = 0.15;
const double ALPHA = 3;

/**
 * This function compares two numbers with given precision. i.e. if
 * y - eps <= x <= y + eps holds true.
 */
inline bool eq(double x, double y, double eps);

/**
 * This is extending functionality of overlap::Overlap by providing some
 * useful functions that are used in layout step.
 */
class BetterOverlap {
  typedef overlap::ReadSet* ReadSetPtr;
  typedef overlap::Overlap* OverlapPtr;

 public:
  /**
   * Default constructor. Creates an empty overlap.
   */
  BetterOverlap();

  /**
   * Constructs an overlap from an old overlap.
   * Most of the times, you want to use this constructor.
   */
  BetterOverlap(OverlapPtr overlap, ReadSetPtr read_set);

  /**
   * Similar to a copy constructor, but takes a pointer instead.
   */
  explicit BetterOverlap(BetterOverlap* other);

  /**
   * Destructor.
   */
  virtual ~BetterOverlap();

  /**
   * Returns a pointer to the overlap::Overlap.
   */
  const OverlapPtr& overlap() const { return overlap_; }

  /**
   * Returns the overlap::ReadSet the reads from this overlap are from.
   */
  const ReadSetPtr& readSet() const { return read_set_; }

  /**
   * Set the overlap::ReadSet the reads from this overlap are from.
   */
  void setReadSet(ReadSetPtr read_set) { read_set_ = read_set; }

  /**
   * Gets the first read.
   */
  const overlap::Read* one() const;

  /**
   * Gets the second read.
   */
  const overlap::Read* two() const;

  /**
   * General method for getting a read.
   */
  const overlap::Read* get(uint32_t read) const;

  /**
   * This shows what is considered to be a hanigng part. :
   *  ============= Read A ==========================>
   *  \________f.hang(A)_________/|||| Overlap f |||||
   *                              =============== Read B =================>
   *                                                   \____f.hang(B)_____/
   */
  const uint32_t Hang(uint32_t read) const;

  /**
   * Length of the overlap.
   */
  const uint32_t Length() const;

  /**
   * Length of the read.
   */
  const uint32_t Length(uint32_t read) const;

  /**
   * Returns which end of the read this overlap is using.
   * 0 if it is using begining, 1 if it using the end.
   */
  const uint32_t Suf(uint32_t read) const;

  /**
   * Returns the ID of the other read.
   */
  const uint32_t Other(uint32_t read) const;

 private:
  OverlapPtr overlap_;
  ReadSetPtr read_set_;
};

typedef BetterOverlap* BetterOverlapPtr;

class BetterOverlapSet;

/**
 * Helper class for iterating over BetterOverlapSet.
 */
class BetterOverlapSetIter {
 public:
  BetterOverlapSetIter(const BetterOverlapSet*, int);
  bool operator!=(const BetterOverlapSetIter&) const;
  bool operator==(const BetterOverlapSetIter&) const;
  BetterOverlapPtr operator*() const;
  const BetterOverlapSetIter& operator++();
  const BetterOverlapSetIter& operator++(int);
 private:
  const BetterOverlapSet* better_overlap_set_;
  int position_;
};

/**
 * A collection of BetterOverlaps.
 *
 * Use this to store the overlaps.
 */
class BetterOverlapSet {
  typedef BetterOverlap* BetterOverlapPtr;

 public:
  /**
   * Creates an empty overlap set.
   */
  explicit BetterOverlapSet(overlap::ReadSet* read_set);

  /**
   * Creates a BetterOverlapSet from a given overlap::OverlapSet.
   */
  BetterOverlapSet(
      overlap::ReadSet* read_set,
      overlap::OverlapSet* overlap_set);

  /**
   * Destructor.
   */
  virtual ~BetterOverlapSet();

  /**
   * Returns the read set for this overlap set.
   */
  const overlap::ReadSet* readSet() const { return read_set_; }

  /**
   * Returns the i-th overlap.
   */
  BetterOverlapPtr& operator[](size_t);

  /**
   * Returns the i-th overlap (const version).
   */
  const BetterOverlapPtr& operator[](size_t) const;

  /**
   * Returns the number of overlaps in this overlap set.
   */
  const size_t size() const;

  /**
   * Returns an iterator pointing to the beginning of the overlap set.
   */
  BetterOverlapSetIter begin() const;

  /**
   * Returns an iterator pointing to the end of the overlap set.
   */
  BetterOverlapSetIter end() const;

  /**
   * Adds an overlap to the overlap set.
   */
  void Add(BetterOverlapPtr);
  void Add(overlap::Overlap*);

 private:
  overlap::ReadSet* read_set_;
  std::vector< BetterOverlapPtr > overlap_set_;
};

};  // namespace layout

#endif  // LAYOUT_BETTER_OVERLAP_H_
