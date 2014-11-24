// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_BETTER_READ_H_
#define LAYOUT_BETTER_READ_H_

#include <overlap/read.h>
#include <layout/better_overlap.h>
#include <memory>
#include <vector>
#include <utility>

namespace layout {

/**
 * An improvement over overlap::Read that adds some additional functionality.
 */
class BetterRead {
  typedef overlap::Read* ReadPtr;
  typedef layout::BetterOverlap* BetterOverlapPtr;
  typedef std::vector< std::pair< uint32_t, BetterOverlapPtr > > OverlapVector;

 public:
  /**
   * Constructs a read that has a given ID.
   */
  explicit BetterRead(uint32_t id);

  /**
   * Constructs a read from overlap::Read.
   */
  BetterRead(uint32_t id, ReadPtr read);

  /**
   * Destructor.
   */
  virtual ~BetterRead();

  /**
   * Gets the overlap::Read().
   */
  const ReadPtr read() const { return read_; }

  /**
   * Gets the ID of the current read.
   */
  const uint32_t id() const { return id_; }

  /**
   * Adds an overlap to this read.
   */
  void AddOverlap(BetterOverlapPtr overlap);

  /**
   * Call this once you finished adding overlaps.
   */
  void Finalize();

  /**
   * Retuns overlaps that are connected to this read.
   */
  const OverlapVector& overlaps() const;

 private:
  uint32_t id_;
  ReadPtr read_;
  OverlapVector overlaps_;
  bool finalized_;
};

class BetterReadSet;

/**
 * A helper iterator class for BetterReadSet.
 */
class BetterReadSetIter {
 public:
  BetterReadSetIter(const BetterReadSet*, int);

  bool operator!=(const BetterReadSetIter&) const;
  bool operator==(const BetterReadSetIter&) const;
  BetterRead* operator*() const;
  const BetterReadSetIter& operator++();
  const BetterReadSetIter& operator++(int);
 private:
  const BetterReadSet* better_read_set_;
  int position_;
};


/**
 * A set of BetterReads. Use it to store a set of reads.
 */
class BetterReadSet {
  typedef overlap::ReadSet* ReadSetPtr;
  typedef BetterRead* BetterReadPtr;

 public:
  /**
   * Constructs an empty set and reserves space for the given number of elements.
   */
  explicit BetterReadSet(size_t);

  /**
   * Constructs a BetterReadSet from the given overlap::ReadSet.
   */
  BetterReadSet(ReadSetPtr, bool);

  /**
   * Destructor.
   */
  virtual ~BetterReadSet();

  /**
   * Returns the read with the given ID.
   */
  BetterReadPtr& operator[](size_t);

  /**
   * Returns the read with the given ID (const version).
   */
  const BetterReadPtr& operator[](size_t) const;

  /**
   * Size of the BetterReadSet (number of reads).
   */
  const size_t size() const;

  /**
   * Iterator pointing to the begining of the BetterReadSet.
   */
  BetterReadSetIter begin() const;

  /**
   * Iterator pointing to the begining of the BetterReadSet.
   */
  BetterReadSetIter end() const;

  /**
   * Call this when you finish creating the read set.
   */
  void Finalize();

 private:
  std::vector< BetterReadPtr > read_set_;
  bool master_;
};

};  // namespace layout

#endif  // LAYOUT_BETTER_READ_H_
