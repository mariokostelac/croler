// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_LABEL_H_
#define LAYOUT_LABEL_H_

#include <layout/better_read.h>

#include <string>

namespace layout {

/**
 * Used to store labels in a string graph.
 */
class Label {
 public:
  /**
   * Direction of the edge - it can either be from read one to read two,
   * or from read two to read one.
   */
  enum Direction {
    FROM_ONE_TO_TWO,
    FROM_TWO_TO_ONE
  };

  /**
   * Creates a label that will represent the given overlap and direction.
   */
  Label(BetterOverlapPtr overlap, Direction direction);

  /**
   * Destructor.
   */
  virtual ~Label();

  /**
   * Returns the string representation of the label.
   */
  std::string get() const;

 private:
  BetterOverlapPtr overlap_;
  Direction direction_;
  bool reverse_complemented_;

  std::string getRawLabel() const;
};

};  // namespace layout

#endif  // LAYOUT_LABEL_H_
