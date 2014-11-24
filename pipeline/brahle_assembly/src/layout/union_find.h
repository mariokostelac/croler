// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_UNION_FIND_H_
#define LAYOUT_UNION_FIND_H_

namespace layout {

/**
 * An implementation of the disjoint-set datastructure.
 */
class UnionFind {
 public:
  /**
   * Creates a UnionFind object that stores n objects.
   */
  explicit UnionFind(int n);

  /**
   * Destructor.
   */
  ~UnionFind();

  /**
   * Finds the group element x belongs to.
   */
  int find(int x);

  /**
   * Joins groups elements x and y belong to. Returns the new group both
   * belong to after the operation, which is guarnateed to be either x or y.
   */
  int join(int x, int y);

 private:
  struct Node {
    int value_;
    int count_;
    int parent_;
    Node() : Node(0, 1, -1) {}
    Node(int value, int parent) : Node(value, 1, parent) {}
    Node(int value, int count, int parent) :
      value_(value), count_(count), parent_(parent) {}
  };

  Node *data_;
};

}  // namespace layout

#endif  // LAYOUT_UNION_FIND_H_

