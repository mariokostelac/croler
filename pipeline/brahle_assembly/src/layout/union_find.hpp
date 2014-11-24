#ifndef OVERLAP_UNIONFIND_H_
#define OVERLAP_UNIONFIND_H_

namespace layout {

class UnionFind {
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
public:
  UnionFind(int n) {
    data_ = new Node[n];
  }

  ~UnionFind() { delete [] data_; }

  int find(int x);
  bool join(int x, int y);
};

}

#endif

