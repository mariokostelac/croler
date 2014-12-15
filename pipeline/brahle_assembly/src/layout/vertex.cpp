// Copyright 2014 Bruno Rahle

#include <string>

#include "layout/vertex.h"
#include "layout/string_graph.h"

namespace layout {

Vertex::Vertex(BetterRead* read, uint32_t id, const Graph& graph) :
    data_(read->read()), id_(id), graph_(graph), marked_(false) {
        edges_dir1_.reserve(10);
        edges_dir2_.reserve(10);
}

Vertex::~Vertex() {
}

const std::string Vertex::getName() const {
  static char tmp[1000];
  snprintf(tmp, sizeof(tmp), "%d", id_);
  return std::string(tmp);
}

void Vertex::AddEdge(std::shared_ptr< Edge > edge, DIR dir) {
  if (dir == DIR::FROM_ONE_TO_TWO) {
    edges_dir1_.emplace_back(edge);
  } else {
    edges_dir2_.emplace_back(edge);
  }
}

void Vertex::markEdges() {
    if (edges_dir1_.size() == 0) {
      for (auto const& edge: edges_dir2_) {
        edge->mark();
      }
      edges_dir2_.clear();
    } else {
      for (auto const&edge: edges_dir1_) {
        edge->mark();
      }
      edges_dir1_.clear();
    }
}

};  // namespace layout
