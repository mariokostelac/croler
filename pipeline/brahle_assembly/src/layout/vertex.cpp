// Copyright 2014 Bruno Rahle

#include <string>

#include "layout/vertex.h"

namespace layout {

Vertex::Vertex(BetterRead* read, uint32_t id, const Graph& graph) :
    data_(read->read()), id_(id), graph_(graph) {
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

};  // namespace layout
