// Copyright 2014 Bruno Rahle

#include <string>

#include "layout/vertex.h"

namespace layout {

Vertex::Vertex(BetterRead* read, uint32_t id, const Graph& graph) :
    data_(read->read()), id_(id), graph_(graph) {
}

Vertex::~Vertex() {
}

const std::string Vertex::getName() const {
  static char tmp[1000];
  snprintf(tmp, sizeof(tmp), "%d", id_);
  return std::string(tmp);
}

};  // namespace layout
