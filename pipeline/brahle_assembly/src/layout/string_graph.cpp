// Copyright 2014 Bruno Rahle

#include <algorithm>
#include <string>

#include "layout/string_graph.h"

namespace layout {

Edge::Edge(
    BetterOverlap* overlap,
    uint32_t id,
    uint32_t read_id,
    const Graph& graph)
    : id_(id),
      graph_(graph),
      A_(graph.getVertex(read_id)),
      B_(graph.getVertex(overlap->Other(read_id))),
      label_(
          overlap,
          overlap->overlap()->read_one == read_id ?
          layout::Label::Direction::FROM_ONE_TO_TWO :
          layout::Label::Direction::FROM_TWO_TO_ONE) {
}

Edge::~Edge() {
}

const std::string Edge::getName() const {
  return label_.get();
}

const std::string Edge::getFormatedName() const {
  std::string ret = getName();
  if (ret.size() > 15) {
    return ret.substr(0, 4) + "(...)" + ret.substr(ret.size()-3);
  }
  return ret;
}

Graph::Graph() : vertices_(), edges_(), finalized_(false) {
}

Graph::~Graph() {
}

bool Graph::EdgeComparator::operator()(
    const std::shared_ptr< Edge > &x,
    const std::shared_ptr< Edge > &y) {
  if (x->A()->id() != y->A()->id()) return x->A()->id() < y->A()->id();
  if (x->B()->id() != y->B()->id()) return x->B()->id() < y->B()->id();
  return false;
}

void Graph::finalize() {
  std::sort(edges_.begin(), edges_.end(),  EdgeComparator());
  finalized_ = true;
}

std::shared_ptr< Vertex > Graph::getVertex(uint32_t id) const {
  return vertices_[id_to_vertex_map_.at(id)];
}

void Graph::printToGraphviz(FILE* file) const {
  fprintf(file, "digraph G {\n");
  for (auto edge : edges_) {
    fprintf(
        file,
        "\"%s\" -> \"%s\" [ label = \"%s\" ];\n",
        edge->A()->getName().c_str(),
        edge->B()->getName().c_str(),
        edge->getFormatedName().c_str());
  }
  fprintf(file, "};\n");
}

};  // namespace layout
