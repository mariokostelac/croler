// Copyright 2014 Bruno Rahle

#include <string>
#include <vector>

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

void Vertex::eraseEdgeDir1(std::shared_ptr< Edge > pair_edge) {
    size_t i = 0;
    bool found = false;
    std::vector<std::shared_ptr< Edge >> temp_edges;
    for (i = 0; i < edges_dir1_.size(); ++i) {
        // if (edges_dir1_[i]->A()->id() == pair_edge->A()->id() &&
        //    edges_dir1_[i]->B()->id() == pair_edge->B()->id()) {
        if (edges_dir1_[i]->label().overlap() == pair_edge->label().overlap()) {
            found = true;
            // fprintf(stderr, "\tFound\n");
        } else {
            temp_edges.push_back(edges_dir1_[i]);
        }
    }
    edges_dir1_.clear();
    edges_dir1_ = temp_edges;
}

void Vertex::eraseEdgeDir2(std::shared_ptr< Edge > pair_edge) {
    size_t i = 0;
    bool found = false;
    std::vector<std::shared_ptr< Edge >> temp_edges;
    for (i = 0; i < edges_dir2_.size(); ++i) {
        // if (edges_dir2_[i]->A()->id() == pair_edge->A()->id() &&
        //    edges_dir2_[i]->B()->id() == pair_edge->B()->id()) {
        if (edges_dir2_[i]->label().overlap() == pair_edge->label().overlap()) {
            found = true;
            // fprintf(stderr, "\tFound\n");
        } else {
            temp_edges.push_back(edges_dir2_[i]);
        }
    }
    edges_dir2_.clear();
    edges_dir2_ = temp_edges;
}

void Vertex::markEdges() {
    if (edges_dir1_.size() == 0) {
      for (auto const& edge: edges_dir2_) {
        edge->mark();
        // find and mark overlap pair edge
        const auto &edges_dir1_B = edge->B()->getEdgesDir1();
        for (auto const& pair_edge: edges_dir1_B) {
            if (pair_edge->label().overlap() == edge->label().overlap()) {
            // if (pair_edge->B()->id() == this->id()) {
                pair_edge->mark();
                pair_edge->A()->eraseEdgeDir1(pair_edge);
            }
        }
      }
      edges_dir2_.clear();
    } else {
      for (auto const&edge: edges_dir1_) {
        edge->mark();
         // find and mark overlap pair edge
        const auto &edges_dir2_B = edge->B()->getEdgesDir2();
        for (auto const& pair_edge: edges_dir2_B) {
            if (pair_edge->label().overlap() == edge->label().overlap()) {
            // if (pair_edge->B()->id() == this->id()) {
                pair_edge->mark();
                pair_edge->A()->eraseEdgeDir2(pair_edge);
            }
        }
      }
      edges_dir1_.clear();
    }
}

const std::vector<std::shared_ptr< Edge >>& Vertex::getEdges(Label::Direction dir) {
    if (dir == Label::Direction::FROM_ONE_TO_TWO)
        return edges_dir1_;
    else
        return edges_dir2_;
}

};  // namespace layout
