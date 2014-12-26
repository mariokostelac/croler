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
        edges_B_.reserve(10);
        edges_E_.reserve(10);
}

Vertex::~Vertex() {
}

const std::string Vertex::getName() const {
  static char tmp[1000];
  snprintf(tmp, sizeof(tmp), "%d", id_);
  return std::string(tmp);
}

void Vertex::AddEdge(std::shared_ptr< Edge > edge) {
    if (edge->label().overlap()->Suf(this->id()) == 0) {
        edges_B_.emplace_back(edge);
    } else {
        edges_E_.emplace_back(edge);
    }
}

void Vertex::AddEdge(std::shared_ptr< Edge > edge, DIR dir) {
  if (dir == DIR::FROM_ONE_TO_TWO) {
    edges_dir1_.emplace_back(edge);
  } else {
    edges_dir2_.emplace_back(edge);
  }
}

void Vertex::eraseEdgeB(std::shared_ptr< Edge > edge) {
    size_t i = 0;
    bool found = false;
    std::vector<std::shared_ptr< Edge >> temp_edges;
    for (i = 0; i < edges_B_.size(); ++i) {
        if (edges_B_[i]->A()->id() == edge->B()->id() && edges_B_[i]->B()->id() == edge->A()->id()) {
        // if (edges_B_[i]->label().overlap() == edge->label().overlap()) {
            found = true;
        } else {
            temp_edges.push_back(edges_B_[i]);
        }
    }
    edges_B_.clear();
    edges_B_ = temp_edges;
}

void Vertex::eraseEdgeE(std::shared_ptr< Edge > edge) {
    size_t i = 0;
    bool found = false;
    std::vector<std::shared_ptr< Edge >> temp_edges;
    for (i = 0; i < edges_E_.size(); ++i) {
        if (edges_E_[i]->A()->id() == edge->B()->id() && edges_E_[i]->B()->id() == edge->A()->id()) {
        // if (edges_E_[i]->label().overlap() == edge->label().overlap()) {
            found = true;
        } else {
            temp_edges.push_back(edges_E_[i]);
        }
    }
    edges_E_.clear();
    edges_E_ = temp_edges;
}

void Vertex::eraseEdgeDir1(std::shared_ptr< Edge > edge) {
    size_t i = 0;
    bool found = false;
    std::vector<std::shared_ptr< Edge >> temp_edges;
    for (i = 0; i < edges_dir1_.size(); ++i) {
        if (edges_dir1_[i]->A()->id() == edge->B()->id() && edges_dir1_[i]->B()->id() == edge->A()->id()) {
        // if (edges_dir1_[i]->label().overlap() == edge->label().overlap()) {
            found = true;
        } else {
            temp_edges.push_back(edges_dir1_[i]);
        }
    }
    edges_dir1_.clear();
    edges_dir1_ = temp_edges;
    if (edge->label().overlap()->Suf(this->id()) == 0) {
        eraseEdgeB(edge);
    } else {
        eraseEdgeE(edge);
    }
}

void Vertex::eraseEdgeDir2(std::shared_ptr< Edge > edge) {
    size_t i = 0;
    bool found = false;
    std::vector<std::shared_ptr< Edge >> temp_edges;
    for (i = 0; i < edges_dir2_.size(); ++i) {
        if (edges_dir2_[i]->A()->id() == edge->B()->id() && edges_dir2_[i]->B()->id() == edge->A()->id()) {
        // if (edges_dir2_[i]->label().overlap() == edge->label().overlap()) {
            found = true;
        } else {
            temp_edges.push_back(edges_dir2_[i]);
        }
    }
    edges_dir2_.clear();
    edges_dir2_ = temp_edges;
    if (edge->label().overlap()->Suf(this->id()) == 0) {
        eraseEdgeB(edge);
    } else {
        eraseEdgeE(edge);
    }
}

void Vertex::markEdges() {
    if (edges_B_.size() == 0) {
      for (auto const& edge: edges_E_) {
        // fprintf(stderr, "Prvi: %d %d\n", edge->A()->id(), edge->B()->id());
        edge->mark();
        // find and mark overlap pair edge
        uint32_t dir = 1;
        auto &edges_opposite = edge->B()->getEdgesDir2();
        if (edge->label().direction() == Label::Direction::FROM_TWO_TO_ONE) {
            dir = 2;
            edges_opposite = edge->B()->getEdgesDir1();
        }
        for (auto const& pair_edge: edges_opposite) {
            if (pair_edge->A()->id() == edge->B()->id() && pair_edge->B()->id() == edge->A()->id()) {
            // if (pair_edge->label().overlap() == edge->label().overlap()) {
                // fprintf(stderr, "Par: %d %d\n", pair_edge->A()->id(), pair_edge->B()->id());
                pair_edge->mark();
                if (dir == 1) {
                    pair_edge->A()->eraseEdgeDir1(pair_edge);
                    edge->A()->eraseEdgeDir2(edge);
                } else {
                    pair_edge->A()->eraseEdgeDir2(pair_edge);
                    edge->A()->eraseEdgeDir1(edge);
                }
            }
        }
      }
      edges_E_.clear();
    } else {
      for (auto const&edge: edges_B_) {
        // fprintf(stderr, "Prvi: %d %d\n", edge->A()->id(), edge->B()->id());
        edge->mark();
         // find and mark overlap pair edge
        uint32_t dir = 1;
        auto &edges_opposite = edge->B()->getEdgesDir2();
        if (edge->label().direction() == Label::Direction::FROM_TWO_TO_ONE) {
            dir = 2;
            edges_opposite = edge->B()->getEdgesDir1();
        }
        for (auto const& pair_edge: edges_opposite) {
            if (pair_edge->A()->id() == edge->B()->id() && pair_edge->B()->id() == edge->A()->id()) {
            // if (pair_edge->label().overlap() == edge->label().overlap()) {
                // fprintf(stderr, "Par: %d %d\n", pair_edge->A()->id(), pair_edge->B()->id());
                pair_edge->mark();
                if (dir == 1) {
                    pair_edge->A()->eraseEdgeDir1(pair_edge);
                    edge->A()->eraseEdgeDir2(edge);
                } else {
                    pair_edge->A()->eraseEdgeDir2(pair_edge);
                    edge->A()->eraseEdgeDir1(edge);
                }
            }
        }
      }
      edges_B_.clear();
    }
}

std::vector<std::shared_ptr< Edge >> Vertex::getEdges() {
    std::vector<std::shared_ptr< Edge >> all_edges;
    all_edges.insert(all_edges.end(), edges_B_.begin(), edges_B_.end());
    all_edges.insert(all_edges.end(), edges_E_.begin(), edges_E_.end());
    return all_edges;
}

const std::vector<std::shared_ptr< Edge >>& Vertex::getEdges(uint32_t dir) {
    if (dir == 0)
        return edges_B_;
    else
        return edges_E_;
}

};  // namespace layout
