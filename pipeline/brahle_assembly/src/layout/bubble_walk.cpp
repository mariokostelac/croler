// Copyright 2014 Marko Culinovic
#include <layout/bubble_walk.h>
#include <layout/string_graph.h>

#include <set>

namespace layout {

BubbleWalk::BubbleWalk(std::shared_ptr< Vertex > fst) {
    first = fst;
    read_ids = new std::set<uint32_t>;
}

BubbleWalk::BubbleWalk(const BubbleWalk& other) {
    first = other.first;
    edges = other.edges;
    if (other.read_ids == nullptr) {
        read_ids = nullptr;
    } else {
        read_ids = new std::set<uint32_t>(*other.read_ids);
    }
}

BubbleWalk::~BubbleWalk() {
    if (read_ids != nullptr) {
        delete read_ids;
        read_ids = nullptr;
    }
}

BubbleWalk& BubbleWalk::operator=(const BubbleWalk& other) {
    if (&other == this) return *this;
    first = other.first;
    edges = other.edges;
    if (read_ids != nullptr) delete read_ids;
    if (other.read_ids != nullptr) read_ids = new std::set<uint32_t>(*other.read_ids);
    return *this;
}

void BubbleWalk::addEdge(std::shared_ptr<Edge> edge) {
    edges.emplace_back(edge);
    read_ids->insert(edge->B()->id());
}

};  // namespace layout
