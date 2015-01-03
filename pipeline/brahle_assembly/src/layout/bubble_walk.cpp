// Copyright 2014 Marko Culinovic
#include <layout/bubble_walk.h>
#include <layout/string_graph.h>

#include <set>
#include <string>

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

std::string BubbleWalk::getSequence() {
    // check if this works
    // add first read data to sequence
    std::string sequence(reinterpret_cast<char *> (
                const_cast<uint8_t *> (first->data()->data())));
    // reverse read data if prefix of first read is part of first overlap
    bool isReverse = !edges.empty() &&
            edges[0]->label().overlap()->Suf(first->id()) == 0;
    if (isReverse) sequence = std::string(sequence.rbegin(), sequence.rend());

    for (auto const& edge: edges) {
        // check if reverseal and complements are handled correctly
        // in class label method get() already implements reverse complement
        // but reversal for prefix and suffix isn't handled so reverse needed
        std::string label = edge->label().get();
        if (isReverse) {
            label = std::string(label.rbegin(), label.rend());
        }
        sequence.append(label);
    }

    if (isReverse) sequence = std::string(sequence.rbegin(), sequence.rend());
    return sequence;
}

bool BubbleWalk::containsRead(uint32_t id) {
    return read_ids->count(id) > 0;
}

};  // namespace layout
