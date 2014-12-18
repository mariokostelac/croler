// Copyright 2014 Marko Culinovic
#ifndef LAYOUT_BUBBLE_WALK_H
#define LAYOUT_BUBBLE_WALK_H

#include <memory>
#include <vector>
#include <set>

#include "layout/string_graph.h"
#include "layout/vertex.h"

namespace layout {

class BubbleWalk {
    public:
        explicit BubbleWalk(std::shared_ptr< Vertex > first);
        ~BubbleWalk();
    private:
        std::shared_ptr< Vertex > first;
        std::vector<std::shared_ptr< Edge >> edges;
        std::set<uint32_t> read_ids;
};

};  // namespace layout

#endif
