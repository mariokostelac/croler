// Copyright 2014 Marko Culinovic
#ifndef LAYOUT_BUBBLE_WALK_H
#define LAYOUT_BUBBLE_WALK_H

#include <layout/vertex.h>

#include <memory>
#include <vector>
#include <set>

namespace layout {

class Edge;

class BubbleWalk {
    public:
        explicit BubbleWalk(std::shared_ptr< Vertex > first);
        ~BubbleWalk();
        BubbleWalk(const BubbleWalk& other);
        BubbleWalk& operator=(const BubbleWalk& other);

        void addEdge(std::shared_ptr<Edge> edge);
        std::vector<std::shared_ptr< Edge >>& Edges() { return edges; }
        std::string getSequence();
    private:
        std::shared_ptr< Vertex > first;
        std::vector<std::shared_ptr< Edge >> edges;
        std::set<uint32_t> *read_ids;
};

};  // namespace layout

#endif
