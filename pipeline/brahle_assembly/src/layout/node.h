// Copyright 2014 Marko Culinovic

#ifndef LAYOUT_NODE_H
#define LAYOUT_NODE_H

#include <layout/vertex.h>

#include <deque>
#include <memory>

namespace layout {

class Edge;

// BFS Wrapper for Vertex
class Node {
  public:
    Node(std::shared_ptr<Vertex> vertex, uint32_t dir, Node *parent,
         std::shared_ptr<Edge> edge_from_parent, uint32_t distance);
    ~Node();
    uint32_t expand(std::deque<Node*>& expand_queue);

    const std::shared_ptr< Vertex >& vertex() const { return vertex_; }

  private:
    std::shared_ptr<Vertex> vertex_;
    uint32_t direction_;
    Node *parent_;
    std::shared_ptr<Edge> edge_from_parent_;
    uint32_t num_children_;
    uint64_t distance_;
};

};  // namespace layout

#endif  // LAYOUT_NODE_H
