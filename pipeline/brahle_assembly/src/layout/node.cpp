// Copyright 2014 Marko Culinovic


#include <layout/node.h>
#include <layout/string_graph.h>

#include <deque>
#include <vector>

namespace layout {

Node::Node(std::shared_ptr<Vertex> vertex, uint32_t dir, Node *parent,
         std::shared_ptr<Edge> edge_from_parent, uint32_t distance) {
    vertex_ = vertex;
    direction_ = dir;
    parent_ = parent;
    edge_from_parent_ = edge_from_parent;
    if (parent == nullptr) {
        distance_ = 0;
    } else {
        distance_ = parent_->distance_ + distance;
    }
    num_children_ = 0;
}

Node::~Node() {
    assert(num_children_ == 0);
    if (parent_ != NULL) parent_->num_children_--;
}

uint32_t Node::expand(std::deque<Node*>& expand_queue) {
    assert(num_children_ == 0);

    // fprintf(stderr, "Id: %d, Dir: %d\n", vertex_->id() + 159, direction_);
    std::vector<std::shared_ptr< Edge >> edges = vertex_->getEdges(direction_);
    for (auto const &edge: edges) {
        uint32_t child_expand_dir = direction_;  // !direction_;
        Node *child = new Node(edge->B(), child_expand_dir, this,
                              edge, edge->label().get().length());
        expand_queue.emplace_back(child);
        // fprintf(stderr, "\tchild id: %d, dir: %d\n", (child->vertex()->id() + 159), child_expand_dir);
        ++num_children_;
    }
    return edges.size();
}

};  // namespace layout
