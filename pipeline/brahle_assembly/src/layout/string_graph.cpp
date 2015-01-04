// Copyright 2014 Bruno Rahle

#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <set>

#include "layout/string_graph.h"
#include "lib/edlib/src/edlib.h"

namespace layout {

Edge::Edge(
    BetterOverlap* overlap,
    uint32_t id,
    uint32_t read_id,
    const Graph& graph)
    : id_(id),
      graph_(graph),
      marked_(false),
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

void Graph::deleteMarked() {
  std::vector< std::shared_ptr< Vertex > > vertices_temp;
  std::vector< std::shared_ptr< Edge > > edges_temp;
  std::map< uint32_t, uint32_t > id_to_vertex_map_temp;

  // delete vertices from graph
  size_t num_vertices = vertices_.size();
  for (size_t i = 0, j = 0; i < num_vertices; ++i) {
    if (vertices_[i]->isMarked()) {
      vertices_[i].reset();
    } else {
      vertices_temp.emplace_back(vertices_[i]);
      id_to_vertex_map_temp[vertices_[i]->id()] = j++;
    }
  }

  // delete edges from graph
  size_t num_edges = edges_.size();
  for (size_t i = 0; i < num_edges; ++i) {
    if (edges_[i]->isMarked()) {
      edges_[i].reset();
    } else {
      edges_temp.emplace_back(edges_[i]);
    }
  }

  vertices_.clear();
  vertices_ = vertices_temp;
  edges_.clear();
  edges_ = edges_temp;
  id_to_vertex_map_.clear();
  id_to_vertex_map_ = id_to_vertex_map_temp;
}

void Graph::trim(const uint32_t trimSeqLenThreshold) {
  uint32_t disconnected_ctr = 0;
  uint32_t tips_ctr = 0;

  // fprintf(stderr, "vertices: %d\n", vertices_.size());
  // fprintf(stderr, "map: %d\n", id_to_vertex_map_.size());
  // fprintf(stderr, "edges: %d\n", edges_.size());
  fprintf(stderr, "Trimmming threshold: %d\n", trimSeqLenThreshold);


  for (auto const& vertex: vertices_) {

    // check threshold
    if (vertex->data()->size() > trimSeqLenThreshold)
      continue;

    // check if disconnected
    if (vertex->count_edges_B() == 0 && vertex->count_edges_E() == 0) {
      vertex->mark();
      ++disconnected_ctr;
      continue;
    }

    // check if tip
    if (vertex->count_edges_B() == 0 || vertex->count_edges_E() == 0) {
      vertex->mark();
      vertex->markEdges();
      ++tips_ctr;
    }
  }

  if (!(disconnected_ctr == 0 && tips_ctr == 0)) {
    deleteMarked();
  }

  // fprintf(stderr, "vertices: %d\n", vertices_.size());
  // fprintf(stderr, "map: %d\n", id_to_vertex_map_.size());
  // fprintf(stderr, "edges: %d\n", edges_.size());

  fprintf(stderr, "Removed %d tips and %d disconnected vertices\n",
                   tips_ctr, disconnected_ctr);
}

overlap::ReadSet* Graph::extractReads() {
  read_set_ = new overlap::ReadSet(vertices_.size());
  for (auto const& vertex: vertices_) {
    read_set_->Add((overlap::Read *)vertex->data());
  }
  return read_set_;
}

Unitigging::BetterOverlapSetPtr Graph::extractOverlaps() {
  assert(read_set_ != nullptr);
  overlap_set_ = Unitigging::BetterOverlapSetPtr(
                    new BetterOverlapSet(read_set_));
  // add all overlaps from edges
  for (auto &edge: edges_) {
    bool not_in_set = true;
    size_t  overlaps_size = overlap_set_->size();
    for (size_t i = 0; i < overlaps_size; ++i) {
      auto overlap = (*overlap_set_)[i];
      if (overlap->one()->id() == edge->label().overlap()->one()->id() &&
          overlap->two()->id() == edge->label().overlap()->two()->id()) {
        not_in_set = false;
        break;
      }
    }
    if (not_in_set) {
      overlap_set_->Add(new BetterOverlap(edge->label().overlap()));
    }
  }
  return overlap_set_;
}

void Graph::removeBubbles() {
  fprintf(stderr, "Bubble popping\n");
  for (auto const& vertex: vertices_) {
    // std::shared_ptr< Vertex > vertex = getVertex(269 - 159);
    // fprintf(stderr, "Edges of vertex: %d %d\n", vertex->count_edges_B(), vertex->count_edges_E());
    // skip vertices already marked for removal
    if (vertex->isMarked()) continue; // continue;

    // check overlaps where read part of overlap
    // is prefix (dir = 0) or suffix (dir = 1)
    for (size_t dir = 0; dir < 2; ++dir) {
      const std::vector<std::shared_ptr< Edge >> &edges = (dir == 0) ?
                    vertex->getEdgesB() : vertex->getEdgesE();
      if (edges.size() <= 1) continue;

      // check for every edge if vertex B is marked
      bool skip = false;
      for (auto const& edge: edges) {
        if (edge->B()->isMarked()) {
          skip = true;
          break;
        }
      }
      if (skip) continue;

      std::vector<BubbleWalk> bubble_walks;
      getBubbleWalks(vertex, dir, bubble_walks);

      if (bubble_walks.size() == 0) continue; // continue;
      else fprintf(stderr, "Bubble walks: %d\n", bubble_walks.size());

      uint32_t selected_walk = -1;
      double selected_coverage = 0;
      bool is_transitive = false;  // exits walk with only one edge

      // minimum overlaps at start and end of bubble walks
      uint32_t overlap_start = std::numeric_limits<uint32_t>::max();
      uint32_t overlap_end = std::numeric_limits<uint32_t>::max();

      size_t i = 0;
      for (auto bubble_walk: bubble_walks) {
        if (bubble_walk.Edges().size() <= 1) {
          puts("Transitive bubble");
          is_transitive = true;
          break;
        }

        double curr_coverage = 0;
        for (auto const& walk_edge: bubble_walk.Edges()) {
          curr_coverage += walk_edge->B()->coverage();
        }

        if (curr_coverage > selected_coverage || selected_coverage == 0) {
          selected_walk = i;
          selected_coverage = curr_coverage;
        }
        fprintf(stderr, "Coverage: %d\n", curr_coverage);

        std::shared_ptr< Edge > first = bubble_walk.Edges().front();
        std::shared_ptr< Edge > last = bubble_walk.Edges().back();

        if (first->label().overlap()->Length() < overlap_start) {
          overlap_start = first->label().overlap()->Length();
        }
        if (last->label().overlap()->Length() < overlap_end) {
          overlap_end = last->label().overlap()->Length();
        }

        ++i; 
      }

      // bubble is transitive so it's not valid for removal
      if (is_transitive)
        continue;

      // extract sequences from walks
      std::vector< std::string > bubble_sequences;
      for (auto &bubble_walk: bubble_walks) {
        std::shared_ptr< Vertex > start = bubble_walk.Edges().front()->A();
        std::shared_ptr< Vertex > end = bubble_walk.Edges().back()->B();

        std::string sequence = bubble_walk.getSequence();
        std::string matching_sequence;
        uint32_t start_idx = 0;
        uint32_t end_idx = 0;
        if (dir == 0) {
          // prefix of first read is part of first overlap in bubble walk
          // it means sequence direction is from end to start
          start_idx = end->data()->size() - overlap_end;
          end_idx = sequence.size() - (start->data()->size() - overlap_start);
        } else {
          // suffix of first read is part of first overlap in bubble walk
          // it means sequence direction is from start to end
          start_idx = start->data()->size() - overlap_start;
          end_idx = sequence.size() - (end->data()->size() - overlap_end);
        }

        if (end_idx > start_idx) {
          matching_sequence = sequence.substr(start_idx, end_idx - start_idx);
        }
        bubble_sequences.emplace_back(matching_sequence);
      }
      puts("Sequences generated");

      BubbleWalk& walk = bubble_walks[selected_walk];

      // ubaciti šošića

      std::string selected_sequence = walk.getSequence();
      for (size_t j = 0; j < bubble_walks.size(); ++j) {
        if (j == selected_walk) continue;
        BubbleWalk& curr_walk = bubble_walks[j];
        auto &walk_edges = curr_walk.Edges();
        // puts("Novi walk");
        for (size_t k = 0; k < walk_edges.size() - 1; ++k) {
          uint32_t id = walk_edges[k]->B()->id();
          // fprintf(stderr, "Brid: %d %d\n", walk_edges[k]->A()->id() + 159, walk_edges[k]->B()->id() + 159);
          if (!walk.containsRead(id)) {
            fprintf(stderr, "Mark vertex id: %d\n", id);
            getVertex(id)->mark();
            getVertex(id)->markEdges();
          }
        }
      }

      // count bubbles
    }
  }
}


void Graph::getBubbleWalks(const std::shared_ptr<Vertex>& vertex_root,
                            size_t dir,
                            std::vector<BubbleWalk> &bubble_walks) {
  puts("Get Bubble Walks");
  uint32_t reads_cnt = 0;
  uint64_t distance = 0;

  // puts("BFS");
  // breadth-first search graph
  Node *root = new Node(vertex_root, dir, nullptr, nullptr, 0);
  opened_queue.emplace_back(root);
  ++reads_cnt;
  while (!opened_queue.empty()) {
    if (reads_cnt > MAX_NODES) {
      closed_queue.insert(closed_queue.end(),
                          opened_queue.begin(),
                          opened_queue.end());
      opened_queue.clear();
      break;
    }
    std::deque< Node* > expand_queue;
    while (!opened_queue.empty()) {
      Node *node = opened_queue.front();
      opened_queue.pop_front();
      if (distance > MAX_DISTANCE) {
        closed_queue.emplace_back(node);
      } else {
        // expand current vertex
        uint32_t expansion_reads_cnt = node->expand(expand_queue);
        reads_cnt += expansion_reads_cnt;
        if (expansion_reads_cnt == 0)
          closed_queue.emplace_back(node);
      }
    }
    opened_queue = expand_queue;
    std::shared_ptr<Vertex> end_vertex;
    if (bubbleFound(root, &end_vertex)) {
      puts("Bubble found");
      fprintf(stderr, "End vertex id: %d\n", end_vertex->id() + 159);
      generateBubbleWalks(vertex_root, end_vertex, bubble_walks);
      break;
    }
  }
  opened_queue.clear();
  closed_queue.clear();

  if (bubble_walks.size() <= 1 || bubble_walks.size() > MAX_WALKS) {
    bubble_walks.clear();
    return;
  }

  puts("Check for orientation");

  // check that last vertex has same orientation for all walks
  // it's possible that last direction is not the same because of
  // overlap type EE - innnie
  std::shared_ptr<Edge> last_edge = bubble_walks.front().Edges().back();
  uint32_t last_direction = last_edge->label().overlap()->Suf(last_edge->B()->id());

  std::set<uint32_t> vertices_ids;
  for (size_t i = 0; i < bubble_walks.size(); ++i) {
    last_edge = bubble_walks[i].Edges().back();
    if (last_edge->label().overlap()->Suf(last_edge->B()->id()) != last_direction) {
      bubble_walks.clear();
      return;
    }
    // add all vertices from bubble to set vertices_ids
    std::vector<uint32_t> walk_vertices;
    walk_vertices.emplace_back(vertex_root->id());
    for (auto const& edge: bubble_walks[i].Edges()) {
      walk_vertices.emplace_back(edge->B()->id());
    }
    for (auto id: walk_vertices) {
      vertices_ids.insert(id);
    }
  }

  std::shared_ptr< Vertex > end_vertex = bubble_walks.front().Edges().back()->B();
  // check that all links are only to vertices in set vertices_ids
  for (auto id: vertices_ids) {
    std::vector<std::shared_ptr< Edge >> v_edges;
    if (id == vertex_root->id()) {
      v_edges = getVertex(id)->getEdges(dir);
    } else if (id == end_vertex->id()) {
      v_edges = end_vertex->getEdges(last_direction);
    } else {
      v_edges = getVertex(id)->getEdges();
    }

    for (auto const& edge: v_edges) {
      // fprintf(stderr, "Edge: %d %d\n", edge->A()->id() + 159, edge->B()->id() + 159);
      if (vertices_ids.find(edge->B()->id()) == vertices_ids.end()) {
        bubble_walks.clear();
        puts("All vertices not in bubble vertices set");
        return;
      }
    }
  }
}


bool Graph::bubbleFound(Node* root, std::shared_ptr<Vertex>* end) {
  std::deque< Node* > nodes;
  nodes.insert(nodes.end(), opened_queue.begin(), opened_queue.end());
  nodes.insert(nodes.end(), closed_queue.begin(), closed_queue.end());

  for (auto const& end_node: opened_queue) {
    if (end_node->vertex()->id() == root->vertex()->id())
      continue;

    bool isBubbleEnd = true;
    for (auto const& node: nodes) {
      if (!isEndVertex(end_node->vertex(), node, root)) {
        isBubbleEnd = false;
        break;
      }
    }
    if (isBubbleEnd) {
      *end = end_node->vertex();
      return true;
    }
  }
  end->reset();
  return false;
}

bool Graph::isEndVertex(std::shared_ptr<Vertex> end,
                        Node *node,
                        Node *root) {
  if (node == nullptr) return false;
  if (node->vertex()->id() == end->id() && node != root) return true;
  return isEndVertex(end, node->parent(), root);
}


void Graph::generateBubbleWalks(std::shared_ptr< Vertex > start_vertex,
                                std::shared_ptr< Vertex > end_vertex,
                                std::vector<BubbleWalk> &bubble_walks) {
  std::deque< Node* > nodes;
  nodes.insert(nodes.end(), opened_queue.begin(), opened_queue.end());
  nodes.insert(nodes.end(), closed_queue.begin(), closed_queue.end());

  std::set< Node *> end_nodes;
  for (auto const &node: nodes) {
    Node* end_node = nullptr;
    findEndNode(node, end_vertex, end_node);
    assert(end_node != nullptr);
    end_nodes.insert(end_node);
  }

  nodes.clear();
  nodes.insert(nodes.end(), end_nodes.begin(), end_nodes.end());

  for (Node* node: nodes) {
    std::vector<std::shared_ptr< Edge >> walk_edges;
    while (node->parent() != nullptr) {
      walk_edges.emplace_back(node->edge_from_parent());
      node = node->parent();
    }

    // create walk from walk edges
    BubbleWalk* new_walk = new BubbleWalk(start_vertex);
    for (std::vector<std::shared_ptr< Edge >>::reverse_iterator it = walk_edges.rbegin();
          it != walk_edges.rend(); ++it) {
      new_walk->addEdge(*it);
    }
    bubble_walks.emplace_back(*new_walk);
    delete new_walk;
    new_walk = nullptr;
  }
}

void Graph::findEndNode(Node *node, std::shared_ptr< Vertex > end_vertex,
                        Node*& end_node) {
  if (node == nullptr) {
    end_node = nullptr;
    return;
  }
  if (node->vertex()->id() == end_vertex->id()) {
    end_node = node;
    return;
  }
  return findEndNode(node->parent(), end_vertex, end_node);
}

};  // namespace layout
