// Copyright 2014 Bruno Rahle
#ifndef LAYOUT_STRING_GRAPH_H_
#define LAYOUT_STRING_GRAPH_H_

#include <overlap/read.h>

#include <layout/better_read.h>
#include <layout/label.h>
#include <layout/unitigging.h>
#include <layout/vertex.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace layout {

class Graph;

/**
 * An edge in a string graph.
 */
class Edge {
 public:
  /**
   * An edge can have two types:
   *    - NORMAL -> direction is the same as the overlap it represents
   *    - REVERSED -> direction is opposite to the one the overlap represents
   */
  enum Type {
    NORMAL,
    REVERSED
  };

  /**
   * Creates an edge from the given overlap with the first read being the one
   * with ID equal to read_id.
   */
  Edge(
      BetterOverlap* overlap,
      uint32_t id,
      uint32_t read_id,
      const Graph& graph);

  /**
   * Destructor.
   */
  virtual ~Edge();

  /**
   * Getter for the vertex A of the edge.
   */
  std::shared_ptr< Vertex >& A() { return A_; }

  /**
   * Getter for the vertex A of the edge (const version).
   */
  const std::shared_ptr< Vertex >& A() const { return A_; }

  /**
   * Getter for the vertex B of the edge.
   */
  std::shared_ptr< Vertex >& B() { return B_; }

  /**
   * Getter for the vertex B of the edge (const version).
   */
  const std::shared_ptr< Vertex >& B() const { return B_; }

  /**
   * Getter for the label (const version).
   */
  const Label& label() const { return label_; }

  /**
   * Full user-friendly name of the edge. Usually just the label.
   */
  virtual const std::string getName() const;

  /**
   * Reformated user-friendly name of the edge. Has a "finite" number of
   * characters.
   */
  virtual const std::string getFormatedName() const;

 private:
  uint32_t id_;
  const Graph& graph_;
  std::shared_ptr< Vertex > A_, B_;
  Label label_;
  Type type_;
};

/**
 * String Graph. Use Graph::create() to create the graph.
 */
class Graph {
  /**
   * Comparator for two edges.
   *
   * Return true if (x.A.ID, x.B.ID) is before (y.A.ID, y.B.ID).
   */
  class EdgeComparator {
   public:
    bool operator()(
        const std::shared_ptr< Edge > &x,
        const std::shared_ptr< Edge > &y);
  };

  /**
   * Used to store some information about vertices.
   */
  struct VertexInfo {
   public:
    uint32_t order_;
    uint32_t first_edge_, last_edge_;
  };

  std::vector< std::shared_ptr< Vertex > > vertices_;
  std::vector< std::shared_ptr< Edge > > edges_;
  std::map< uint32_t, uint32_t > id_to_vertex_map_;
  bool finalized_;

  /**
   * Default constructor is private (by design). Use Graph::create() instead.
   */
  Graph();

 public:
  /**
   * Destructor.
   */
  virtual ~Graph();

  /**
   * Finishes construction of the graph.
   */
  void finalize();

  /**
   * Returns a Vertex with the given ID.
   */
  std::shared_ptr< Vertex > getVertex(uint32_t) const;

  /**
   * Create a string graph from a given set of reads and overlaps.
   */
  static Graph create(
      Unitigging::BetterReadSetPtr reads_,
      Unitigging::BetterOverlapSetPtr overlaps_) {
    Graph g;
    g.vertices_.reserve(reads_->size());
    g.edges_.reserve(overlaps_->size());
    for (auto read : *reads_) {
      g.id_to_vertex_map_[read->id()] = g.vertices_.size();
      g.vertices_.push_back(
          std::shared_ptr< Vertex >(new Vertex(read, read->id(), g)));
    }
    for (auto overlap : *overlaps_) {
      g.edges_.push_back(
          std::shared_ptr< Edge >(
              new Edge(
                  overlap,
                  g.edges_.size(),
                  overlap->overlap()->read_one,
                  g)));
      g.edges_.push_back(
          std::shared_ptr< Edge >(
              new Edge(
                  overlap,
                  g.edges_.size(),
                  overlap->overlap()->read_two,
                  g)));
    }
    g.finalize();
    return g;
  }

  /**
   * Prints the graph to the given file.
   */
  void printToGraphviz(FILE* file) const;
};

};  // namespace layout

#endif  // LAYOUT_STRING_GRAPH_H_

