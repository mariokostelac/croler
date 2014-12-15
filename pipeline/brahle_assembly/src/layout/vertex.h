// Copyright 2014 Bruno Rahle

#ifndef LAYOUT_VERTEX_H_
#define LAYOUT_VERTEX_H_

#include <overlap/read.h>

#include <layout/better_read.h>
#include <layout/label.h>

#include <string>
#include <vector>

namespace layout {

class Graph;
class Edge;

typedef layout::Label::Direction DIR;

/**
 * A vertex in a string graph.
 */
class Vertex {
  typedef overlap::String* StringPtr;
  StringPtr data_;
  uint32_t id_;
  const Graph& graph_;
  bool marked_;

  // edges with label direction FROM_ONE_TO_TWO
  std::vector <std::shared_ptr< Edge > > edges_dir1_;
  // edges with label direction FROM_TWO_TO_ONE
  std::vector <std::shared_ptr< Edge > > edges_dir2_;

 public:
  /**
   * Constructor. Given a read, creates a vertex.
   */
  Vertex(BetterRead*, uint32_t, const Graph&);

  /**
   * Destructor.
   */
  virtual ~Vertex();

  /**
   * Getter for the read data this vertex is representing.
   */
  StringPtr& data() { return data_; }

  /**
   * Getter for the read data this vertex is representing (const version).
   */
  const StringPtr& data() const { return data_; }

  /**
   * Getter for the vertex ID.
   */
  uint32_t& id() { return id_; }

  /**
   * Getter for the vertex ID (const version).
   */
  const uint32_t& id() const { return id_; }

  /**
   * User-friendly name of the vertex. Usually just ID.
   */
  virtual const std::string getName() const;

  /**
    * Adds edge to edge vector in given direction
    */
  void AddEdge(std::shared_ptr< Edge > edge, DIR dir);

  /**
   * Returns number of edges of vertex with direction Label::Direction::FROM_ONE_TO_TWO
   */
  const uint32_t count_edges_dir1() const { return edges_dir1_.size() ;}

  /**
   * Returns number of edges of vertex with direction Label::Direction::FROM_TWO_TO_ONE
   */
  const uint32_t count_edges_dir2() const { return edges_dir2_.size() ;}

  /**
   * Mark vertex for removal
   */
  void mark() {marked_ = true ;}

  /**
   * Check if vertex is marked for removal
   */
  bool isMarked() {return marked_ ;}

  /**
   * Mark correspondig edges for removal
   */
  void markEdges();
};

};  // namespace layout

#endif  // LAYOUT_VERTEX_H_
