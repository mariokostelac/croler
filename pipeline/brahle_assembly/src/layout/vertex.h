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

  // overlaps where read is prefix
  std::vector <std::shared_ptr< Edge > > edges_B_;
  // overlaps where read is suffix
  std::vector <std::shared_ptr< Edge > > edges_E_;

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
   * Adds edge to vertex vector in accordance with read
   * being prefix or suffix of overlap
   */
  void AddEdge(std::shared_ptr< Edge > edge);

  /**
    * Adds edge to edge vector in given direction
    */
  void AddEdge(std::shared_ptr< Edge > edge, DIR dir);

  /**
   * Return edges with given direction - read is prefix or suffix in overlap
   */
  const std::vector<std::shared_ptr< Edge >>& getEdges(uint32_t dir);

  /**
   * Returns edges of vertex with direction Label::Direction::FROM_ONE_TO_TWO
   */
  std::vector<std::shared_ptr< Edge >>& getEdgesDir1() { return edges_dir1_ ;}

  /**
   * Returns edges of vertex with direction Label::Direction::FROM_TWO_TO_ONE
   */
  std::vector<std::shared_ptr< Edge >>& getEdgesDir2() { return edges_dir2_ ;}

  /**
   * Returns number of edges of vertex with direction Label::Direction::FROM_ONE_TO_TWO
   */
  const uint32_t count_edges_dir1() const { return edges_dir1_.size() ;}

  /**
   * Returns number of edges of vertex with direction Label::Direction::FROM_TWO_TO_ONE
   */
  const uint32_t count_edges_dir2() const { return edges_dir2_.size() ;}


  /**
   * Return all edges of vertex
   */
  std::vector<std::shared_ptr< Edge >> getEdges();

    /**
   * Returns edges of vertex where read of overlap is prefix
   */
  std::vector<std::shared_ptr< Edge >>& getEdgesB() { return edges_B_ ;}

  /**
   * Returns edges of vertex where read of overlap is suffix
   */
  std::vector<std::shared_ptr< Edge >>& getEdgesE() { return edges_E_ ;}

  /**
   * Number of overlaps where read is prefix
   */
  const uint32_t count_edges_B() const { return edges_B_.size() ;}

  /**
   * Number of overlaps where read is suffix
   */
  const uint32_t count_edges_E() const { return edges_E_.size() ;}

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

  /**
   * Erase given edge in direction one for removal 
   */
  void eraseEdgeDir1(std::shared_ptr< Edge > pair_edge);

  /**
   * Erase given edge in direction two for removal 
   */
  void eraseEdgeDir2(std::shared_ptr< Edge > pair_edge);

  /**
   * Erase given edge in edges_B
   */
  void eraseEdgeB(std::shared_ptr< Edge > pair_edge);

  /**
   * Erase given edge in edges_E 
   */
  void eraseEdgeE(std::shared_ptr< Edge > pair_edge);

  /**
   * Read coverage
   */
   double coverage() { return ((overlap::Read *) data_)->coverage(); }
};

};  // namespace layout

#endif  // LAYOUT_VERTEX_H_
