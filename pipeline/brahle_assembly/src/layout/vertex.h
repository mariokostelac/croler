// Copyright 2014 Bruno Rahle

#ifndef LAYOUT_VERTEX_H_
#define LAYOUT_VERTEX_H_

#include <overlap/read.h>

#include <layout/better_read.h>

#include <string>

namespace layout {

class Graph;

/**
 * A vertex in a string graph.
 */
class Vertex {
  typedef overlap::String* StringPtr;
  StringPtr data_;
  uint32_t id_;
  const Graph& graph_;

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
};

};  // namespace layout

#endif  // LAYOUT_VERTEX_H_
