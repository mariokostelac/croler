// Copyright 2014 Marko Culinovic
#include "layout/bubble_walk.h"

namespace layout {

BubbleWalk::BubbleWalk(std::shared_ptr< Vertex > fst) {
    first = fst;
}

BubbleWalk::~BubbleWalk() {}

};  // namespace layout
