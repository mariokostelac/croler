// Copyright 2014 Bruno Rahle

#include <overlap/overlap.h>

#include <algorithm>
#include <string>

#include "layout/label.h"

namespace layout {

Label::Label(BetterOverlapPtr overlap, Direction direction) :
    overlap_(overlap),
    direction_(direction),
    reverse_complemented_(
        overlap->overlap()->type != overlap::Overlap::Type::EB) {
}

Label::~Label() {
}

std::string Label::getRawLabel() const {
  uint32_t idx;
  if (direction_ == FROM_ONE_TO_TWO) {
    idx = overlap_->overlap()->read_two;
  } else {
    idx = overlap_->overlap()->read_one;
  }
  auto read = overlap_->get(idx);
  auto data = read->data();
  int len, start, end;

  if (overlap_->Suf(idx)) {
    len = overlap_->Hang(idx);
    start = 0;
    end = len;
  } else {
    len = read->size() - overlap_->Hang(idx) + 1;
    start = overlap_->Hang(idx);
    end = read->size() + 1;
  }

  std::string ret;
  ret.reserve(len);
  for (int i = start; i < end; ++i) {
    ret.push_back(data[i]);
  }
  return ret;
}

std::string Label::get() const {
  std::string ret = getRawLabel();
  if (reverse_complemented_) {
    std::reverse(ret.begin(), ret.end());
    for (auto &c : ret) {
      switch (c) {
        case 'A': c = 'T'; break;
        case 'C': c = 'G'; break;
        case 'G': c = 'C'; break;
        case 'T': c = 'A'; break;
      }
    }
  }
  return ret;
}

}  // namespace layout


