// Copyright 2014 Bruno Rahle
#include "layout/better_overlap.h"

namespace layout {

BetterOverlap::BetterOverlap() : BetterOverlap(nullptr, nullptr) {
}

BetterOverlap::BetterOverlap(
    OverlapPtr overlap,
    ReadSetPtr read_set) :
    overlap_(overlap),
    read_set_(read_set) {
}

BetterOverlap::BetterOverlap(BetterOverlap* other) :
    overlap_(other->overlap_),
    read_set_(other->read_set_) {
}

BetterOverlap::~BetterOverlap() {
}

const overlap::Read* BetterOverlap::one() const {
  return read_set_->Get(overlap_->read_one);
}

const overlap::Read* BetterOverlap::two() const {
  return read_set_->Get(overlap_->read_two);
}

const overlap::Read* BetterOverlap::get(uint32_t read) const {
  return read_set_->Get(read);
}

const uint32_t BetterOverlap::Hang(uint32_t read) const {
  assert(overlap_ != nullptr);
  assert(read_set_ != nullptr);
  assert(read == overlap_->read_one || read == overlap_->read_two);
  if (read == overlap_->read_one) {
    return read_set_->Get(overlap_->read_one)->size() - overlap_->len_one;
  }
  return read_set_->Get(overlap_->read_two)->size() - overlap_->len_two;
}

const uint32_t BetterOverlap::Length() const {
  assert(overlap_ != nullptr);
  return (overlap_->len_one + overlap_->len_two) / 2;
}

const uint32_t BetterOverlap::Length(uint32_t read) const {
  assert(overlap_ != nullptr);
  assert(read == overlap_->read_one || read == overlap_->read_two);
  if (read == overlap_->read_one) {
    return overlap_->len_one;
  }
  return overlap_->len_two;
}

const uint32_t BetterOverlap::Suf(uint32_t read) const {
  assert(overlap_ != nullptr);
  assert(read == overlap_->read_one || read == overlap_->read_two);
  if (read == overlap_->read_one) {
    return
        overlap_->type == overlap::Overlap::Type::EB ||
        overlap_->type == overlap::Overlap::Type::EE;
  }
  return
      overlap_->type == overlap::Overlap::Type::BE ||
      overlap_->type == overlap::Overlap::Type::EE;
}

const uint32_t BetterOverlap::Other(uint32_t read) const {
  assert(overlap_ != nullptr);
  assert(read == overlap_->read_one || read == overlap_->read_two);
  if (read == overlap_->read_one) {
    return overlap_->read_two;
  }
  return overlap_->read_one;
}

BetterOverlapSetIter::BetterOverlapSetIter(
    const BetterOverlapSet* better_overlap_set,
    int position) :
    better_overlap_set_(better_overlap_set),
    position_(position) {
}

bool BetterOverlapSetIter::operator!=(const BetterOverlapSetIter& other) const {
  return
      better_overlap_set_ != other.better_overlap_set_ ||
      position_ != other.position_;
}

bool BetterOverlapSetIter::operator==(const BetterOverlapSetIter& other) const {
  return
      better_overlap_set_ == other.better_overlap_set_ &&
      position_ == other.position_;
}

BetterOverlapPtr BetterOverlapSetIter::operator*() const {
  return (*better_overlap_set_)[position_];
}

const BetterOverlapSetIter& BetterOverlapSetIter::operator++() {
  position_++;
  return *this;
}

const BetterOverlapSetIter& BetterOverlapSetIter::operator++(int ignorable) {
  position_++;
  return *this;
}

BetterOverlapSet::BetterOverlapSet(overlap::ReadSet* read_set) :
    read_set_(read_set),
    overlap_set_() {
}

BetterOverlapSet::BetterOverlapSet(
    overlap::ReadSet* read_set,
    overlap::OverlapSet* overlap_set) :
    read_set_(read_set),
    overlap_set_() {
  for (size_t i = 0; i < overlap_set->size(); ++i) {
    overlap_set_.push_back(new BetterOverlap((*overlap_set)[i], read_set_));
  }
}

BetterOverlapSet::~BetterOverlapSet() {
  for (size_t i = 0; i < overlap_set_.size(); ++i) {
    delete overlap_set_[i];
  }
  overlap_set_.clear();
}

BetterOverlapPtr& BetterOverlapSet::operator[](size_t i) {
  return overlap_set_[i];
}

const BetterOverlapPtr& BetterOverlapSet::operator[](size_t i) const {
  return overlap_set_[i];
}

const size_t BetterOverlapSet::size() const {
  return overlap_set_.size();
}

BetterOverlapSetIter BetterOverlapSet::begin() const {
  return BetterOverlapSetIter(this, 0);
}

BetterOverlapSetIter BetterOverlapSet::end() const {
  return BetterOverlapSetIter(this, size());
}

void BetterOverlapSet::Add(BetterOverlapPtr overlap) {
  overlap_set_.emplace_back(overlap);
}

void BetterOverlapSet::Add(overlap::Overlap* overlap) {
  overlap_set_.emplace_back(new BetterOverlap(overlap, read_set_));
}

};  // namespace layout
