#include <algorithm>

#include "overlap.h"


namespace overlap {


Overlap::Overlap(
    uint32_t r1,
    uint32_t r2,
    uint32_t l1,
    uint32_t l2,
    Type t,
    int32_t s)
    : read_one(r1),
      read_two(r2),
      len_one(l1),
      len_two(l2),
      type(t),
      score(s) {}

bool Overlap::operator<(const Overlap& rhs) const {
  if (read_one != rhs.read_one) return read_one < rhs.read_one;
  if (read_two != rhs.read_two) return read_two < rhs.read_two;
  if (type != rhs.type) return type < rhs.type;
  if (len_one != rhs.len_one) return len_one < rhs.len_one;
  if (len_two != rhs.len_two) return len_two < rhs.len_two;
  return score > rhs.score;
}

bool OverlapCmp::operator()(const Overlap* lhs, const Overlap* rhs) const {
  return *lhs < *rhs;
}

OverlapSet::OverlapSet(size_t capacity) {
  overlaps_.reserve(capacity);
}

OverlapSet::~OverlapSet() {
  for (Overlap* o : overlaps_) {
    delete o;
  }
}

void OverlapSet::Add(Overlap* overlap) {
 overlaps_.push_back(overlap);
}

Overlap* OverlapSet::Get(uint32_t idx) const {
  return overlaps_[idx];
}

void OverlapSet::Sort() {
  std::sort(overlaps_.begin(), overlaps_.end(), OverlapCmp());
}

size_t OverlapSet::size() const {
  return overlaps_.size();
}

Overlap* const& OverlapSet::operator[](uint32_t idx) const {
  return overlaps_[idx];
}

}  // namespace overlap
