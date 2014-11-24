#include "layout/contig.h"

namespace layout {

Contig::Contig(BetterReadPtr starting, BetterReadSetPtr read_set) :
    reads_(), read_set_(read_set), valid_(false), alive_(true) {
  reads_.emplace_back(starting);
}

Contig::~Contig() {
}

// TODO(brahle): ovo je spora metoda spajanja
void Contig::Join(BetterOverlapPtr better_overlap, Contig* contig) {
  assert(alive_);
  assert(contig->alive_);
  if (IsLeftOverlap(better_overlap)) {
    if (contig->IsLeftOverlap(better_overlap)) {
      reads_.insert(
          reads_.begin(),
          contig->reads_.rbegin(),
          contig->reads_.rend());
    } else {
      reads_.insert(
          reads_.begin(),
          contig->reads_.begin(),
          contig->reads_.end());
    }
  } else {
    if (contig->IsLeftOverlap(better_overlap)) {
      reads_.insert(reads_.end(), contig->reads_.begin(), contig->reads_.end());
    } else {
      reads_.insert(
          reads_.end(),
          contig->reads_.rbegin(),
          contig->reads_.rend());
    }
  }
  valid_ = true;
  contig->Kill();
}

void Contig::SetValid(bool value = true) {
  valid_ = value;
}

bool Contig::IsUsable() const {
  return valid_ && alive_;
}

bool Contig::IsLeftOverlap(BetterOverlapPtr better_overlap) const {
  auto overlap = better_overlap->overlap();
  auto left_read = reads_.front()->id();
  return overlap->read_one == left_read || overlap->read_two == left_read;
}

void Contig::Kill() {
  alive_ = false;
  reads_.clear();
}

size_t Contig::size() const {
  return reads_.size();
}

ContigSet::ContigSet(BetterReadSet* read_set) : contigs_(read_set->size()) {
  for (size_t i = 0; i < read_set->size(); ++i) {
    contigs_[i] = new Contig((*read_set)[i], read_set);
  }
}

ContigSet::~ContigSet() {
  for (auto contig : contigs_) {
    delete contig;
  }
}

size_t ContigSet::size() const {
  return contigs_.size();
}

ContigSet::ContigPtr& ContigSet::operator[](int i) {
  return contigs_[i];
}

const ContigSet::ContigPtr& ContigSet::operator[](int i) const {
  return contigs_[i];
}

};  // namespace layout
