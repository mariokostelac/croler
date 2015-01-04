#include "layout/contig.h"

using std::deque;
using std::reverse;

namespace layout {

Contig::Contig(BetterReadPtr starting, BetterReadSetPtr read_set) :
    reads_(), read_set_(read_set), valid_(false), alive_(true) {

  reads_.emplace_back(starting);
}

Contig::~Contig() {
}

bool Contig::ForwardOriented() {
  const auto& f_read = reads_[0]->read();
  const auto& f_overlap = overlaps_[0]->overlap();
  if (f_read->id() == f_overlap->read_one) {
    //   ------>
    // ------>
    if (f_overlap->a_hang < 0) {
      return false;
    }
  } else if (f_read->id() == f_overlap->read_two) {
    if (f_overlap->type == overlap::Overlap::Type::EB) {
      // ------>
      //   ------>
      if (f_overlap->b_hang > 0) {
        return false;
      }
    } else {
      //   ------>
      // <------
      if (f_overlap->b_hang < 0) {
        return false;
      }
    }
  }
  return true;
}

// TODO(brahle): ovo je spora metoda spajanja
void Contig::Join(BetterOverlapPtr better_overlap, Contig* contig) {
  assert(alive_);
  assert(contig->alive_);
  if (IsLeftOverlap(better_overlap)) {
    if (contig->IsLeftOverlap(better_overlap)) {
      reads_.insert(reads_.begin(), contig->reads_.rbegin(), contig->reads_.rend());
      overlaps_.push_front(better_overlap);
      overlaps_.insert(overlaps_.begin(), contig->overlaps_.rbegin(), contig->overlaps_.rend());
    } else {
      reads_.insert(reads_.begin(), contig->reads_.begin(), contig->reads_.end());
      overlaps_.push_front(better_overlap);
      overlaps_.insert(overlaps_.begin(), contig->overlaps_.begin(), contig->overlaps_.end());
    }
  } else {
    if (contig->IsLeftOverlap(better_overlap)) {
      reads_.insert(reads_.end(), contig->reads_.begin(), contig->reads_.end());
      overlaps_.push_back(better_overlap);
      overlaps_.insert(overlaps_.end(), contig->overlaps_.begin(), contig->overlaps_.end());
    } else {
      reads_.insert(reads_.end(), contig->reads_.rbegin(), contig->reads_.rend());
      overlaps_.push_back(better_overlap);
      overlaps_.insert(overlaps_.end(), contig->overlaps_.rbegin(), contig->overlaps_.rend());
    }
  }
  valid_ = true;
  contig->Kill();

  assert(reads_.size() == overlaps_.size() + 1);
}

void Contig::SetValid(bool value = true) {
  valid_ = value;
}

bool Contig::IsUsable() const {
  return valid_ && alive_;
}

//bool Contig::IsLeftOverlap(BetterOverlapPtr better_overlap) const {
  //const auto& overlap = better_overlap->overlap();
  //const auto& contig_start = reads_.front()->id();
  //const auto& overlap_end = better_overlap->GoesFromFirst() ? overlap->read_two : overlap->read_one;
  //return overlap_end == contig_start;
//}

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

const deque< BetterRead* >& Contig::getReads() {
  return reads_;
}

const deque< BetterOverlapPtr >& Contig::getOverlaps() {
  return overlaps_;
}

ContigSet::ContigSet(BetterReadSet* read_set) : contigs_(read_set->size()) {
  for (size_t i = 0; i < read_set->size(); ++i) {
    contigs_[i] = new Contig((*read_set)[i], read_set);
    // make contigs of only usable reads
    if (!((*read_set)[i]->read()->isUsable())) {
      contigs_[i]->Kill();
    } else {
      // puts("usable");
    }
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
