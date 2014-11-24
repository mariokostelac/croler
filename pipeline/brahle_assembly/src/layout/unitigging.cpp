// Copyright 2014 Bruno Rahle

#include <layout/union_find.h>

#include <cstring>

#include <algorithm>
#include <vector>

#include "layout/unitigging.h"

namespace layout {

inline bool eq(double x, double y, double eps) {
  return y <= x + eps && x <= y + eps;
}

Unitigging::Unitigging(
    overlap::ReadSet* reads,
    overlap::OverlapSet* overlaps)
    : reads_(reads),
      orig_overlaps_(overlaps),
      overlaps_(reads, overlaps),
      no_contains_(nullptr),
      no_transitives_(nullptr),
      contigs_(nullptr) {
}

Unitigging::~Unitigging() {
}

void Unitigging::start() {
  removeContainmentEdges();
  removeTransitiveEdges();
  makeContigs();
}

Unitigging::ContigSetPtr& Unitigging::contigs() {
  return contigs_;
}

const Unitigging::BetterOverlapSetPtr& Unitigging::noContains() const {
  return no_contains_;
}

const Unitigging::BetterOverlapSetPtr& Unitigging::noTransitives() const {
  return no_transitives_;
}

const Unitigging::BetterReadSetPtr& Unitigging::readSet() {
  if (better_read_set_ == nullptr) {
    better_read_set_ = BetterReadSetPtr(new BetterReadSet(reads_, false));
  }
  return better_read_set_;
}

void Unitigging::removeContainmentEdges() {
  bool *erased = new bool[reads_->size()];
  memset(erased, 0, sizeof(bool) * reads_->size());
  for (size_t i = 0; i < overlaps_.size(); ++i) {
    BetterOverlap* better_overlap = overlaps_[i];
    overlap::Overlap* overlap = better_overlap->overlap();
    if (better_overlap->one()->size() == overlap->len_one) {
      erased[overlap->read_one] = true;
    } else if (better_overlap->two()->size() == overlap->len_two) {
      erased[overlap->read_two] = true;
    }
  }
  no_contains_ = BetterOverlapSetPtr(new BetterOverlapSet(reads_));
  for (size_t i = 0; i < overlaps_.size(); ++i) {
    BetterOverlap* better_overlap = overlaps_[i];
    overlap::Overlap* overlap = better_overlap->overlap();
    if (erased[overlap->read_one]) continue;
    if (erased[overlap->read_two]) continue;
    no_contains_->Add(new BetterOverlap(better_overlap));
  }
  int erased_count = std::count(erased, erased + reads_->size(), true);
  fprintf(
      stderr,
      "Contained reads = %d (%.2lf%%)\n",
      erased_count,
      static_cast< double >(erased_count * 100.0) / reads_->size());
  fprintf(
      stderr,
      "Edges after removing contained reads = %d (%.2lf%%)\n",
      no_contains_->size(),
      static_cast< double >(no_contains_->size()* 100.0) / overlaps_.size());
  delete [] erased;
}

bool Unitigging::isTransitive(
    BetterOverlap* o1,
    BetterOverlap* o2,
    BetterOverlap* o3) const {
  auto A = o1->overlap()->read_one;
  auto B = o1->overlap()->read_two;
  auto C = o2->overlap()->read_one;
  if (C == A) {
    C = o2->overlap()->read_two;
  }
  if (o2->Suf(C) == o3->Suf(C)) return false;
  if (o1->Suf(A) != o2->Suf(A)) return false;
  if (o1->Suf(B) != o3->Suf(B)) return false;
  if (!eq(
          o2->Hang(A) + o3->Hang(C),
          o1->Hang(A),
          EPSILON * o1->Length() + ALPHA)) {
    return false;
  }
  if (!eq(
          o2->Hang(C) + o3->Hang(B),
          o1->Hang(B),
          EPSILON * o1->Length() + ALPHA)) {
    return false;
  }
  return true;
}

inline void Unitigging::removeTransitiveEdges() {
  layout::BetterReadSet brs(reads_, 1);
  for (size_t i = 0; i < no_contains_->size(); ++i) {
    auto better_overlap = (*no_contains_)[i];
    auto overlap = better_overlap->overlap();
    brs[overlap->read_one]->AddOverlap(better_overlap);
    brs[overlap->read_two]->AddOverlap(better_overlap);
  }
  brs.Finalize();

  std::vector< size_t > erased;
  erased.reserve(no_contains_->size());
  for (size_t i = 0; i < no_contains_->size(); ++i) {
    auto better_overlap = (*no_contains_)[i];
    auto overlap = better_overlap->overlap();
    // TODO(brahle): izdvoji ovo u zasebnu klasu iterator
    auto v1 = brs[overlap->read_one]->overlaps();
    auto v2 = brs[overlap->read_two]->overlaps();
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    bool done = false;
    while (!done && it1 != v1.end() && it2 != v2.end()) {
      if (it1->first == overlap->read_one || it1->first == overlap->read_two) {
        ++it1;
        continue;
      }
      if (it2->first == overlap->read_one || it2->first == overlap->read_two) {
        ++it2;
        continue;
      }
      if (it1->first == it2->first) {
        if (isTransitive(better_overlap, it1->second, it2->second)) {
          erased.emplace_back(i);
          done = true;
        }
        ++it1;
        ++it2;
      } else if (it1->first < it2->first) {
        ++it1;
      } else {
        ++it2;
      }
    }
  }

  no_transitives_ = BetterOverlapSetPtr(new BetterOverlapSet(reads_));
  size_t idx = 0;
  for (size_t i = 0; i < no_contains_->size(); ++i) {
    if (idx < erased.size() && i == erased[idx]) {
      ++idx;
      continue;
    }
    auto better_overlap = (*no_contains_)[i];
    no_transitives_->Add(new BetterOverlap(better_overlap));
  }
  int transitive_edge_count = no_contains_->size() - no_transitives_->size();
  fprintf(
      stderr,
      "Transitive edges = %d (%.2lf%%)\n",
      transitive_edge_count,
      (transitive_edge_count * 100.0) / no_contains_->size());
}

void Unitigging::makeContigs() {
  uint32_t** degrees = new uint32_t*[reads_->size()];
  for (size_t i = 0; i < reads_->size(); ++i) {
    degrees[i] = new uint32_t[2]();
  }
  auto add_degree =
      [] (uint32_t** degrees,
          uint32_t read,
          layout::BetterOverlap* better_overlap) {
    uint32_t suf = better_overlap->Suf(read);
    degrees[read][suf] += 1;
  };
  for (size_t i = 0; i < no_transitives_->size(); ++i) {
    auto better_overlap = (*no_transitives_)[i];
    auto overlap = better_overlap->overlap();
    add_degree(degrees, overlap->read_one, better_overlap);
    add_degree(degrees, overlap->read_two, better_overlap);
  }

  UnionFind uf(reads_->size());
  BetterReadSet brs(reads_, 1);
  contigs_ = ContigSetPtr(new ContigSet(&brs));

  for (size_t i = 0; i < no_transitives_->size(); ++i) {
    auto better_overlap = (*no_transitives_)[i];
    auto overlap = better_overlap->overlap();
    auto read_one = overlap->read_one;
    auto read_two = overlap->read_two;
    if (degrees[read_one][better_overlap->Suf(read_one)] == 1 &&
        degrees[read_two][better_overlap->Suf(read_two)] == 1) {
      auto contig_one = uf.find(read_one);
      auto contig_two = uf.find(read_two);
      // printf("Spajam %d i %d\n", contig_one, contig_two);
      auto larger = uf.join(read_one, read_two);
      if (larger == contig_one) {
        (*contigs_)[contig_one]->Join(better_overlap, (*contigs_)[read_two]);
      } else {
        (*contigs_)[contig_two]->Join(better_overlap, (*contigs_)[read_one]);
      }
    }
  }

  for (size_t i = 0; i < reads_->size(); ++i) {
    delete [] degrees[i];
  }
  delete [] degrees;
}

};  // namespace layout

