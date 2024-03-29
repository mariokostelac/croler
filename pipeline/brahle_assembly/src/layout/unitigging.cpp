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
    makeContigs(no_transitives_, reads_);
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
        (*reads_)[overlap->read_two]->addCoverage(1);
        fprintf(stderr, "Removing read %d as containment read (%d, %d)\n",
            (*reads_)[overlap->read_one]->orig_id(),
            (*reads_)[overlap->read_one]->orig_id(),
            (*reads_)[overlap->read_two]->orig_id()
            );
      } else if (better_overlap->two()->size() == overlap->len_two) {
        erased[overlap->read_two] = true;
        (*reads_)[overlap->read_one]->addCoverage(1);
        fprintf(stderr, "Removing read %d as containment read (%d, %d)\n",
            (*reads_)[overlap->read_two]->orig_id(),
            (*reads_)[overlap->read_one]->orig_id(),
            (*reads_)[overlap->read_two]->orig_id()
            );
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

    // for o1(x, y) and o2(x, a), o3(y, a)
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
            (*reads_)[overlap->read_one]->addCoverage(
                static_cast<double> (better_overlap->Length()) /
                (*reads_)[overlap->read_one]->size());
            (*reads_)[it1->first]->addCoverage(
                static_cast<double> (better_overlap->Length()) /
                (*reads_)[it1->first]->size());
            erased.emplace_back(i);

            const auto& o = better_overlap->overlap();
            const auto& o1 = it1->second->overlap();
            const auto& o2 = it2->second->overlap();
            fprintf(stderr,
                "Removing overlap (%d, %d) as transitive of (%d, %d), (%d, %d)\n",
                (*reads_)[o->read_one]->orig_id(),
                (*reads_)[o->read_two]->orig_id(),
                (*reads_)[o1->read_one]->orig_id(),
                (*reads_)[o1->read_two]->orig_id(),
                (*reads_)[o2->read_one]->orig_id(),
                (*reads_)[o2->read_two]->orig_id()
                );
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

  void Unitigging::makeContigs(BetterOverlapSetPtr& c_overlaps, overlap::ReadSet*& read_set) {
    fprintf(stderr, "Reads: %d, Overlaps: %d\n", read_set->size(), c_overlaps->size());
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
    for (size_t i = 0; i < c_overlaps->size(); ++i) {
      auto better_overlap = (*c_overlaps)[i];
      auto overlap = better_overlap->overlap();
      add_degree(degrees, overlap->read_one, better_overlap);
      add_degree(degrees, overlap->read_two, better_overlap);
    }

    UnionFind uf(reads_->size());
    // first mark all reads as unusable
    for (size_t i = 0; i < reads_->size(); ++i) {
      (*reads_)[i]->usable(false);
    }
    // mark reads as usable
    for (size_t i = 0; i < read_set->size(); ++i) {
      size_t id = (*read_set)[i]->id();
      (*reads_)[id]->usable(true);
    }
    BetterReadSet brs(reads_, 1);
    // @mculinovic adding overlaps
    for (size_t i = 0; i < c_overlaps->size(); ++i) {
      auto better_overlap = (*c_overlaps)[i];
      auto overlap = better_overlap->overlap();
      brs[overlap->read_one]->AddOverlap(better_overlap);
      brs[overlap->read_two]->AddOverlap(better_overlap);
    }
    brs.Finalize();
    // end adding overlaps

    // create contigs from reads
    contigs_ = ContigSetPtr(new ContigSet(&brs));

    for (size_t i = 0; i < c_overlaps->size(); ++i) {
      auto better_overlap = (*c_overlaps)[i];
      auto overlap = better_overlap->overlap();
      auto read_one = overlap->read_one;
      auto read_two = overlap->read_two;
      if (degrees[read_one][better_overlap->Suf(read_one)] == 1 &&
          degrees[read_two][better_overlap->Suf(read_two)] == 1) {
        auto contig_one = uf.find(read_one);
        auto contig_two = uf.find(read_two);
        if (contig_one == contig_two) {
          fprintf(stderr, "Detected circular contig %d; reads:", contig_one);
          auto& reads = (*contigs_)[contig_one]->getReads();
          for (auto it = reads.begin(); it != reads.end(); ++it) {
            const auto& read = (*it)->read();
            fprintf(stderr, "%d(size: %d) ", read->orig_id(), read->hi() - read->lo());
          }
          fprintf(stderr, "\n");
          continue;
        }

        fprintf(stderr, "Joining contigs %d and %d by reads (%d, %d)\n", contig_one, contig_two, read_one, read_two);
        auto larger = uf.join(read_one, read_two);
        if (larger == contig_one) {
          (*contigs_)[contig_one]->Join(better_overlap, (*contigs_)[contig_two]);
        } else {
          (*contigs_)[contig_two]->Join(better_overlap, (*contigs_)[contig_one]);
        }
      } else {
        fprintf(stderr,
            "(%d, %d) is not a candidate for contig because of degrees d[%d][%d] = %d; d[%d][%d] = %d\n",
            read_one,
            read_two,
            read_one,
            better_overlap->Suf(read_one),
            degrees[read_one][better_overlap->Suf(read_one)],
            read_two,
            better_overlap->Suf(read_two),
            degrees[read_two][better_overlap->Suf(read_two)]
            );
      }
    }
    for (size_t i = 0; i < reads_->size(); ++i) {

      delete [] degrees[i];
    }
    delete [] degrees;
  }

};  // namespace layout

