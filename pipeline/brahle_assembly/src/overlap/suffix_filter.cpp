#include <cassert>
#include <cmath>
#include <algorithm>
#include <memory>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "fm_index.h"
#include "overlap.h"
#include "read.h"
#include "suffix_filter.h"
#include "util.hpp"


namespace overlap {


namespace {

void FilterPass(
    const std::vector<Overlap*>& invec,
    std::vector<Overlap*>& outvec) {

  Overlap* prev = nullptr;
  for (Overlap* curr : invec) {
    if (prev == nullptr) {
      outvec.push_back(curr);
      prev = curr;
    } else if (prev->read_one == curr->read_one and
               prev->read_two == curr->read_two and
               prev->type == curr->type) {
      if ((int32_t)abs(curr->len_one - prev->len_one) > (prev->score - curr->score) + 1) {
        outvec.push_back(curr);
        prev = curr;
      }
    } else {
      outvec.push_back(curr);
      prev = curr;
    }
  }
}

}  // unnamed namespace

SuffixFilter::SuffixFilter(double error_rate, size_t min_overlap_size)
    : min_overlap_size_(min_overlap_size),
      factor_size_(SuffixFilter::FactorSize(error_rate, min_overlap_size)) {
}

SuffixFilter::~SuffixFilter() {
}

size_t SuffixFilter::FactorSize(double error_rate, size_t min_overlap_size) {
  size_t final_size = 1000;
  for(size_t size = min_overlap_size; size < 1000; ++size) {
    size_t new_size = (size_t)(ceil(size / ceil(error_rate * size) + 1) + 1e-9);
    final_size = std::min(final_size, new_size);
  }
  return final_size;
}

BFSSuffixFilter::BFSSuffixFilter(double error_rate, size_t min_overlap_size)
    : SuffixFilter(error_rate, min_overlap_size) {
}

BFSSuffixFilter::~BFSSuffixFilter() {
}

OverlapSet* BFSSuffixFilter::FindCandidates(const ReadSet& reads,
    const UintArray& read_order, const FMIndex& fmi) {

  const size_t num_reads = reads.size();
  std::unique_ptr<OverlapSet> candidates(
      new OverlapSet(num_reads * num_reads / 2));

  for (uint32_t read_idx = 0; read_idx < num_reads; ++read_idx) {
    const Read& read = *reads[read_idx];
    const size_t read_size = read.size();

    for (int32_t pos = read_size - 1;
         pos + 1 >= (int32_t)factor_size_ * 2;
         pos -= factor_size_) {

      std::unique_ptr<BFSContext> bfs(
        new BFSContext(read, read_order, fmi, factor_size_,
                       min_overlap_size_, candidates.get()));

      bfs->Start(pos, (size_t)pos == read_size - 1);
      //bfs->Clear();
    }
  }

  return candidates.release();
}

OverlapSet* BFSSuffixFilter::FilterCandidates(const OverlapSet& candidates) {
  std::vector<Overlap*> cont;

  for (uint32_t idx = 0; idx < candidates.size(); ++idx) {
    Overlap* curr = candidates[idx];
    if (curr->read_one != curr->read_two) {
      cont.push_back(curr);
    }
  }

  std::vector<Overlap*> fwd, bwd;
  FilterPass(cont, fwd);

  std::reverse(fwd.begin(), fwd.end());
  FilterPass(fwd, bwd);

  std::reverse(bwd.begin(), bwd.end());
  std::unique_ptr<OverlapSet> filtered(new OverlapSet(bwd.size()));

  for (Overlap* o : bwd) {
    filtered->Add(new Overlap(*o));
  }

  return filtered.release();
}

BFSSuffixFilter::BFSContext::BFSContext(const Read& read, const UintArray& read_order,
    const FMIndex& fmi, size_t factor_size, size_t min_overlap_size,
    OverlapSet* results)
    : read_(read),
      read_order_(read_order),
      fmi_(fmi),
      factor_size_(factor_size),
      min_overlap_size_(min_overlap_size),
      results_(results),
      states_(100000) {
}

void BFSSuffixFilter::BFSContext::Start(uint32_t start_pos, uint32_t error) {
  Queue(0, fmi_.size(), 0, error, queue_[0], false);

  for (size_t qid = 0; !queue_[qid].empty(); qid = 1 - qid) {
    std::queue<State>& curr = queue_[qid];
    std::queue<State>& next = queue_[1 - qid];

    uint32_t low, high, pos;
    while(!curr.empty()) {
      low = curr.front().low;
      high = curr.front().high;
      pos = curr.front().pos;
      error = states_[curr.front()];

      CheckOverlaps(low, high, start_pos, pos, error);

      if (error > 0 && pos <= start_pos) {
        Queue(low, high, pos + 1, error - 1, next, true);
      }

      uint32_t newlow, newhigh;
      for (uint8_t cix = 1; cix <= fmi_.max_val(); ++cix) {
        newlow = fmi_.Less(cix) + fmi_.Rank(cix, low);
        newhigh = fmi_.Less(cix) + fmi_.Rank(cix, high);
        if (newlow >= newhigh) continue;

        if (pos <= start_pos) {
          if (cix == read_[start_pos - pos]) {
            Queue(newlow, newhigh, pos + 1, error, curr, true);
          } else if (error > 0) {
            Queue(newlow, newhigh, pos + 1, error - 1, next, true);
          }
        } else if (error > 0) {
          Queue(newlow, newhigh, pos, error - 1, next, false);
        }
      }

      curr.pop();
    }
  }
}

void BFSSuffixFilter::BFSContext::Clear() {
  assert(queue_[0].empty() and queue_[1].empty());
  states_.clear();
}

void BFSSuffixFilter::BFSContext::CheckOverlaps(uint32_t low, uint32_t high,
    uint32_t start, uint32_t pos, uint32_t error) {

  size_t overlap_size = pos + read_.size() - start - 1;
  if (pos >= factor_size_ && overlap_size >= min_overlap_size_) {
    low = fmi_.Rank(0, low);
    high = fmi_.Rank(0, high);

    for (uint32_t idx = low; idx < high; ++idx) {
      results_->Add(
          new Overlap(read_.id(), read_order_[idx], overlap_size,
                      overlap_size, Overlap::Type::EB, error));
    }
  }
}

void BFSSuffixFilter::BFSContext::Queue(uint32_t low, uint32_t high,
    uint32_t pos, uint32_t error, StateQueue& queue, bool can_inc) {
  State new_state = {low, high, pos};
  if (states_.find(new_state) == states_.end()) {
    queue.push(new_state);
    states_[new_state] = error + (can_inc && !(pos % factor_size_) ? 1 : 0);
  }
}

std::size_t BFSSuffixFilter::BFSContext::StateHash::operator()(
    const BFSSuffixFilter::BFSContext::State& k) const {
  /*
  return (std::hash<uint32_t>()(std::get<0>(k)) ^
          std::hash<uint32_t>()(std::get<1>(k)) ^
          std::hash<uint32_t>()(std::get<2>(k)));
  */
  register size_t state;
  state = ((k.low << 5) + k.low) ^ k.high;
  return ((state << 5) + state) ^ k.pos;
}

bool BFSSuffixFilter::BFSContext::StateEqual::operator()(
    const BFSSuffixFilter::BFSContext::State& lhs,
    const BFSSuffixFilter::BFSContext::State& rhs) const {
  /*
  return (std::get<0>(lhs) == std::get<0>(rhs) &&
          std::get<1>(lhs) == std::get<1>(rhs) &&
          std::get<2>(lhs) == std::get<2>(rhs));
          */
  return (lhs.low == rhs.low and
          lhs.high == rhs.high and
          lhs.pos == lhs.pos);
}

}  // namespace overlap
