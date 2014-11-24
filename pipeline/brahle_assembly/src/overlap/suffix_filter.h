#ifndef OVERLAP_SUFFIX_FILTER_H_
#define OVERLAP_SUFFIX_FILTER_H_

#include <stdint.h>
#include <sys/types.h>
#include <queue>
#include <unordered_map>

#include <overlap/util.hpp>


namespace overlap {


class FMIndex;
class OverlapSet;
class Read;
class ReadSet;


class SuffixFilter {
 public:
  SuffixFilter(
      double error_rate,
      size_t min_overlap_size);

  virtual ~SuffixFilter();

  virtual OverlapSet* FindCandidates(
      const ReadSet& reads,
      const UintArray& read_order,
      const FMIndex& fmi) = 0;

  virtual OverlapSet* FilterCandidates(
      const OverlapSet& candidates) = 0;

  static size_t FactorSize(double error_rate, size_t min_overlap_size);

 protected:
  const size_t min_overlap_size_;
  const size_t factor_size_;
};


class BFSSuffixFilter : public SuffixFilter {
 public:
  BFSSuffixFilter(
      double error_rate,
      size_t min_overlap_size);

  ~BFSSuffixFilter();

  OverlapSet* FindCandidates(
      const ReadSet& reads,
      const UintArray& read_order,
      const FMIndex& fmi);

  OverlapSet* FilterCandidates(
      const OverlapSet& candidates);

 private:
  class BFSContext {
   public:
    BFSContext(
        const Read& read,
        const UintArray& read_order,
        const FMIndex& fmi,
        size_t factor_size,
        size_t min_overlap_size,
        OverlapSet* results);

    void Start(
        uint32_t start_pos,
        uint32_t max_error);

    void Clear();

   private:
//    typedef std::tuple<uint32_t, uint32_t, uint32_t> State;

    struct State {
      uint32_t low;
      uint32_t high;
      uint32_t pos;
    };

    struct StateHash {
      std::size_t operator()(const State& k) const;
    };

    struct StateEqual {
      bool operator()(const State& lhs, const State& rhs) const;
    };

    typedef std::queue<State> StateQueue;
    typedef std::unordered_map<State, uint32_t, StateHash, StateEqual> StateMap;

    void CheckOverlaps(
        uint32_t low,
        uint32_t high,
        uint32_t start,
        uint32_t pos,
        uint32_t error);

    void Queue(
        uint32_t low,
        uint32_t high,
        uint32_t pos,
        uint32_t error,
        StateQueue& queue,
        bool can_inc);

    const Read& read_;
    const UintArray& read_order_;
    const FMIndex& fmi_;
    const size_t factor_size_;
    const size_t min_overlap_size_;

    OverlapSet* results_;

    StateQueue queue_[2];
    StateMap states_;
  };
};


}  // namespace overlap

#endif  // OVERLAP_SUFFIX_FILTER_H_
