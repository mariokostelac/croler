#include <cstdint>
#include <algorithm>
#include <map>
#include <memory>

#include "myers.h"
#include "overlap.h"
#include "read.h"
#include "validate.h"


namespace overlap {


namespace {

uint32_t OverlapConfidence(Overlap* o) {
  if (o) return 1;
  return 0;
}

}  // unnamed namespace

OverlapSet* ValidateCandidates(
    const ReadSet& reads,
    const OverlapSet& candidates) {

  typedef std::pair<uint32_t, uint32_t> Pair;
  std::map<Pair, Overlap*> best;

  std::unique_ptr<OverlapSet> overlaps(new OverlapSet(candidates.size()));

  const double error_rate = 0.01;
  const double extra_ratio = 1.3;

  for (uint32_t oid = 0; oid < candidates.size(); ++oid) {
    Overlap* o = candidates[oid];

    Read* read_one = reads[o->read_one];
    Read* read_two = reads[o->read_two];

    int len_two, score;
    int ret = MyersEditDistance(
        read_one->data() + read_one->size() - o->len_one,
        o->len_one,
        read_two->data(),
        std::min((int)(o->len_one * extra_ratio), (int)read_two->size()),
        5,
        (int)(o->len_one * error_rate * 2),
        MYERS_MODE_SHW,
        &score,
        &len_two);

    if (ret == MYERS_STATUS_OK and len_two != -1 and score != -1) {
      o->len_two = len_two + 1;
      o->score = score;

      auto elem = best.find(std::make_pair(o->read_one, o->read_two));
      if (elem == best.end() or elem->second->score < o->score) {
        best[std::make_pair(o->read_one, o->read_two)] = o;
        best[std::make_pair(o->read_two, o->read_one)] = o;
      }
    }
  }

  for (auto it : best) {
    const Pair& p = it.first;
    Overlap* o = it.second;

    if (o->read_one == p.first and o->read_two == p.second) {
      overlaps->Add(new Overlap(*o));
    }
  }

  return overlaps.release();
}

}  // namespace overlap
