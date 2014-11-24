#ifndef OVERLAP_VALIDATE_H_
#define OVERLAP_VALIDATE_H_


namespace overlap {

class OverlapSet;
class ReadSet;

OverlapSet* ValidateCandidates(
    const ReadSet& reads,
    const OverlapSet& candidates);

}  // namespace overlap

#endif  // OVERLAP_VALIDATE_H_
