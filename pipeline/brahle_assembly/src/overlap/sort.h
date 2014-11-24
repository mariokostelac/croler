#ifndef OVERLAP_SORT_H_
#define OVERLAP_SORT_H_

#include <cstdint>
#include <vector>


namespace overlap {


class ReadSet;
class String;

std::vector<uint32_t> STLStringOrder(const ReadSet& strings);
std::vector<uint32_t> RadixStringOrder(const ReadSet& strings, const size_t max_val);


}  // namespace overlap

#endif  // OVERLAP_SORT_H
