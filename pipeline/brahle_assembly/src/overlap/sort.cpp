#include <cstring>
#include <algorithm>
#include <functional>

#include "read.h"
#include "sort.h"


namespace overlap {


namespace {


class ReadCompare : public std::binary_function<uint32_t, uint32_t, bool> {
 public:
  ReadCompare(const ReadSet& reads) : reads_(reads) {}
  bool operator() (uint32_t one, uint32_t two) const {
    return (*reads_[one]) < (*reads_[two]);
  }

 private:
  const ReadSet& reads_;
};


class ReadLenCompare : public std::binary_function<uint32_t, uint32_t, bool> {
 public:
  ReadLenCompare(const ReadSet& reads) : reads_(reads) {}
  bool operator() (uint32_t one, uint32_t two) const {
    return reads_[one]->size() < reads_[two]->size();
  }

 private:
  const ReadSet& reads_;
};

}  // unnamed namespace


std::vector<uint32_t> STLStringOrder(const ReadSet& reads) {
  const uint32_t num_reads = reads.size();
  std::vector<uint32_t> order(num_reads);

  for (uint32_t idx = 0; idx < num_reads; ++idx) {
    order[idx] = idx;
  }
  std::sort(order.begin(), order.end(), ReadCompare(reads));

  return order;
}

std::vector<uint32_t> RadixStringOrder(const ReadSet& reads, const size_t max_val) {
  const uint32_t num_reads = reads.size();
  std::vector<uint32_t> order(num_reads);

  for (uint32_t idx = 0; idx < num_reads; ++idx) {
    order[idx] = idx;
  }
  std::sort(order.begin(), order.end(), ReadLenCompare(reads));

  const uint32_t max_size = reads[order[num_reads - 2]]->size();
  uint32_t start_idx = num_reads - 2;

  uint32_t* char_count = new uint32_t[max_val + 1];
  std::vector<uint32_t> tmpord(num_reads);

  for (int32_t pos = max_size - 1; pos >= 0; --pos) {
    // Move the index.
    while (start_idx > 0 && reads[order[start_idx - 1]]->size() > (uint32_t)pos) {
      --start_idx;
    }
    // Calculate character bucket sizes.
    memset(char_count, 0, sizeof(uint32_t) * (max_val + 1));
    for (uint32_t rix = start_idx; rix < num_reads; ++rix) {
      uint8_t val = (*reads[order[rix]])[pos];
      char_count[val]++;
    }
    // Calculate character bucket start positions.
    uint32_t sum = 0;
    for (uint32_t cix = 0; cix <= max_val; ++cix) {
      uint32_t tmp = char_count[cix];
      char_count[cix] = sum;
      sum += tmp;
    }
    // Iterate the values again and put them to their buckets.
    for (uint32_t rix = start_idx; rix < num_reads; ++rix) {
      uint8_t val = (*reads[order[rix]])[pos];
      tmpord[start_idx + char_count[val]++] = order[rix];
    }
    // Copy everything back to original order array.
    std::copy(tmpord.begin() + start_idx, tmpord.end(), order.begin() + start_idx);
  }

  delete[] char_count;
  return order;
}

}  // namespace overlap
