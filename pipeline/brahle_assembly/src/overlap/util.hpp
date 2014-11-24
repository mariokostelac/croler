#ifndef OVERLAP_UTIL_H_
#define OVERLAP_UTIL_H_

#include <cstdint>
#include <vector>

namespace overlap {

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&) = delete; \
  TypeName& operator=(const TypeName&) = delete

typedef std::vector<uint32_t> UintArray;
typedef std::vector<int32_t> IntArray;

}  // namespace overlap

#endif  // OVERLAP_UTIL_H_
