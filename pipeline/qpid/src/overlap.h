#include <vector>
#include "./minimizer/minimizer.h"
#include "../lib/ThreadPool/ThreadPool.h"
using std::vector;

struct offset_t {
    unsigned int index;
    int lo_offset;
    int hi_offset;
    offset_t() {}
    offset_t(unsigned int index, int offset) : index(index), lo_offset(offset), hi_offset(offset) {}
};

typedef HashList<nstring_t, std::pair<unsigned int, unsigned int>> minimizers_t;
