#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <cstdlib>
#include <vector>
#include <utility>
#include <unordered_map>
#include "../nucleo_buffer/nucleo_buffer.h"
#include "../fixed_min_queue/fixed_min_queue.cpp"
#include "../hash_list/hash_list.cpp"

struct minimizer_t {
    unsigned int pos;
    nstring_t str;
    nstring_t ordstr;
    minimizer_t(nstring_t str, int pos) : pos(pos), str(str) {
        // invert some bits
        ordstr = str ^ 0x33333333;
    }
};

class Minimizer {
 public:
     explicit Minimizer(int minimizer_len, int window_len);
     ~Minimizer();
     void calculate_and_store(int str_index, const char *str);
     void calculate_and_get(std::vector<minimizer_t>& container, const char *str);
     const HashList<nstring_t, std::pair<unsigned int, unsigned int> >& get_minimizers() const;
     void shrink();
 private:
     int minimizer_len;
     int window_len;
     NucleoBuffer* buff;
     FixedMinQueue<minimizer_t>* q;
     HashList<nstring_t, std::pair<unsigned int, unsigned int> > minimizers;
};
#endif
