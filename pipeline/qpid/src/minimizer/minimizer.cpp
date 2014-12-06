#include "./minimizer.h"
#include <cstring>
#include <cstdio>
#include <cassert>
#include <vector>
#include <utility>
#include <algorithm>
#include "../hash_list/hash_list.cpp"

bool operator==(const minimizer_t& lhs, const minimizer_t& rhs) {
    return lhs.pos == rhs.pos && lhs.str == rhs.str;
}

bool operator<(const minimizer_t& lhs, const minimizer_t& rhs) {
    if (lhs.str == rhs.str)
        return lhs.pos < rhs.pos;
    return lhs.ordstr < rhs.ordstr;
}

bool operator>(const minimizer_t& lhs, const minimizer_t& rhs) {
    if (lhs.str == rhs.str)
        return lhs.pos > rhs.pos;
    return lhs.ordstr > rhs.ordstr;
}

bool operator!=(const minimizer_t& lhs, const minimizer_t& rhs) {
    return !operator==(lhs, rhs);
}

Minimizer::Minimizer(int mlen, int wlen) {
    assert(mlen <= wlen);
    minimizer_len = mlen;
    window_len = wlen;
    buff = new NucleoBuffer(mlen);
    q = new FixedMinQueue<minimizer_t>(wlen - mlen + 1);
}

Minimizer::~Minimizer() {
    if (buff != 0)  delete(buff);
    if (q != 0)     delete(q);
}

void  Minimizer::calculate_and_store(int str_index, const char *str) {
    int len = strlen(str);
    assert(len >= window_len);

    q->clear();

    for (int i = 0; i < minimizer_len; ++i) {
        buff->write(str[i]);
    }
    q->push(minimizer_t(buff->get_content(), 0));

    for (int i = minimizer_len; i < window_len; ++i) {
        buff->write(str[i]);
        q->push(minimizer_t(buff->get_content(), i - minimizer_len + 1));
    }

    minimizer_t prev = q->min();
    this->minimizers.add(prev.str, std::make_pair(str_index, prev.pos));

    for (int i = window_len; i < len; ++i) {
        buff->write(str[i]);
        q->push(minimizer_t(buff->get_content(), i - minimizer_len + 1));
        if (q->min() != prev) {
            prev = q->min();
            this->minimizers.add(prev.str, std::make_pair(str_index, prev.pos));
        }
    }
}

void Minimizer::calculate_and_get(std::vector<minimizer_t>& container, const char *str) {
    int len = strlen(str);
    assert(len >= window_len);

    // we need it here to have a thread safe implementation
    FixedMinQueue<minimizer_t> q(window_len - minimizer_len + 1);
    NucleoBuffer buff(minimizer_len);

    for (int i = 0; i < minimizer_len; ++i) {
        buff.write(str[i]);
    }
    q.push(minimizer_t(buff.get_content(), 0));

    for (int i = minimizer_len; i < window_len; ++i) {
        buff.write(str[i]);
        q.push(minimizer_t(buff.get_content(), i - minimizer_len + 1));
    }

    minimizer_t prev = q.min();
    container.push_back(prev);

    for (int i = window_len; i < len; ++i) {
        buff.write(str[i]);
        q.push(minimizer_t(buff.get_content(), i - minimizer_len + 1));
        if (q.min() != prev) {
            prev = q.min();
            container.push_back(prev);
        }
    }
}

void Minimizer::shrink() {
    minimizers.shrink();
}

const HashList<nstring_t, std::pair<unsigned int, unsigned int> >& Minimizer::get_minimizers() const {
    return minimizers;
}
