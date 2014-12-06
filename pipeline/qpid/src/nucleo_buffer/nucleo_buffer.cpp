#include "./nucleo_buffer.h"
#include <cassert>
#include <cstdio>

NucleoBuffer::NucleoBuffer(int size): size(size) {
    assert(size > 0 && size <= 16);
    buffer = 0;
    content_mask = 0;
    for (int i = 0; i < size; ++i) {
        content_mask = content_mask << 2;
        content_mask |= 0x3;
    }
}

inline char NucleoBuffer::get_repr(char ch) {
    ch = toupper(ch);
    if (ch == 'A') return 0x0;
    if (ch == 'C') return 0x1;
    if (ch == 'G') return 0x2;
    if (ch == 'T') return 0x3;
    return 0x0;
}

nstring_t NucleoBuffer::write(char ch) {
    buffer = buffer << 2;
    buffer |= get_repr(ch);
    return buffer;
}

nstring_t NucleoBuffer::get_content() {
    return buffer & content_mask;
}

int NucleoBuffer::get_size() {
    return size;
}
