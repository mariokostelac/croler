#include <ctype.h>
#include <cstdint>

typedef unsigned int nstring_t;

class NucleoBuffer {
 public:
        explicit NucleoBuffer(int size);
        nstring_t write(char ch);
        nstring_t get_content();
        int get_size();
 private:
        int size;
        nstring_t buffer;
        nstring_t content_mask;
        char get_repr(char ch);
};


