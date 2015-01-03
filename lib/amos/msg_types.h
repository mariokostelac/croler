
#ifndef _AMOS_MSG_TYPES_H
#define _AMOS_MSG_TYPES_H

#include <cstdint>

namespace AMOS {

  typedef struct Read {
    uint32_t iid;
    uint32_t clr_lo;
    uint32_t clr_hi;
    const char *seq;
    const char *qlt;

    Read() : seq(nullptr), qlt(nullptr) {}

    Read(uint32_t iid, uint32_t clr_lo, uint32_t clr_hi, const char *seq, const char *qlt)
      :iid(iid), clr_lo(clr_lo), clr_hi(clr_hi), seq(seq), qlt(qlt)
    {}

    ~Read() {
      if (seq != nullptr) {
        delete[] seq;
      }
      if (qlt != nullptr) {
        delete[] qlt;
      }
    }
  } Read;
}
#endif
