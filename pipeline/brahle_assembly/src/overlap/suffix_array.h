#ifndef OVERLAP_BWT_H_
#define OVERLAP_BWT_H_

#include <stdint.h>
#include <sys/types.h>


namespace overlap {


class ReadSet;
class String;

class SACA {
  public:
    SACA();
    virtual ~SACA();

    virtual String* BuildBWT(const ReadSet& reads, size_t depth) = 0;
};


class SaisSACA : public SACA {
  public:
    SaisSACA();

    String* BuildBWT(const ReadSet& reads, size_t depth);
};


}  // namespace overlap

#endif  // OVERLAP_BWT_H_
