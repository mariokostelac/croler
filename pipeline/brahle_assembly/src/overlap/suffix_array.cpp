#include <stddef.h>
#include <algorithm>

#include "read.h"
#include "sais.h"
#include "suffix_array.h"


namespace overlap {


SACA::SACA() {}

SACA::~SACA() {}

SaisSACA::SaisSACA() {}

String* SaisSACA::BuildBWT(
    const ReadSet& reads,
    size_t depth) {
  size_t num_strings = reads.size();
  size_t bwt_size = num_strings;

  for (uint32_t read_idx = 0; read_idx < num_strings; ++read_idx) {
    bwt_size += reads[read_idx]->size();
  }

  int* super_string = new int[bwt_size + 1];
  uint32_t ss_pos = 0;
  for (uint32_t read_idx = 0; read_idx < num_strings; ++read_idx) {
    const Read* read = reads[read_idx];
    super_string[ss_pos++] = read_idx;
    for (uint32_t chr_idx = 0; chr_idx < read->size(); ++chr_idx) {
      super_string[ss_pos++] = (*read)[chr_idx] + num_strings - 1;
    }
  }
  super_string[bwt_size] = '\0';

  int* bwt32 = new int[bwt_size + 1];
  bwt32[bwt_size] = '\0';

  int32_t* A = new int32_t[bwt_size + 1];
  int32_t ret = sais_int_bwt(super_string, bwt32, A, (int)bwt_size, num_strings + depth);
  if (ret < 0) {
    fprintf(stderr, "SAIS failed (%d).\n", ret);
    return nullptr;
  }

  uint8_t* bwt = new uint8_t[bwt_size + 1];
  for (uint32_t bwt_pos = 0; bwt_pos <= bwt_size; ++bwt_pos) {
    int c = bwt32[bwt_pos];
    bwt[bwt_pos] = ((size_t)c < num_strings ? 0 : c - num_strings + 1);
  }

  delete[] super_string;
  delete[] bwt32;
  delete[] A;

  return new String(bwt, bwt_size);
}


}  // namespace overlap
