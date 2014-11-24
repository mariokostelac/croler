#ifndef OVERLAP_FM_INDEX_H_
#define OVERLAP_FM_INDEX_H_

#include <stdint.h>
#include <sys/types.h>
#include <memory>

namespace overlap {


class String;

class FMIndex {
 public:
  FMIndex(const String& bwt, size_t max_val);
  virtual ~FMIndex();

  // Accessors.
  size_t size() const;
  size_t max_val() const;
  // Return count of values less than 'chr' in complete BWT.
  virtual uint32_t Less(uint8_t chr) const = 0;
  // Return count of values equal to 'chr' in first 'pos' values of BWT.
  virtual uint32_t Rank(uint8_t chr, uint32_t pos) const = 0;

 protected:
  const size_t size_;
  const size_t max_val_;
};

class BucketedFMIndex : public FMIndex {
 public:
  BucketedFMIndex(const String& bwt, size_t max_val, size_t bucket_size);
  ~BucketedFMIndex();

  uint32_t Less(uint8_t chr) const;
  uint32_t Rank(uint8_t chr, uint32_t pos) const;

 private:
  void Init();

  const uint8_t* bwt_data_;
  // Cumulative sum of total char counts.
  uint32_t* const char_counts_;
  // Bucket count data.
  const size_t bucket_size_;
  const size_t num_buckets_;
  uint32_t* const buckets_;
};


}  // namespace overlap

#endif  // OVERLAP_FM_INDEX_H_
