#ifndef OVERLAP_READ_H_
#define OVERLAP_READ_H_

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <vector>

#include "util.hpp"


namespace overlap {


class String {
public:
  String(const uint8_t* data, size_t size);
  virtual ~String();

  const uint8_t* data() const;
  size_t size() const;

  const uint8_t& operator[](const uint32_t idx) const;
  bool operator< (const String& other) const;

protected:
  const uint8_t* data_;
  const size_t size_;
};

class Read : public String {
public:
  Read(const uint8_t* data, size_t size, uint32_t id, uint32_t orig_id);
  ~Read();

  uint32_t id() const;
  uint32_t orig_id() const;

private:
  const uint32_t id_;
  const uint32_t orig_id_;
};

const uint8_t* ReverseComplement(const uint8_t* data, size_t size);
Read* ReverseComplement(const Read& read);

void PrintRead(FILE* fd, const Read& read);

class ReadSet {
public:
  ReadSet(size_t capacity);
  ~ReadSet();

  void Add(Read* read);
  const Read* Get(uint32_t read_idx) const;

  size_t size() const;

  Read* const& operator[](const uint32_t idx) const;
  Read*& operator[](const uint32_t idx);

private:
  std::vector<Read*> reads_;

  DISALLOW_COPY_AND_ASSIGN(ReadSet);
};

ReadSet* ReadFasta(FILE* fd, size_t min_read_size);

}  // namespace overlap

#endif  // OVERLAP_READ_H_
