#include <cassert>
#include <cstdio>
#include <cstring>

#include "read.h"


namespace overlap {


String::String(const uint8_t* data, size_t size)
    : data_(data), size_(size) {
}

String::~String() {
  delete[] data_;
}

const uint8_t* String::data() const {
  return data_;
}

size_t String::size() const {
  return size_;
}

const uint8_t& String::operator[] (const uint32_t idx) const {
  return data_[idx];
}

bool String::operator< (const String& other) const {
  return strcmp((const char*)data_, (const char*)other.data_) < 0;
}

Read::Read(const uint8_t* data, size_t size, uint32_t id, uint32_t orig_id)
    : String(data, size),
      id_(id),
      orig_id_(orig_id) {
}

Read::~Read() {
}

uint32_t Read::id() const {
  return id_;
}

uint32_t Read::orig_id() const {
  return orig_id_;
}

const uint8_t* ReverseComplement(const uint8_t* data, size_t size) {
  uint8_t* data_rc = new uint8_t[size + 1];
  data_rc[size] = '\0';
  for (uint32_t idx = 0; idx < size; ++idx) {
    data_rc[size - idx - 1] = 4 - data[idx];
  }
  return data_rc;
}

Read* ReverseComplement(const Read& read) {
  return new Read(
      ReverseComplement(read.data(), read.size()),
      read.size(),
      read.id(),
      read.orig_id());
}

void PrintRead(FILE* fd, const Read& read) {
  for (uint32_t idx = 0; idx < read.size(); ++idx) {
    fprintf(fd, "%c", (char)read[idx] + 'A');
  }
  fprintf(fd, "\n");
}

ReadSet::ReadSet(size_t capacity) {
  reads_.reserve(capacity);
}

ReadSet::~ReadSet() {
  for (Read* r : reads_) {
    delete r;
  }
}

void ReadSet::Add(Read* read) {
  reads_.push_back(read);
}

const Read* ReadSet::Get(uint32_t read_idx) const {
  return reads_[read_idx];
}

size_t ReadSet::size() const {
  return reads_.size();
}

Read* const& ReadSet::operator[] (const uint32_t idx) const {
  return reads_[idx];
}

Read*& ReadSet::operator[] (const uint32_t idx) {
  return reads_[idx];
}

ReadSet* ReadFasta(FILE* fd, size_t min_read_size) {
  char buff[100000];
  ReadSet* read_set = new ReadSet(1 << 16);

  uint32_t id = 0;
  for (uint32_t orig_id = 0; fgets(buff, sizeof buff, fd); ++orig_id) {
    assert(fgets(buff, sizeof buff, fd));
    size_t size = strlen(buff);

    if (size >= min_read_size) {
      uint8_t* read_data = new uint8_t[size];
      buff[size - 1] = '\0';
      strcpy((char*)read_data, buff);

      for (uint32_t pos = 0; pos < size - 1; ++pos) {
        uint8_t fc = 0;
        switch (read_data[pos]) {
        case 'A':
          fc = 1;
          break;
        case 'C':
          fc = 2;
          break;
        case 'G':
          fc = 3;
          break;
        case 'T':
          fc = 4;
          break;
        default:
          assert(0);
          break;
        }
        read_data[pos] = fc;
      }

      read_set->Add(new Read(read_data, size - 1, id++, orig_id));
      assert(read_set->Get(id - 1)->id() == id - 1);
    }
  }

  return read_set;
}


}  // namespace overlap
