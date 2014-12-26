#ifndef READ_H
#define READ_H

typedef struct Read {
  int id;
  char* sequence;

  Read() {}
  Read(int idx, char* seq): id(idx), sequence(seq) {}
} Read;
#endif
