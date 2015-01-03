
#include "reader.h"

#include <cassert>
#include <cstring>

#define BUFF_SIZE 4096
#define MAX_SEQ 32768

namespace AMOS {

  enum ReaderState {
    OUT,
    IN_READ,
    IN_SEQ,
    IN_QLT,
  };

  int get_reads(std::vector<Read*>& container, const char* afg_filename) {
    int records = 0;

    FILE *fd = fopen(afg_filename, "r");
    if (fd == nullptr) {
      return -1;
    }

    char line[BUFF_SIZE] = {0};

    ReaderState state = OUT;
    Read* curr_read = new Read();
    char curr_seq[MAX_SEQ] = {0};
    int curr_seq_len = 0;
    int line_len = 0;
    long last_pos = ftell(fd);
    const char **dst_str = nullptr;
    while (!feof(fd)) {
      // read next line and update pointers
      fgets(line, BUFF_SIZE, fd);
      long curr_pos = ftell(fd);
      assert(curr_pos - last_pos < BUFF_SIZE);
      line_len = curr_pos - last_pos;
      last_pos = curr_pos;

      // strip \n
      line[line_len - 1] = 0;
      line_len--;

      switch (state) {
        case OUT:
          if (strstr(line, "{RED") != nullptr) {
            state = IN_READ;
          }
          break;
        case IN_READ:
          if (sscanf(line, " iid: %d ", &curr_read->iid)) {
            state = IN_READ;
          } else if (strstr(line, "seq:") != nullptr) {
            dst_str = &(curr_read->seq);
            state = IN_SEQ;
          } else if (strstr(line, "qlt:") != nullptr) {
            dst_str = &(curr_read->qlt);
            state = IN_QLT;
          } else if (sscanf(line, " clr: %d, %d ", &curr_read->clr_lo, &curr_read->clr_hi)) {
            state = IN_READ;
          } else if (line[0] == '}') {
            state = OUT;
            container.push_back(curr_read);
            curr_read = new Read();
            records++;
          }
          break;
        case IN_SEQ:
        case IN_QLT:
          if (line[0] == '.') {
            assert(dst_str != nullptr);
            char *cpy = new char[curr_seq_len + 1];
            strncpy(cpy, curr_seq, curr_seq_len + 1);
            curr_seq_len = 0;
            *dst_str = cpy; // assign copied string to the destination
            dst_str = nullptr;
            state = IN_READ;
          } else {
            assert(line_len < MAX_SEQ - curr_seq_len);
            strcpy(curr_seq + curr_seq_len, line);
            curr_seq_len += line_len;
          }
          break;
        default:
          assert(false);
      }
    }
    delete curr_read;

    fclose(fd);
    return records;
  }
}
