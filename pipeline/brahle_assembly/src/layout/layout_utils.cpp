// Copyright 2014 Bruno Rahle

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <utility>
#include <zlib.h>
#include <string>
#include <sstream>
#include <unordered_map>
using std::stringstream;
using std::string;
using std::unordered_map;

#include "layout/layout_utils.h"
#include "fastq_parser/parser.h"
#include "AMOS/src/foundation_AMOS.hh"

// init fasta/fastq reader
KSEQ_INIT(gzFile, gzread)

  namespace layout {


    /**
     * Since ids in ReadSet start with 0 and real ids (reads read from afg file/AMOS bank) can start with arbitrary number,
     * we have to map real_id -> internal_id (sequence that starts with 0).
     * That's why introduced this type.
     */
    typedef std::unordered_map<int, int> ReadIdMap;

    ReadIdMap _MapIds(overlap::ReadSet* reads) {
      ReadIdMap mapped;

      int reads_len = reads->size();
      for (int i = 0; i < reads_len; ++i) {
        auto& read = (*reads)[i];
        if (mapped.count(read->orig_id())) {
          fprintf(stderr, "Read with orig_id '%d' already seen\n", read->orig_id());
          exit(2);
        }
        mapped[read->orig_id()] = read->id();
      }

      return mapped;
    }

    /**
     * Reads data for one read from the .afg file.
     */
    uint8_t* _ReadDataForOneRead(FILE *fd) {
      static char buff[1 << 20];
      char *buff_free = buff;
      fscanf(fd, " seq: ");
      buff[0] = 0;
      do {
        // We remove one character for '\n'. However, in the first step, strlen is
        // 0, so we don't want to move in the negative direction.
        buff_free += std::max(0, static_cast<int>(strlen(buff_free) - 1));
        fgets(buff_free, sizeof(buff) - sizeof(char) * (buff_free - buff), fd);
      } while (buff_free[0] != '.');
      int length = buff_free - buff;
      char* data = new char[length+1];
      strncpy(data, buff, length);
      data[length] = 0;
      return reinterpret_cast<uint8_t*>(data);
    }

    overlap::Read* ReadOneReadAfg(FILE *fd, int pos) {
      static char buff[1 << 20];
      char *buff_free = buff;
      int id;
      int length;
      fscanf(fd, " iid:%d eid:%*s", &id);
      uint8_t *data = _ReadDataForOneRead(fd);
      do {
        fgets(buff, sizeof(buff), fd);
      } while (buff[0] != '.');
      fscanf(fd, " frg:%*d clr:%*d,%d", &length);
      return new overlap::Read(data, length, pos, id);
    }

    overlap::ReadSet* ReadReadsAfg(FILE *fd) {
      clock_t start = clock();
      auto read_set = new overlap::ReadSet(10000);
      int i = 0;
      char buff[1 << 20];
      while (fscanf(fd, " %s", buff) == 1) {
        if (!strcmp(buff, "{RED")) {
          read_set->Add(ReadOneReadAfg(fd, i));
          ++i;
        }
      }
      printf(
          "Reads read in %.2lfs\n",
          (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));
      return read_set;
    }

    overlap::ReadSet* ReadReadsSeq(char* filename) {
      gzFile fp = gzopen(filename, "r");
      kseq_t *seq = kseq_init(fp);

      clock_t start = clock();
      auto read_set = new overlap::ReadSet(10000);

      int len = 0;
      int id = 0;
      while ((len = kseq_read(seq)) > 0) {
        char *read_string = new char[len + 1];
        strcpy(read_string, seq->seq.s);
        auto read = new overlap::Read(reinterpret_cast<uint8_t*>(read_string), len, id, id);
        read_set->Add(read);
        id++;
      }

      kseq_destroy(seq);
      gzclose(fp);

      printf(
          "Reads read in %.2lfs\n",
          (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));
      return read_set;
    }

    /**
     * Reads all reads from AMOS bank.
     */
    overlap::ReadSet* ReadReadsAmos(const char *bank_name) {
      clock_t start = clock();
      auto read_set = new overlap::ReadSet(10000);

      AMOS::BankStream_t bank(AMOS::Read_t::NCODE);
      assert(bank.exists(bank_name));

      bank.open(bank_name, AMOS::B_READ);

      AMOS::Read_t read;
      bank.seekg(0, AMOS::BankStream_t::BEGIN);
      int index = 0;
      while (bank >> read) {
        AMOS::Range_t clear = read.getClearRange();
        string seq = read.getSeqString(clear);
        int len = seq.length();

        char *read_string = new char[len + 1];
        strcpy(read_string, seq.c_str());
        auto rd = new overlap::Read(reinterpret_cast<uint8_t*>(read_string), len, index, read.getIID());
        read_set->Add(rd);
        index++;
      }

      bank.close();

      printf(
          "Reads read in %.2lfs\n",
          (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));

      return read_set;
    }

    overlap::OverlapSet* ReadOverlapsAfg(overlap::ReadSet* read_set, FILE *fd) {
      clock_t start = clock();

      ReadIdMap internal_id = _MapIds(read_set);

      overlap::OverlapSet* overlap_set = new overlap::OverlapSet(10000);
      char type;
      int read_one;
      int read_two;
      int score;
      int hang_one;
      int hang_two;
      while (fscanf(
            fd,
            " {OVL adj:%c rds:%d,%d scr:%d ahg:%d bhg:%d }",
            &type,
            &read_one,
            &read_two,
            &score,
            &hang_one,
            &hang_two) == 6) {

        if (!internal_id.count(read_one)) {
          fprintf(stderr, "Read with orig_id '%d' has not been found\n", read_one);
          exit(3);
        }

        if (!internal_id.count(read_two)) {
          fprintf(stderr, "Read with orig_id '%d' has not been found\n", read_two);
          exit(3);
        }

        if (type != 'N' && type != 'I') {
          fprintf(stderr, "Unkown overlap type '%c'\n", type);
          exit(3);
        }

        read_one = internal_id[read_one];
        read_two = internal_id[read_two];

        std::pair<int, int> lenghts = getOverlapLengths(read_set, read_one, read_two, hang_one, hang_two);
        overlap_set->Add(new overlap::Overlap(
              read_one,
              read_two,
              lenghts.first,
              lenghts.second,
              hang_one,
              hang_two,
              type == 'N' ? overlap::Overlap::Type::EB : overlap::Overlap::Type::EE,
              0));
      }
      printf(
          "Overlaps read in %.2lfs\n",
          (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));
      return overlap_set;
    }

    /**
     * Reads all overlaps form the AMOS bank.
     */
    overlap::OverlapSet* ReadOverlapsAmos(overlap::ReadSet* read_set, const char* bank_name) {
      clock_t start = clock();
      overlap::OverlapSet* overlap_set = new overlap::OverlapSet(10000);

      AMOS::Overlap_t amos_overlap;
      AMOS::Message_t msg;
      AMOS::BankStream_t bank(AMOS::Overlap_t::NCODE);

      char type;
      int read_one;
      int read_two;
      int score;
      int hang_one;
      int hang_two;
      int len_one;
      int len_two;

      if(bank.exists(bank_name)) {
        bank.open(bank_name, AMOS::B_READ);
        if (bank.empty()) {
          return overlap_set;
        }

        while(bank >> amos_overlap) {

          auto reads = amos_overlap.getReads();
          read_one = reads.first - 1;
          read_two = reads.second - 1;
          score = amos_overlap.getScore();
          hang_one = amos_overlap.getAhang();
          hang_two = amos_overlap.getBhang();
          type = amos_overlap.getAdjacency();
          assert(type == 'N' || type == 'I');

          std::pair<int, int> lenghts = getOverlapLengths(read_set, read_one, read_two, hang_one, hang_two);

          overlap_set->Add(new overlap::Overlap(
                read_one,
                read_two,
                lenghts.first,
                lenghts.second,
                hang_one,
                hang_two,
                type == 'N' ? overlap::Overlap::Type::EB : overlap::Overlap::Type::EE,
                score));
        }

        bank.close();
      } else {
        fprintf(stderr, "AMOS Overlap bank %s does not exist.\n", bank_name);
        return overlap_set;
      }

      fprintf(stderr, "Overlaps read in %.2lfs\n",
          (clock() - start)/static_cast<double>(CLOCKS_PER_SEC));
      return overlap_set;
    }

    int n50(Unitigging::ContigSetPtr contig_set) {
      std::vector< int > v;
      int cnt = 0;
      v.reserve(contig_set->size()*10);
      for (size_t i = 0; i < contig_set->size(); ++i) {
        if ((*contig_set)[i]->IsUsable()) {
          printf("Contig %d: size = %d\n", ++cnt, (*contig_set)[i]->size());
          for (size_t j = 0; j < (*contig_set)[i]->size(); ++j) {
            v.push_back((*contig_set)[i]->size());
          }
        }
      }
      if (v.size() == 0) {
        return 0;
      }
      std::sort(v.begin(), v.end());
      return v[v.size()/2];
    }

    /**
     * Calculates lengths of overlaps between read_one and read_two.
     */
    std::pair<int, int> getOverlapLengths(const overlap::ReadSet* read_set, const int a, const int b, const int a_hang, const int b_hang) {
      int len_one, len_two;
      int a_size = (*read_set)[a]->size();
      int b_size = (*read_set)[b]->size();
      //
      // -------|--------------> b_hang
      // a_hang ---------------|------>
      //
      // -a_hang -------------|------->
      // --------|------------> -b_hang
      //
      // -------|-------------|------->
      // a_hang --------------> -b_hang
      //
      // -a_hang --------------> b_hang
      // --------|-------------|------>
      //
      if (a_hang >= 0 && b_hang >= 0) {
        len_one = a_size - a_hang;
        len_two = b_size - b_hang;
      } else if (a_hang <= 0 && b_hang <= 0) {
        len_one = a_size + b_hang;
        len_two = b_size + a_hang;
      } else if (a_hang >= 0 && b_hang <= 0) {
        len_one = a_size + b_hang - a_hang;
        len_two = b_size;
      } else if (a_hang <= 0 && b_hang >= 0) {
        len_one = a_size;
        len_two = b_size + a_hang - b_hang;
      } else {
        // case not covered
        assert(false);
      }

      assert(len_one > 0 && len_two > 0);
      return std::make_pair(len_one, len_two);
    }

    /**
     * Writes all usable contigs from contig set to a file with given filename.
     * If file does not exist, it will be created.
     * Returns the number of written contigs.
     */
    int ContigsToFile(std::shared_ptr<layout::ContigSet> contigs, const char *contigs_filename) {
      FILE *contigs_file = fopen(contigs_filename, "w");
      if (contigs_file == nullptr) {
        return -1;
      }

      int written = 0;
      int contigs_size = contigs->size();
      for (int i = 0; i < contigs_size; ++i) {
        // skip non-usable contigs
        if (!((*contigs)[i]->IsUsable())) continue;

        fprintf(contigs_file, "{LAY\n");
        uint32_t offset = 0;
        const std::deque< layout::BetterRead* > &reads = (*contigs)[i]->getReads();

        int num_reads = reads.size();
        for (int j = 0; j < num_reads - 1; ++j) {
          layout::BetterRead* read1 = reads[j];
          layout::BetterRead* read2 = reads[j + 1];
          read1->Finalize();
          const std::vector< std::pair< uint32_t, layout::BetterOverlap* >> &overlaps = read1->overlaps();

          // find overlap between first and second read
          for (const auto& overlap: overlaps) {
            if (overlap.first == read2->id() && overlap.second != nullptr) {
              fprintf(contigs_file, "{TLE\n");
              fprintf(contigs_file, "clr:%u,%u\n", 0, read1->read()->size());
              fprintf(contigs_file, "off:%u\n", offset);
              fprintf(contigs_file, "src:%d\n}\n", read1->read()->orig_id());
              offset += read1->read()->size() - overlap.second->Length();
              break;
            }
          }
        }

        // output last read
        layout::BetterRead* read = reads[num_reads - 1];
        fprintf(contigs_file, "{TLE\n");
        fprintf(contigs_file, "clr:%u,%u\n", 0, read->read()->size());
        fprintf(contigs_file, "off:%u\n", offset);
        fprintf(contigs_file, "src:%d\n}\n", read->read()->orig_id());

        fprintf(contigs_file, "}\n");
        written++;
      }

      fclose(contigs_file);
      return written;
    }

    std::string dot_graph(overlap::ReadSet* reads, overlap::OverlapSet* overlaps) {
      BetterReadSet brs(reads, false);
      BetterOverlapSet bos(reads, overlaps);
      return dot_graph(&brs, &bos);
    }

    std::string dot_graph(const BetterReadSet* reads, const BetterOverlapSet* overlaps) {
      stringstream graph;
      graph << "digraph overlaps {\n";
      int overlaps_size = overlaps->size();
      for (int i = 0; i < overlaps_size; ++i) {
        const auto& overlap = (*overlaps)[i];
        int read1 = (*reads)[overlap->overlap()->read_one]->read()->orig_id();
        int read2 = (*reads)[overlap->overlap()->read_two]->read()->orig_id();
        if (overlap->Suf(overlap->overlap()->read_one)) {
          graph << read1 << " -> " << read2;
        } else {
          graph << read2 << " -> " << read1;
        }
        if (overlap->overlap()->type == overlap::Overlap::Type::EB) {
          graph << " [color=green] ";
        } else {
          graph << " [color=red] ";
        }
        graph << ";\n";
      }
      graph << "}\n";
      return graph.str();
    }
  };  // namespace layout
