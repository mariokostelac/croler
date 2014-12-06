CC = g++
CFLAGS = -g -Wall -std=c++11 -O2 -I vendor/
LDFLAGS = -lz -pthread

align = src/align/align.h src/align/align.cpp
minq = src/fixed_min_queue/fixed_min_queue.cpp
nucleo_buffer = src/nucleo_buffer/nucleo_buffer.cpp src/nucleo_buffer/nucleo_buffer.h
hash_list = src/hash_list/hash_list.cpp
minimizer = src/minimizer/minimizer.h src/minimizer/minimizer.cpp $(hash_list) $(nucleo_buffer)
overlap = src/overlap.cpp
parser = src/parser/parser.h
timer = src/timer/timer.h src/timer/timer.cpp
parsero = src/parsero/parsero.h

all: bin/test_nucleo_buffer bin/test_align bin/test_minq bin/test_hash_list bin/test_minimizer bin/overlap

bin/test_nucleo_buffer: $(nucleo_buffer) src/nucleo_buffer/test_nucleo_buffer.cpp
	@echo [CC] $@
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@

bin/test_minq: $(minq) src/fixed_min_queue/test_minq.cpp
	@echo [CC] $@
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@

bin/test_hash_list: $(hash_list) src/hash_list/test_hash_list.cpp
	@echo [CC] $@
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@

bin/test_minimizer: $(minimizer) src/minimizer/test_minimizer.cpp
	@echo [CC] $@
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

obj/align.o: src/align/align.cpp
	@echo [CC] $@
	@$(CC) $(CFLAGS) -c -o $@ $^

obj/nucleo_buffer.o: src/nucleo_buffer/nucleo_buffer.cpp
	@echo [CC] $@
	@$(CC) $(CFLAGS) -c -o $@ $^

obj/overlap.o: src/overlap.cpp
	@echo [CC] $@
	@$(CC) $(CFLAGS) -c -o $@ $^

obj/minq.o: src/fixed_min_queue/fixed_min_queue.cpp
	@echo [CC] $@
	@$(CC) $(CFLAGS) -c -o $@ $^

obj/minimizer.o: src/minimizer/minimizer.cpp
	@echo [CC] $@
	@$(CC) $(CFLAGS) -c -o $@ $^

obj/timer.o: src/timer/timer.cpp
	@echo [CC] $@
	@$(CC) $(CFLAGS) -c -o $@ $^

bin/overlap: obj/overlap.o obj/align.o obj/nucleo_buffer.o obj/minq.o obj/minimizer.o obj/timer.o
	@echo [LD] $@
	@$(CC) $(LDFLAGS) -o $@ obj/*.o \
		vendor/AMOS/obj/Read_AMOS.o \
	  vendor/AMOS/obj/Sequence_AMOS.o \
	  vendor/AMOS/obj/Universal_AMOS.o \
	  vendor/AMOS/obj/Message_AMOS.o \
	  vendor/AMOS/obj/datatypes_AMOS.o \
	  vendor/AMOS/obj/BankStream_AMOS.o \
	  vendor/AMOS/obj/Bank_AMOS.o \
	  vendor/AMOS/obj/IDMap_AMOS.o \
	  vendor/AMOS/obj/utility_AMOS.o

bin/test_align: $(align) src/align/test_align.cpp
	@echo [CC] $@
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean cleanall:
	rm -rf bin/*
