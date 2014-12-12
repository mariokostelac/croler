CC = g++
CFLAGS = -g -Wall -std=c++11 -O2 -I vendor/
CFLAGS_AMOS = -Wno-deprecated -O2
LDFLAGS = -lz -pthread

VPATH=obj:bin:src

.PHONY = amos

align = src/align/align.h src/align/align.cpp
minq = src/fixed_min_queue/fixed_min_queue.cpp
nucleo_buffer = src/nucleo_buffer/nucleo_buffer.cpp src/nucleo_buffer/nucleo_buffer.h
hash_list = src/hash_list/hash_list.cpp
minimizer = src/minimizer/minimizer.h src/minimizer/minimizer.cpp $(hash_list) $(nucleo_buffer)
overlap = src/overlap.cpp
parser = src/parser/parser.h
timer = src/timer/timer.h src/timer/timer.cpp
parsero = src/parsero/parsero.h
AMOS_OBJ_ALL = vendor/AMOS/obj/Overlap_AMOS.o \
							 vendor/AMOS/obj/Read_AMOS.o \
							 vendor/AMOS/obj/Sequence_AMOS.o \
							 vendor/AMOS/obj/Universal_AMOS.o \
							 vendor/AMOS/obj/Message_AMOS.o \
							 vendor/AMOS/obj/datatypes_AMOS.o \
							 vendor/AMOS/obj/BankStream_AMOS.o \
							 vendor/AMOS/obj/Bank_AMOS.o \
							 vendor/AMOS/obj/IDMap_AMOS.o \
							 vendor/AMOS/obj/utility_AMOS.o

default: amos prepare bin/overlap
all: prepare bin/test_nucleo_buffer bin/test_align bin/test_minq bin/test_hash_list bin/test_minimizer bin/overlap

prepare:
	@test -d bin || mkdir bin
	@test -d obj || mkdir obj

amos:
	@$(MAKE) -C vendor/AMOS

bin/overlap: overlap.o align.o nucleo_buffer.o minq.o minimizer.o timer.o
	@/bin/echo -e "\e[34m  LD $@ \033[0m"
	@$(CC) -o $@ $(AMOS_OBJ_ALL) $(^:%=obj/%) $(LDFLAGS)

bin/test_nucleo_buffer: $(nucleo_buffer) src/nucleo_buffer/test_nucleo_buffer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@

bin/test_minq: $(minq) src/fixed_min_queue/test_minq.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@

bin/test_hash_list: $(hash_list) src/hash_list/test_hash_list.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@

bin/test_minimizer: $(minimizer) src/minimizer/test_minimizer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

align.o: src/align/align.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $<

nucleo_buffer.o: src/nucleo_buffer/nucleo_buffer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $<

overlap.o: src/overlap.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $<

minq.o: src/fixed_min_queue/fixed_min_queue.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $<

minimizer.o: src/minimizer/minimizer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $^

timer.o: src/timer/timer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $<

bin/test_align: $(align) src/align/test_align.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean cleanall:
	@test -d bin && rm -r bin
	@test -d obj && rm -r obj
