CC = g++
CFLAGS = -g -Wall -std=c++11 -O2 -I ./vendor -I ./src
CFLAGS_AMOS = -Wno-deprecated -O2
LDFLAGS = -lz -pthread

VPATH=obj:bin:src:vendor

.PHONY = amos

align = align/align.h align/align.cpp
minq = fixed_min_queue/fixed_min_queue.cpp
nucleo_buffer = nucleo_buffer/nucleo_buffer.cpp nucleo_buffer/nucleo_buffer.h
hash_list = hash_list/hash_list.cpp
minimizer = minimizer/minimizer.h minimizer/minimizer.cpp $(hash_list) $(nucleo_buffer)
overlap = overlap.cpp
parser = parser/parser.h
timer = timer/timer.h timer/timer.cpp
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

bin/overlap: $(addprefix obj/,overlap.o align.o nucleo_buffer.o minq.o minimizer.o timer.o)
	@/bin/echo -e "\e[34m  LD $@ \033[0m"
	@$(CC) -o $@ $(AMOS_OBJ_ALL) $^ $(LDFLAGS)

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

obj/align.o: src/align/align.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o $@ $<

obj/nucleo_buffer.o: src/nucleo_buffer/nucleo_buffer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o $@ $<

obj/overlap.o: src/overlap.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o $@ $<

obj/minq.o: src/fixed_min_queue/fixed_min_queue.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o $@ $<

obj/minimizer.o: src/minimizer/minimizer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o $@ $^

obj/timer.o: src/timer/timer.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o $@ $<

bin/test_align: $(align) src/align/test_align.cpp
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@rm -f $@
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean cleanall:
	@test -d bin && rm -r bin
	@test -d obj && rm -r obj
