default: amos prepare run

INCLUDE_SRC=./src
INCLUDE_VENDOR=./vendor

VPATH=src:src/layout:src/overlap:src/test:bin

# these headers are different. first, they don't exists
# but *.hpp files exists. it means that they contain template
# implementation and everything must be recompiled when
# they are changed.
HPP=util.h

OBJ=fm_index.o overlap.o better_overlap.o better_read.o read.o sort.o suffix_array.o suffix_filter.o unitigging.o util.o union_find.o contig.o layout_utils.o vertex.o string_graph.o label.o
OBJ_SPECIAL=sais.o
ALL_OBJ=$(OBJ) $(OBJ_SPECIAL)

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

EXE=main_layout

include buildnumber.mak

.PHONY: run valgrind amos

OPTIMIZATION_FLAGS=-flto -finline-limit=200
ASORTED_FLAGS=-std=c++11
WARINIG_FLAGS=-Wall -Wextra -pedantic -Werror
FLAGS=$(OPTIMIZATION_FLAGS) $(ASORTED_FLAGS) $(WARNING_FLAGS) $(BUILD_NUMBER_LDFLAGS)
DEBUG_FLAGS=-g -ggdb -DDEBUG
NODEBUG_FLAGS=-s -O3 -DNDEBUG

CC=g++ $(FLAGS) $(DEBUG_FLAGS) -I$(INCLUDE_SRC) -I$(INCLUDE_VENDOR)
CCE=g++ $(FLAGS) $(DEBUG_FLAGS) -I$(INCLUDE_SRC) -I$(INCLUDE_VENDOR)

TEST=unitigging_test

prepare:
	@test -d bin || mkdir bin

amos:
	@$(MAKE) -C vendor/AMOS

clean:
	@/bin/echo -e "\e[34m  Building AMOS \033[0m"
	rm bin/*

$(HPP): %.h: %.hpp
	@touch bin/$@

$(OBJ): %.o: %.cpp %.h $(patsubst %,%pp,$(HPP))
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) -c $< -o bin/$@

sais.o:	src/overlap/sais.c
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CCE) -fomit-frame-pointer -DNDEBUG  -c -o bin/$@ $^

# build rule for test executables.
$(TEST): %: src/test/%.cpp $(ALL_OBJ)
	@/bin/echo -e "\e[34m  LINK $@ \033[0m"
	@$(CC) $< $(patsubst %,bin/%,$(ALL_OBJ)) -lz -o bin/$@

# build rule for main executables.
$(EXE): %: src/%.cpp $(ALL_OBJ) $(BUILD_NUMBER_FILE)
	@/bin/echo -e "\e[34m  LINK $@ \033[0m"
	@$(CC) $< $(patsubst %,bin/%,$(ALL_OBJ)) $(AMOS_OBJ_ALL) -lz -o bin/$@

valgrind: main_layout
	valgrind bin/main_layout sample/small

run: $(EXE) $(BUILD_NUMBER_FILE)
	@bin/main_layout

lint:
	cpplint --root=src --filter=-legal/copyright,-readability/braces src/*.cpp src/layout/*.cpp src/layout/*.h

test: $(TEST)
	@for t in $(TEST) ; do \
		/bin/echo -e "\e[34m$$t\033[0m" ; \
		bin/$$t ; \
	done

valgrind-test: $(TEST)
	@for t in $(TEST) ; do \
		valgrind bin/$$t ; \
	done
