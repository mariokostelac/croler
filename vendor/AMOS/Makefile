default: prepare build

ALL_OBJ = Overlap_AMOS.o \
					Read_AMOS.o \
					Sequence_AMOS.o \
				 	Universal_AMOS.o \
					Message_AMOS.o \
					datatypes_AMOS.o \
					BankStream_AMOS.o \
					Bank_AMOS.o \
					IDMap_AMOS.o \
					utility_AMOS.o

CC = g++
CFLAGS = -Wno-deprecated -O2
LDFLAGS =

prepare:
	@test -d obj || mkdir obj

build: $(ALL_OBJ)

clean:
	@test -d obj && rm -r obj/*

$(ALL_OBJ): %.o: src/%.cc
	@/bin/echo -e "\e[34m  CC $@ \033[0m"
	@$(CC) $(CFLAGS) -c -o obj/$@ $<
