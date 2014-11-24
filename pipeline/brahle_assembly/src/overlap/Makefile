#!/usr/bin/env bash
CXX=g++
CFLAGS=-Wall -Wextra -Werror -pedantic -O3 -std=c++11
LDFLAGS=
# Sources.
SRCDIR=.
SOURCES=main.cpp read.cpp fm_index.cpp overlap.cpp sort.cpp suffix_array.cpp suffix_filter.cpp
# Objects.
OBJDIR=../../obj
OBJECTS=$(SOURCES:.cpp=.o) sais.o
# Binaries.
BINDIR=../../bin
BINARIES=olapper

all : $(addprefix $(SRCDIR)/, $(SOURCES)) $(BINARIES)

$(BINARIES) : $(addprefix $(OBJDIR)/, $(OBJECTS))
	$(CXX) $^ -o $(BINDIR)/$@ $(LDFLAGS)

$(OBJDIR)/sais.o : $(SRCDIR)/sais.c
	$(CXX) $(CFLAGS) -fomit-frame-pointer -DNDEBUG  -c -o $@ $^

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) $< -c -o $@

clean :
	rm -f $(BINDIR)/* $(OBJDIR)/*
