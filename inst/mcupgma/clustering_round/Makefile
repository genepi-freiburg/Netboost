##########################################################################
# MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA #
#             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.#
# 									 #
# Makefile - compiles multi-round clustering code compilation		 #
#									 #
# Copyright (C) Elon Portugaly, Yaniv Loewenstein                        #              
#                School of Computer Science And Engineering              # 
#                Hebrew University of Jerusalem                          #
#                                                                        #
#      All Rights Reserved                                               #
#                                                                        #
#      This source code is distributed under the terms of the            #
#      GNU General Public License. See the file LICENSE                  #
#      for details.                                                      #
#									 #
##########################################################################

# When using this makefile for a new project, all one has to do is:
#	Define EXEC to be the list of possible target executables
#	Define SOURCES to be the list of the rest of the .cpp files
#	Write the explicit dependencies of each target executable. These include only the .cpp files the executable depends on.
# Do not use makedepend with this makefile. This makefile automaticly generates a .d file for each .o files. The .d file lists all dependencies of the .o file. This makefile then includes all the .d files.
# This makefile does not support dynamic libraries and library creation



#######################################################################################
#  SETTINGS         ###################################################################
#######################################################################################
# Define compiler and flags
CXX = g++
CC = gcc
COMPILER = ${CXX}

# compilation flags (e.g. -I paths - search paths for <> include clauses)
CPPFLAGS_DEBUG = -g -ftemplate-depth-64 -I./
CPPFLAGS = -O3 -ftemplate-depth-64 -I./
CFLAGS_DEBUG = -g -I./
CFLAGS = -O3 -I./
LDFLAGS_DEBUG = -g
LDFLAGS = -O3


#ARCH := $(shell arch)_$(shell uname -s)
ARCH := $(shell uname -m)_$(shell uname -s)

OBJDIR := obj.$(ARCH)
DEPDIR := dep.$(ARCH)
BINDIR := bin.$(ARCH)

# executable names
PROGRAMS = hierarchical_clustering test_heap

# make program file (EXEC) list from program name list (PROGRAMS)
EXEC = $(patsubst %, $(BINDIR)/%, $(PROGRAMS))
EXEC_DEBUG = $(patsubst %, $(BINDIR).debug/%, $(PROGRAMS))


# source lists...
SOURCES = $(wildcard *.cpp) $(wildcard */*.cpp)
GGOFILES = $(wildcard *.ggo)
GGOHFILES = $(patsubst %.ggo, %.cmdline.h, $(GGOFILES))
GGOCFILES = $(patsubst %.ggo, %.cmdline.c, $(GGOFILES))
OFILESWITHGGO = $(patsubst %.ggo, $(OBJDIR)/%.o, $(GGOFILES))
OFILESWITHGGO_DEBUG = $(patsubst %.ggo, $(OBJDIR).debug/%.o, $(GGOFILES))

# make .o file lists from source lists - using pattern substitution
OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o,$(SOURCES)) $(patsubst %.ggo, $(OBJDIR)/%.cmdline.o,$(GGOFILES))
OBJECTS_DEBUG = $(patsubst %.cpp, $(OBJDIR).debug/%.o,$(SOURCES)) $(patsubst %.ggo, $(OBJDIR).debug/%.cmdline.o,$(GGOFILES))

# make .d file list from source lists
# %.d files are makefiles including all the dependecies of %.o
DEPENDENCYFILES = $(patsubst %.cpp,$(DEPDIR)/%.d, $(SOURCES))


######################################################################

# The default target
all: install $(EXEC)
#	$(MAKE)
install: $(OBJDIR) $(BINDIR) $(DEPDIR)

tags: *.cpp *.hpp */*.cpp */*.hpp
	etags --language=c++ $^ 

# cleans up the place
clean:
	-rm -f $(OBJECTS) $(OBJECTS_DEBUG) $(GGOCFILES) $(DEPENDENCYFILES) core

cleanprograms: clean
	-rm -f $(EXEC) $(EXEC_DEBUG)
arch:
	@echo $(ARCH)
PHONY: all clean tags $(PROGRAMS) arch

#######################################################################

# make the necessary directories for the current system (OS and architecture) (lonshy)
$(OBJDIR): 
	mkdir -p $(OBJDIR) 
$(BINDIR):
	mkdir -p $(BINDIR)
$(DEPDIR):
	mkdir -p $(DEPDIR) 

# explicit dependencies of executables
$(BINDIR)/hierarchical_clustering:       $(OBJDIR)/HierarchicalClustering_main.o	$(OBJDIR)/HierarchicalClustering_main.cmdline.o		$(OBJDIR)/HierarchicalClustering_with_unknown_edges.o
$(BINDIR).debug/hierarchical_clustering: $(OBJDIR).debug/HierarchicalClustering_main.o	$(OBJDIR).debug/HierarchicalClustering_main.cmdine.o	$(OBJDIR).debug/HierarchicalClustering_with_unknown_edges.o
$(BINDIR)/test_heap: $(OBJDIR)/test_heap.o $(OBJDIR)/test_heap.cmdline.o
$(BINDIR).debug/test_heap: $(OBJDIR).debug/test_heap.o $(OBJDIR).debug/test_heap.cmdline.o




# include %.d files
-include $(DEPENDENCYFILES)



#######################################################################################
#  TARGETS          ###################################################################
#######################################################################################





# pattern rules for making different kind of files

$(PROGRAMS): %: $(BINDIR)/%

# rules for making %.o files
# the rest of the dependencies of %.o are listed in %.d
$(OBJDIR)/%.o: %.cpp
	${COMPILER} ${CPPFLAGS} -c $< -o $@

$(OFILESWITHGGO): $(OBJDIR)/%.o : %.cmdline.h

$(OFILESWITHGGO_DEBUG): $(OBJDIR).debug/%.o : %.cmdline.h

$(OBJDIR).debug/%.o: %.cpp
	${COMPILER} ${CPPFLAGS_DEBUG} -c $< -o $@

$(OBJDIR)/%.cmdline.o: %.cmdline.c
	${CC} ${CFLAGS} -c $< -o $@

$(OBJDIR).debug/%.cmdline.o: %.cmdline.c
	${CC} ${CFLAGS_DEBUG} -c $< -o $@

# rules for making %.d files
# the %.d file will includes the rest of the dependencies of itself
$(DEPDIR)/%.d: %.cpp $(DEPDIR)
	test -d $(DEPDIR) && $(COMPILER) -MM -MT '$@ $(patsubst $(DEPDIR)/%.d, $(OBJDIR)/%.o, $@) $(patsubst $(DEPDIR)/%.d, $(OBJDIR).debug/%.o, $@)' $(CPPFLAGS) $< > $@


# rules for making %.cmdline.h files
# these files are created using getgetopt from %.ggo files
# the command also generates the %.cmdline.c files
%.cmdline.c %.cmdline.h: %.ggo
	gengetopt -i $< -F $*.cmdline

# rules for making executables from object files and static libraries
$(EXEC): %:
	${COMPILER} $(filter %.o,$^) $(patsubst lib%.a,-l%,$(filter lib%.a,$^)) -o $@  ${LDFLAGS}
$(EXEC_DEBUG): %:
	${COMPILER} $(filter %.o,$^) $(patsubst lib%.a,-l%,$(filter lib%.a,$^)) -o $@  ${LDFLAGS_DEBUG}
