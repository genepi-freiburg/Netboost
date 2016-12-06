#########################################################################
#									 #
# MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA #
#             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.#
#									 #
# Copyright (C) Yaniv Loewenstein				 	 #
#               School of Computer Science And Engineering		 #
#               Hebrew University of Jerusalem				 #
#									 #
#     All Rights Reserved						 #
#									 #
#     This source code is distributed under the terms of the		 #
#     GNU General Public License. See the file LICENSE			 #
#     for details.							 #
#									 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#									 #
# Makefile for compiling clustering utilities for 			 #
#   multi-round clustering (with parallel external merging), and for 	 #
#   single-round MC-UPGMA clustering 		 		         #
#									 #
# Note that only part of the sources are included with the MC-UPGMA      # 
#  distribution, hence some of the unrealted targets can't be made       #
#									 #
# Yaniv Loewenstein, 2007						 #
#									 #
# Revisions: 								 #
#									 #
# Aug 2008 - YL added documentation and MC-UPGMA specific targets	 #
#									 #
#########################################################################



###################################################################
###   configurable compilation options 				  #
###################################################################

# enable flags that help compile problems with STL tempaltes
COMPILE = 1
# implant debugging code for STL 
#RUNTIME = 1
# enable performance flags, and supress some runtime sanity checks by assert macros (-DNDEBUG)
PERFORM = 1
# enable profiling using gprof
#GPROF    = 1

ifdef COMPILE
#CXX = /cs/phd/ofirpele/utils/stl_decrypt/gfilt -Wall
CXXFLAGS +=  -Wall
endif

ifdef RUNTIME
#LDFLAGS  += -lefence  
#added 1/8/07:
CXXFLAGS += -D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG
endif

ifdef PERFORM
CXXFLAGS += -O3
CXXFLAGS += -DNDEBUG
endif

ifdef GPROF
CXXFLAGS += -pg
endif

####################################################################
# we compile C program (gengetopt cmdlines) with C++ compiler      #
####################################################################
CC = $(CXX)
CXX_EXEC = $(CXX) $(CXXFLAGS) $(LDFLAGS)

######################################################################
## seperate binaries for different architectures into different dirs # 
######################################################################
# architecture spec - OS and hardware, 
ARCH := $(shell uname -m)_$(shell uname -s)
BINDIR := bin.$(ARCH)
OBJDIR := obj.$(ARCH)

####################################################################
## Actual targets 						   #
####################################################################
OFILES := ClusteringUtil.o Tree.o
#OFILES = EdgeReader.o ClusteringUtil.o

targets := edges2bin bin2edges \
	test_heap_static test_heap1 test_heap2 \
	test_use_copy1	test_use_copy2 \
	edges2bin_mat bin_mat2edges \
	DiskReader DiskWriter \
	hcluster edges2valid_clusters\
	test_IteratorAdaptor\
	edge_collator
#	DiskEdgeManager.o
#	nodes2parents average_edges 

# add path to OFILES and targets var
targets := $(patsubst %,$(BINDIR)/%,$(targets))
OFILES  := $(patsubst %,$(OBJDIR)/%,$(OFILES)) 


############################################################################################
.PHONY: verbose

#all: $(targets) HierarchicalClustering.o
all: verbose $(OBJDIR) $(BINDIR) $(targets) $(OBJDIR)/HierarchicalClustering.o

verbose:
#	@echo OFILES = $(OFILES)
edit:
	xemacs *pp *ggo Makefile &   

clean:	
	rm -f $(targets) $(OFILES) *~ *.o
uninstall:
	rm -rf $(BINDIR) $(OBJDIR)

###########################################################
### 	MC-UPGMA targets and FILE LIST	###################
###########################################################


###################################################
## files needed for MC-UPGMA realted executables ## 
###################################################
export MCUPGMA_SOURCES := ClusteringUtil.hpp ClusteringUtil.cpp \
Tree.hpp Tree.cpp \
HashingMergerMain.cpp HashingMergerMain.ggo HashingMergerMain.cmdline.h HashingMergerMain.cmdline.c \
EdgeCollatorMain.cpp  EdgeCollatorMain.ggo  EdgeCollatorMain.cmdline.h  EdgeCollatorMain.cmdline.c \
Definitions.hpp Edge.hpp IteratorAdaptor.hpp \
Makefile LICENSE

###################################################
## MC-UPGMA related targets                      ##
###################################################
mcupgma: mcupgma_install
mcupgma_install: install_mcupgma
mcupgma_uninstall: uninstall_mcupgma
install_mcupgma: $(BINDIR) $(OBJDIR) $(BINDIR)/edges2valid_clusters $(BINDIR)/edge_collator
uninstall_mcupgma: unistall
mcupgma_tar: mcupgma_util.tar.gz 
mcupgma_util.tar.gz: $(MCUPGMA_SOURCES)
	tar -czf $@ $^


##########################################################################################
$(OBJDIR):
	mkdir -p $@

$(BINDIR):
	mkdir -p $@


$(OBJDIR)/%.o: %.cpp %.hpp
	$(CXX_EXEC) $(CXXFLAGS) -c -o $@ $<

#################### Object  files  ######################################################

#Heap.o:       Heap.cpp Heap.hpp Edge.hpp
#EdgeReader.o: EdgeReader.cpp EdgeReader.hpp Edge.hpp
#DiskEdgeManager.o: DiskEdgeManager.cpp DiskEdgeManager.hpp Edge.hpp ClusteringUtil.hpp

$(OBJDIR)/Tree.o: Tree.cpp Tree.hpp Definitions.hpp Edge.hpp ClusteringUtil.hpp 
$(OBJDIR)/ClusteringUtil.o: ClusteringUtil.cpp ClusteringUtil.hpp Definitions.hpp Edge.hpp
$(OBJDIR)/HierarchicalClustering.o: HierarchicalClustering.cpp HierarchicalClustering.hpp Edge.hpp Definitions.hpp ClusteringUtil.hpp Heap.hpp

################### Executables ##########################################################

$(BINDIR)/test_IteratorAdaptor: test_IteratorAdaptor.cpp IteratorAdaptor.hpp Edge.hpp
	$(CXX_EXEC)  -o $@ $<  
$(BINDIR)/DiskWriter:DiskEdgeManager.cpp DiskEdgeManager.hpp Edge.hpp ClusteringUtil.hpp 
	$(CXX_EXEC)  -o $@ $<   -DWRITE
$(BINDIR)/DiskReader:DiskEdgeManager.cpp DiskEdgeManager.hpp Edge.hpp ClusteringUtil.hpp
	$(CXX_EXEC) -o $@ $<   -DREAD
$(BINDIR)/test_use_copy1: test_use_copy.cpp 
	$(CXX_EXEC) -o $@ $<   -DBY_PTR
$(BINDIR)/test_use_copy2: test_use_copy.cpp 
	$(CXX_EXEC) -o $@ $<   -UBY_PTR

$(BINDIR)/test_heap_static: test_heap_static.cpp Edge.hpp Heap.hpp
	$(CXX_EXEC) -o $@ $< 
$(BINDIR)/edges2bin:	edges2bin.cpp Edge.hpp 
	$(CXX_EXEC) -o $@ $< -DTOBIN
$(BINDIR)/bin2edges:	edges2bin.cpp Edge.hpp
	$(CXX_EXEC) -o $@ $< -UTOBIN
$(BINDIR)/test_edge_reader1: test_edge_reader.cpp Edge.hpp  
	$(CXX_EXEC) -o $@ $< -UBY_ITER  
$(BINDIR)/test_edge_reader2: test_edge_reader.cpp Edge.hpp  
	$(CXX_EXEC) -o $@ $< -DBY_ITER  
$(BINDIR)/test_heap1: test_heap.cpp  Heap.hpp Definitions.hpp Edge.hpp
	$(CXX_EXEC) -o $@ $<  -UTRY_INPUT_RANGE
$(BINDIR)/test_heap2: test_heap.cpp  Heap.hpp Definitions.hpp Edge.hpp
	$(CXX_EXEC) -o $@ $<  -DTRY_INPUT_RANGE

$(BINDIR)/edges2bin_mat: edges2bin_mat.cpp Edge.hpp
	$(CXX_EXEC) -o $@ $< -DWRITE
$(BINDIR)/bin_mat2edges: edges2bin_mat.cpp Edge.hpp
	$(CXX_EXEC) -o $@ $< -DREAD
$(BINDIR)/hcluster: hcluster.cpp  Heap.hpp 
	$(CXX_EXEC) -o $@ $<

$(BINDIR)/edge_collator: EdgeCollatorMain.cpp $(OBJDIR)/EdgeCollatorMain.cmdline.o $(OFILES)
	$(CXX_EXEC) -o $@ $^ 

$(BINDIR)/edges2valid_clusters:  HashingMergerMain.cpp  $(OBJDIR)/HashingMergerMain.cmdline.o $(OFILES)
	$(CXX_EXEC)  -o $@ $< $(OBJDIR)/HashingMergerMain.cmdline.o $(OFILES)

#nodes2parents: nodes2parents.cpp  Definitions.hpp Edge.hpp ClusteringUtil.o
#	$(CXX_EXEC) -o $@ $< ClusteringUtil.o
#average_edges: average_edges.cpp  Definitions.hpp Edge.hpp
#	$(CXX_EXEC) -o $@ $< 



########################################################
## command line options parsers by gengetopt 	      ##
########################################################

# for cmdline gengetopt c files
$(OBJDIR)/%.cmdline.o: %.cmdline.c %.cmdline.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.cmdline.h %.cmdline.c:%.ggo
	gengetopt < $< --file-name=$*.cmdline

######## EOF ############################################################################
