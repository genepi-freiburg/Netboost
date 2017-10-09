##########################################################################
# MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA #
#             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.#
#                                                                        #
# Copyright (C) 2007, Yaniv Loewenstein                                  #
#                School of Computer Science And Engineering              #
#                Hebrew University of Jerusalem                          #
#                                                                        #
#      All Rights Reserved                                               #
#                                                                        #
#      This source code is distributed under the terms of the            #
#      GNU General Public License. See the file LICENSE                  #
#      for details.                                                      #
#                                                                        #
##########################################################################

include definitions.mk

ifndef INPUT_EDGES_FILES
$(error user must define variable $$INPUT_EDGES_FILES externally)
endif


# max number of iterations , current iteration,next iteration 
T      = 10
t      = 0
# increment (operator++ like):
next_t = $(shell perl -e '{$$a = shift; print ++$$a}' $(t))

# \psi in the paper (missing values are assumed to be that)
export PSI = 1000
export OUTDIR = iteration_$t
TREE := tree
export CUMULATIVE_TREE_FILE = $(OUTDIR)/tree
export TREE_FILE = $(OUTDIR)/tree_part$t
export INPUT_EDGES_FILES
#export OSIZES_FILE = $(OUTDIR)/osizes

##############################################################
#######       TARGETS                                    #####
##############################################################

.PHONY:test_input $(INPUT_EDGES_FILES)

all: verbose test_input $(OUTDIR) $(TREE) $(OUTDIR)/outputs_listing

# Test that directories are not mistakenly given as input files
test_input: $(INPUT_EDGES_FILES)

$(INPUT_EDGES_FILES):
# white space important here:
	@test ! -d $@ || (echo "*** ERROR *** : seems like input file: $@ - is a directory, not a valid input file" >> /dev/stderr ; exit 101)
	@test -e $@   || (echo "*** ERROR *** : seems like input file: $@ - does not exist" >> /dev/stderr ; exit 102)
	@echo "==================================       t = $(t)    ==============================="
	@echo t = $(t)

############################################################################
### the complete tree is built by concatanating tree parts that are newer than the current tree
### How to rebuild the tree from existing tree parts, from a previous run? solutions:
###
### ------------------------- DEPRECATED ---------------------------------
### Now if we remove the overall tree, and rebuild it by concatanating its ready made parts, after the
### the first part, it (tree)  will be newer than the next parts (tree_partX) and thus they will
### falsely not be included in the overall tree.
### To make sure that this does not happen we touch the next part if it exists, so its timestamp becomes newer 
### than the partial tree which are currently constructing, thus forcing the make mechanism to invoke this target again 
### with the next part again and again, when it exists. If it does not exist, it will be remade and there is no
### need to artificially modify its timestamp.
### ----------------------------------------------------------------------
###
### ------------------------- CURRENT SOLUTION ---------------------------
###  We assume that successive tree parts have increasing time stamps. Now if the cumulative 
###  overall tree (i.e. the concatanation of tree parts) is damaged, we want to rebuild it by 
### concatanation of *all* parts in the right order. If the parts already exist, once we create 
### the tree from the first part, it will have a newer timestamp than the next parts, and will 
### thus be missing these parts. To fix this, we modify the timestamp of the overall tree, to 
### be as new as the last concatanated part.
### ----------------------------------------------------------------------
###
######################################################################


$(TREE): $(TREE_FILE)
	@test -s $< || (rm -vf $< ; echo "*** ERROR *** : seems like tree file: $@  is empty (no clustering done in this round). This might happen erroniously if your input includes edges >= $(PSI) (\psi - the user supplied parameter for maximum dissimilarity)." >> /dev/stderr ; exit 104)
	cat $< >> $@ && touch -r $< $@
###	test -s iteration_$(next_t)/tree_part$(next_t) && touch iteration_$(next_t)/tree_part$(next_t)


#############################
#############################
$(TREE_FILE): 
	@echo "---------------------------- start   running clustering_iteration.mk ---------------"
	$(MAKE) -f $(SCRIPTS_DIR)/clustering_iteration.mk 
	sleep $(SLEEP_SECONDS_AFTER_CLUSTERING)
	@echo "---------------------------- finish  running clustering_iteration.mk ---------------"


#############################
#############################
#$(OUTDIR)/outputs_listing:$(TREE)  # changed 29/8/08 to prevent re-running merger.mk
$(OUTDIR)/outputs_listing:$(TREE_FILE)
	rm -fv $(OUTDIR)/*.thin_edges.* # make sure that there are no half baked files 
	@echo "---------------------------- start   running merger.mk ------------------------------"
	$(MAKE) 'TREE_FILE = $(TREE)' -f $(SCRIPTS_DIR)/merger.mk
	@echo "---------------------------- finish  running merger.mk ------------------------------"


#############################
#############################
$(OUTDIR):
	mkdir $(OUTDIR)


#############################
#############################
verbose:
	@echo inputs to clustering process are: $(INPUT_EDGES_FILES)

###########################################################################
##########################################################################
#########################################################################

