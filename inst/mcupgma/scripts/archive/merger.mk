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

###############################################################################################################
## merger.mk -  External computation of valid average-linkage edges after a round of (partial) clustering    ##
##              valid edges are edges between a pair of yet unmerged clusters.				     ##
##													     ##
##              Currently arithmetic average only is supported.						     ##
##													     ##
## Author:  Yaniv Loewenstein , lonshy@cs.huji.ac.il                                                         ##
##													     ##
## Revised: 08/10/2007 - added comments, removed dead code, added cleanup target + .INTERMEDIATE             ##
##                       moved some code to definitions.mk				                     ##
##                                                   							     ##
###############################################################################################################

#####
# Common errors to watch out for:
#  if the implicit rules seems not to exist (i.e. make claims no way to make thick/thin_edges files) 
#    it may be since it can't find the tree file
# 
# a function from a list of items, to a list of uniq items
#
# assumptions:
# currently this make file requires that all inputs are in the same directory (path)
# to overcome this we need to extract rule #1 for each input directory in INPUT_PATHS
####

include definitions.mk

###############################################################################
##  MY FUNCTIONS     ##########################################################
###############################################################################

# a function that converts a list of words, to a list of unique words
define uniq.template
        $(shell echo "$(1)" | awk '{for (i=1 ; i<= NF ;i++){print $$i;}}' | sort | uniq)
endef



###############################################################################
##  EXTERNAL PARAMS  ##########################################################
###############################################################################
# this params are expected to be defined by the caller, but default values are
# given here for ease of debugging without calling from an external makefile.
ifeq (0,${MAKELEVEL})
############################################################################
# to be overridne externally, usefull for debug now                        #
#                                                                          #
t = 1
OUTDIR = iteration_$t
# order is important here, to preserve sorting if we're using head before clustering
#INPUT_EDGES_FILES    = $(wildcard iteration0/*.edges.gz)

TREE_FILE            = $(OUTDIR)/tree.$(t)
SIZES_FILE           = $(OUTDIR)/osizes.$(t)							
# missing edge weight					   	           #
PSI = 1000								   #
#									   #
############################################################################
else
ifndef OUTDIR
$(error $$OUTDIR should have been set externally)
endif
ifndef TREE_FILE
$(error $$TREE_FILE should have been set externally)
endif
ifndef PSI
$(error $$PSI should have been set externally)
endif

endif # endif MAKELEVEL == 0

# a filename for a text file holding a listing of output files of the merger process ("thick edges")
OUTPUT_FILES_LIST_FILE = $(OUTDIR)/outputs_listing
#########################################################################
# number of hashing buckets - when merging an input file of old edges   #
# the new edges are seperated into this number of "new edges" files	#
########################################################################
H = 2

##########################################################################
##   FILES and EXECUTABLES config  ######################################
########################################################################


ifndef TREE_FILE
$(error TREE_FILE undefined)
endif
####################################
#### EXEC ########################
################################


# executable from previous round input edges to valid cluster IDs (i.e. child replaced by parent cluster)
NODES2CLUSTERS := $(MERGER_APP_BIN)/edges2valid_clusters --missing-val=$(PSI) --tree-file $(TREE_FILE)
# Edge merger - collate thin edges and average into one thick "cluster edge" including missing edges
EDGE_AVERAGER  := $(MERGER_APP_BIN)/edge_collator        --missing-val=$(PSI) --tree-file $(TREE_FILE)

###############################################################
##   LOCAL VARIABLES  - CONTROL MAKE TARGETS #################
#############################################################


# one per input temporary intermediate var to define others:
# ( input files in output dir with no extension)
INPUT_PATHS	    := $(dir $(INPUT_EDGES_FILES))
INPUT_PATHS	    := $(call uniq.template,$(INPUT_PATHS))
ifneq ($(words $(INPUT_PATHS)),1)
$(error this makefile currently assumes all inputs are in the same path, paths are $(INPUT_PATHS))
endif
PATHLESS_INPUT_EDGES_FILES                 := $(notdir $(INPUT_EDGES_FILES))
#BASE_NAMES           = $(patsubst %.edges.gz,$(OUTDIR)/%,$(INPUT_EDGES_FILES))
BASE_NAMES           := $(patsubst %.edges.gz,$(OUTDIR)/%,$(PATHLESS_INPUT_EDGES_FILES))
# these are edges whose clusters weren't clustered, hence they are still
# valid edges and are unchanged -  one per input file.
EXISTING_EDGES_FILES := $(patsubst %, %.existing_edges, $(BASE_NAMES))
# a total of H files per input file - intermediates (can be deleted after outputs are ready)
THIN_EDGES_FILES    := $(foreach base,$(BASE_NAMES),$(foreach h,$(shell $(SEQ) 1 $(H)),$(base).thin_edges.$(h)))
# H output files 
THICK_EDGES_FILES    :=                              $(foreach h,$(shell $(SEQ) 1 $(H)),$(OUTDIR)/$(h).thick.edges.gz)

# two vars containing "%" - together they are right hand side of the 
# implicit rule: `input ==> several interamidates`
ALL_EXIST_PAT := $(OUTDIR)/%.existing_edges
WILDCARD_ALL_Hs_THIN := $(foreach h,$(shell $(SEQ) 1 $(H)),$(OUTDIR)/%.thin_edges.$(h))
THIN_PER_ONE_H_PAT := $(foreach f,$(BASE_NAMES),$(f).thin_edges.%)

EXISTING_EDGES_FILE  := $(OUTDIR)/existing.edges.gz
# everything needed for the next round of clustering

############################################################################
### k-ary split of existing edges ##########################################
#K := 5

ifdef K
ONE_TO_K := $(shell $(SEQ) 1 $(K))
# use basename twice to get rid of .edges.gz (two extensions)
EXISTING_EDGES_FILE_SPLIT := $(foreach i,$(ONE_TO_K),$(basename $(basename $(EXISTING_EDGES_FILE))).$(i).edges.gz)
OUTPUTS                := $(THICK_EDGES_FILES) $(EXISTING_EDGES_FILE_SPLIT)
else
OUTPUTS                := $(THICK_EDGES_FILES) $(EXISTING_EDGES_FILE)
endif
############################################################################


##################################################################################################
########## TARGETS #############################################################################
###############################################################################################

########################################################
##### A L L ############################################
########################################################

#all:  verbose $(OUTPUT_FILES_LIST_FILE)
#all:  verbose $(THIN_EDGES_FILES)

#.EXPORT_ALL_VARIABLES:


### to remove intermidates 

#all:  verbose $(OUTDIR) $(IMPLICIT_RULES_MAKEFILE)
#all:  verbose $(THIN_EDGES_FILES)
#all:  verbose $(THICK_EDGES_FILES)
#all:  verbose $(OUTPUT_FILES_LIST_FILE)
#all:  test_tree_exist verbose $(THIN_EDGES_FILES)  $(OUTPUT_FILES_LIST_FILE) 

all:  test_tree_exist verbose $(OUTPUT_FILES_LIST_FILE) cleanup
	@echo "********************************************************************************"
.DELETE_ON_ERROR:
.PHONY: verbose test_tree_exist

ifdef DELETE_INTERMEDIATE_EDGE_FILES
.INTERMEDIATE: $(THIN_EDGES_FILES) $(EXISTING_EDGES_FILES)
cleanup: $(OUTPUTS)
	@echo "cleaning up the intermediate thin_edges & existing_edges files etc." 
	rm -vf $(THIN_EDGES_FILES) $(EXISTING_EDGES_FILES)
else
.SECONDARY: $(THIN_EDGES_FILES) $(EXISTING_EDGES_FILES)
cleanup: $(OUTPUTS)
	@echo "skipping cleanup since DELETE_INTERMEDIATE_EDGE_FILES is not defined"
endif

test_tree_exist:
# must be with @, otherwise false WARNINGS printing will appear because of echoing the command
	@test -e $(TREE_FILE) || $(shell echo "*** WARNING *** : tree file:$(TREE_FILE) - does not exist")




#######################################################
##### implicit rule #1: make thin edges + existing 
##### From 1 input file (+ tree) to H X thin_edges files + 1 X existing_edges_pattern
#######################################################
$(ALL_EXIST_PAT) $(WILDCARD_ALL_Hs_THIN):$(INPUT_PATHS)/%.edges.gz $(TREE_FILE) 
	@echo "*** creating thin edges files $@ .. $(H) ***"
	smart_cat $< |  $(RUN) $(NODES2CLUSTERS) --existing-edges-file $(OUTDIR)/$*.existing_edges.tmp  --new-edges-file $(OUTDIR)/$*.thin_edges  --num-splits $(H); $(CHECK_PIPE_EXIT_CODES)   mv -vf $(OUTDIR)/$*.existing_edges.tmp  $(OUTDIR)/$*.existing_edges
#######################################################
##### implicit rule #2:
##### from H interamidate thin_edges files with same hash key h to one outputfile
#######################################################
$(OUTDIR)/%.thick.edges.gz:$(THIN_PER_ONE_H_PAT) $(TREE_FILE)
#	smart_cat $(filter-out $(TREE_FILE),$^) | $(EDGE_AVERAGER) | tee $@.delme | awk '{print $$1,$$2,$$3}' | gzip -f > $@.tmp; $(CHECK_PIPE_EXIT_CODES)  mv -vf $@.tmp $@
	smart_cat $(filter-out $(TREE_FILE),$^) | $(RUN) $(EDGE_AVERAGER)         | awk '{print $$1,$$2,$$3}' | gzip -f > $@.tmp; $(CHECK_PIPE_EXIT_CODES)  mv -vf $@.tmp $@


$(OUTPUT_FILES_LIST_FILE): $(OUTPUTS)
	@echo $(OUTPUTS) > $@.tmp && mv -vf $@.tmp $@

$(OUTDIR): 
	mkdir $(OUTDIR)

$(EXISTING_EDGES_FILE): $(EXISTING_EDGES_FILES)
	smart_cat $^ |  awk '{print $$1,$$2,$$3}' | gzip -f > $@.tmp; $(CHECK_PIPE_EXIT_CODES) mv -vf $@.tmp $@

### k-ary split of existing edges
$(EXISTING_EDGES_FILE_SPLIT): $(EXISTING_EDGES_FILES)
	smart_cat $^ |  awk '{print $$1,$$2,$$3}' | file_k_split.pl -k $(K) -ogz -base $(basename $(basename $(EXISTING_EDGES_FILE))) -ext edges; $(CHECK_PIPE_EXIT_CODES)

# if something (such as existing edges remain)
verbose: 
ifdef VERBOSE
	echo verbose
	@echo "OUTPUTS =  $(OUTPUTS)"
	@echo "OUTPUT_FILES_LIST_FILE= $(OUTPUT_FILES_LIST_FILE)"
	@echo "PATHLESS_INPUT_EDGES_FILES = $(PATHLESS_INPUT_EDGES_FILES)"
	@echo "INPUT_PATHS = $(INPUT_PATHS)	"
	@echo "INPUT_EDGES_FILES       = $(INPUT_EDGES_FILES)"
	@echo "EXISTING_EDGES_FILES = $(EXISTING_EDGES_FILES)"
	@echo "EXISTING_EDGES_FILE = $(EXISTING_EDGES_FILE)"
	@echo "WILDCARD_ALL_Hs_THIN        = $(WILDCARD_ALL_Hs_THIN)"
	@echo "THICK_EDGES_FILES       = $(THICK_EDGES_FILES)"
	@echo "THIN_EDGES_FILES       = $(THIN_EDGES_FILES)"
	@echo "rule #1 (create thin edges): $(ALL_EXIST_PAT) $(WILDCARD_ALL_Hs_THIN):$(INPUT_PATHS)/%.edges.gz $(TREE_FILE)"
	@echo "rule #2 (create thick edges): $(OUTDIR)/%.thick.edges.gz:$(THIN_PER_ONE_H_PAT) $(TREE_FILE)"
	@echo "==========="
	@echo "K  = $(K)"
	@echo "ONE_TO_K = $(ONE_TO_K)"
	@echo "EXISTING_EDGES_FILE_SPLIT = $(EXISTING_EDGES_FILE_SPLIT)"
endif

#########################################################################################################
#################################### EOF ################################################################
#########################################################################################################
