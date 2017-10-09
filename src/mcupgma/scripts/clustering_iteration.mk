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


##########################################################################
##   NUMERIC PARAMS  #####################################################
##########################################################################
include definitions.mk


#################################################################
# maximum heap size (= M edges) for clustering, depends on the  #
# available memory-size. In general, the larger M the more      #
# merges done per clustering round and the faster the complete  #
# process. 							#
################################################################
M = 10000000

#################################################################
# missing edge weight					   	#
################################################################
PSI = 1000

################################################################
# must be larger than any ID in the input		       #
# TODO: find this automatically if not given by user	       #
###############################################################
ifdef ISIZES_FILE 
FIRST_MERGE_ID = $(shell tail -1 $(ISIZES_FILE) | awk '{print $$1}')
CLUSTERER_FLAGS += --input-cluster-sizes $(ISIZES_FILE)
ifeq ($(FIRST_MERGE_ID),)
$(error cant retrieve FIRST_MERGE_ID from file $(ISIZES_FILE))
endif
else
CLUSTERER_FLAGS += 
ifndef FIRST_MERGE_ID
FIRST_MERGE_ID = 5000000
$(warning first merge ID not defined setting to $(FIRST_MERGE_ID))
endif # not def FIRST_MERGE_ID	
endif # ISIZES_FILE

#filt = head -$(M)
ifndef filt
filt = cat
endif

#ifndef OSIZES_FILE
OSIZES_FILE = $(OUTDIR)/cluster_sizes
#endif

ifndef TREE_FILE
TREE_FILE =   $(OUTDIR)/tree
endif

ifeq ($(strip $(OUTDIR)),)
$(error $$OUTDIR not defined)
endif

ifndef ARCH
ARCH := $(shell uname -m)_$(shell uname -s)
endif


# clustering executable with flags and arguments set
CLUSTERER = $(CLUSTERER_EXEC) --average-type=arithmetic --number-of-input-edges $(M)  --max-distance $(PSI) --max-cluster-index $(FIRST_MERGE_ID) 

##############################################
#             TARGETS                        #
##############################################


$(CUMULATIVE_TREE_FILE):$(TREE_FILE)

#####################################################################
#   Run a round of clustering using the clusterer executable        #
#####################################################################

$(TREE_FILE):$(INPUT_EDGES_FILES)
	@echo making $@ from $^
	smart_cat $(INPUT_EDGES_FILES) 	| $(filt) | $(CLUSTERER)  --output-merges $@.tmp  --output-cluster-sizes $(OSIZES_FILE).tmp $(CLUSTERER_FLAGS) && \
	mv -vf $@.tmp $@ && \
	mv -vf $(OSIZES_FILE).tmp $(OSIZES_FILE)
	@test -s $@ || (echo "*** ERROR *** : seems like tree file: $@  is empty (no clustering done in this round). This might happen erroniously if your input includes edges >= $(PSI) (\psi - the user supplied parameter for maximum dissimilarity)." >> /dev/stderr ; exit 103)
