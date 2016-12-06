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

###############################################
##### MC-UPGMA INSTALLER - MAKEFILE ###########
###############################################

DISTRIBUTION_VERSION  := 1.0.0


################################################
####### VARIABLES ##############################
################################################

### Maintain multiple binaries per OS and architecture ####
ARCH := $(shell uname -m)_$(shell uname -s)

SCRIPTS_DIR             := scripts

# directory with clustering code
CLUSTERING_DIR 		:= clustering_round
# directory with clustering executables
CLUSTERING_DIR_BIN      := $(CLUSTERING_DIR)/bin.$(ARCH)

# directory with clustering utilities code
CLUSTERING_UTIL_DIR     := clustering_util
# directory with clustering utilities executables
CLUSTERING_UTIL_DIR_BIN := $(CLUSTERING_UTIL_DIR)/bin.$(ARCH)

EXECUTABLES := $(CLUSTERING_UTIL_DIR_BIN)/edge_collator $(CLUSTERING_UTIL_DIR_BIN)/edges2valid_clusters
EXECUTABLES += $(CLUSTERING_DIR_BIN)/hierarchical_clustering

MCUPGMA_UTIL_SOURCES := ClusteringUtil.hpp ClusteringUtil.cpp Tree.hpp \
  Tree.cpp HashingMergerMain.cpp HashingMergerMain.ggo HashingMergerMain.cmdline.h HashingMergerMain.cmdline.c \
  EdgeCollatorMain.cpp  EdgeCollatorMain.ggo  EdgeCollatorMain.cmdline.h  EdgeCollatorMain.cmdline.c \
  Definitions.hpp Edge.hpp IteratorAdaptor.hpp \
  Makefile LICENSE CONTENTS

MCUPGMA_UTIL_SOURCES := $(addprefix $(CLUSTERING_UTIL_DIR)/,$(MCUPGMA_UTIL_SOURCES))


DISTRIBUTION_FILES    := $(CLUSTERING_DIR).tar.gz $(CLUSTERING_UTIL_DIR).tar.gz scripts/ \
	LICENSE INSTALL README_USER README_DEVELOPER Makefile examples/

DISTRIBUTION_FILENAME := mcupgma_$(DISTRIBUTION_VERSION).tar.gz 


################################################
######## TARGETS ###############################
################################################
all: 
install: install_util install_clustering 
configure: update_path_definitions
test: test_install 
FORCE:

update_path_definitions: $(SCRIPTS_DIR)/install_path.mk

$(SCRIPTS_DIR)/install_path.mk: FORCE
#	@test -e $(SCRIPTS_DIR)/definitions.mk  || (echo "Can't find the file $(SCRIPTS_DIR)/definitions.mk"; exit 12;)
	echo "export INSTALL_PATH := "`pwd` > $@

install_util: 
	@echo "==== extracting $(CLUSTERING_UTIL_DIR) sources (external merging programs)"
	tar -xvzf  $(CLUSTERING_UTIL_DIR).tar.gz
	@echo "==== compiling $(CLUSTERING_UTIL_DIR) sources using external Makefile ($(CLUSTERING_UTIL_DIR)/Makefile) recursively "
	$(MAKE) -C $(CLUSTERING_UTIL_DIR)  -f Makefile mcupgma_install
install_clustering: 
	@echo "==== extracting $(CLUSTERING_DIR) sources (code for a round of MC-UPGMA clustering)"
	tar -xvzf  $(CLUSTERING_DIR).tar.gz
	@echo "==== compiling $(CLUSTERING__DIR) sources using external Makefile ($(CLUSTERING_DIR)/Makefile) recursively "
#	$(MAKE) -C $(CLUSTERING_DIR)  -f Makefile $(CLUSTERING_DIR_BIN)/hierarchical_clustering
	$(MAKE) -C $(CLUSTERING_DIR)  -f Makefile
test_install:
	@echo "==== testing installation - 'executabilty' of necessary programs."
	@for r in $(EXECUTABLES); do test -x $$r || (echo "Cant execute program $r" ; exit 11;);	done;
	@echo "tested o.k. - all programs ($(EXECUTABLES)) are runnable"
tarball_clustering: $(CLUSTERING_DIR).tar.gz
tarball_util: $(CLUSTERING_UTIL_DIR).tar.gz

$(CLUSTERING_DIR).tar.gz: $(wildcard $(CLUSTERING_DIR)/*)
	@echo "==== creating tarball for clustering code $(CLUSTERING_DIR).tar.gz"
	tar -czf $(CLUSTERING_DIR).tar.gz $(CLUSTERING_DIR)           --exclude '*bin*' --exclude '*obj*' --exclude '*dep*' --exclude '*~'
$(CLUSTERING_UTIL_DIR).tar.gz: $(wildcard $(CLUSTERING_UTIL_DIR)/*)
	@echo "==== creating tarball for clustering utilities code $(CLUSTERING_UTIL_DIR).tar.gz"
#	$(MAKE) -C $(CLUSTERING_UTIL_DIR) -f Makefile $(CLUSTERING_UTIL_DIR).tar.gz && mv $(CLUSTERING_UTIL_DIR)/$(CLUSTERING_UTIL_DIR).tar.gz .
#	tar -czf $(CLUSTERING_UTIL_DIR).tar.gz $(CLUSTERING_UTIL_DIR) --exclude '*bin*' --exclude '*obj*' --exclude '*dep*' --exclude '*~'
	tar -czf $(CLUSTERING_UTIL_DIR).tar.gz  $(MCUPGMA_UTIL_SOURCES)

distribution: $(DISTRIBUTION_FILENAME)

$(DISTRIBUTION_FILENAME): $(DISTRIBUTION_FILES)
	@echo "==== creating distribution tarball $@ with scripts, clustering code, clustering utilities code, etc."
#	@test -e $@ && (echo "refusing to overrun existing copy of $@"; exit 12;) || echo ok
	tar -czf $@ $^  --exclude '*bin*' --exclude '*obj*' --exclude '*dep*' --exclude '*~'
