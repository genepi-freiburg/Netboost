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

export SHELL := /bin/bash

ifndef DEFINITIONS_MK_INCLUDED
DEFINITIONS_MK_INCLUDED := 1


# turn off, if you want the intermediates to be not deleted
export DELETE_INTERMEDIATE_EDGE_FILES := 1

# Turn on to print some variable values when debugging.
# export VERBOSE :=1

# to be defined from command line, if mosix (or other grid env like rsh) is to be used for external merging code.
# export USE_GRID := 1

# This variable is needed to run recursive make processes with the -f file argument which does not 
# rely on the make include path to look for the file. 
ifndef SCRIPTS_DIR
export SCRIPTS_DIR := .
endif


# If a clustering round is done in less than one second (for trivially small inputs, e.g. in test cases), 
# the timestamps of successive parts might be identical on whole seconds resolution, which will cause make to 
# erroniously skip necessary targets, and omitting necessary tree_partX files from the overall cumulative tree. 
# Setting this variable to X > 0, adds up an extra fake slowdown of X seconds after running the clustering 
# to resolve low time resolution of the touch command when updating timestamps of tree parts. 
# You can set this variable to 0, to prevent needless sleeping.  
export SLEEP_SECONDS_AFTER_CLUSTERING := 0

##################################
# grid jobs prefix - you can use rsh or what have you, just be aware that you are going to do a *lot* of IO 
#      between grid nodes. If you are using a Mosix grid, you can use mospipe to circumvent that.
#      Alternatively, you can modify the executables to read directly from multiple gz files, instead of using
#      pipes which are IO expensive, when using a grid.
ifdef USE_GRID
export RUN := mosrun -b -G -m1000 -M$$PWD
else
export RUN := 
endif
# which executable to use (e.g. AMD64bit on a Linux kernel)
export ARCH := $(shell uname -m)_$(shell uname -s)

##################################
# location of the directory with executables for external merging code:
#export INSTALL_PATH := ~lonshy/mc_upgma_distribution/src
include install_path.mk

export MERGER_APP_BIN := $(INSTALL_PATH)/clustering_util/bin.$(ARCH)
export CLUSTERER_EXEC := $(INSTALL_PATH)/clustering_round/bin.$(ARCH)/hierarchical_clustering

ifdef ALLOW_INEXACT_MERGES 
CLUSTERER_EXEC +=  --allow-non-dendrogram
endif

# use unix's seq command, or a perl imitation if it is not installed
SEQ := $(shell test -x /usr/bin/seq && echo seq || echo seq.pl)


# sh pipe exit code forwarding workaround: 
# Normally, in sh shell when one program in the pipeline exits abnormally, you would not know, since you only 
# get the exit code of the last program. This macro will run an sh command to check if some program in the 
# pipeline has exited abnormally. 
# It checks the exit status of all piped programs in the last invocation of sh commandline in foreground
# by using the sh PIPESTATUS array variable.
#
# double $'s are used to escape $'s from makefile interpretation
# sh doesn't like two consecutive ';'  marks, so don't use one after this macro.
#
# example for usage:
# echo hello | tr e a | exit101 | tr a e  > $@.tmp  ; $(CHECK_PIPE_EXIT_CODES) mv -vf $@.tmp $@
# where exit1 is:
# 1:#!/bin/sh
# 2:cat
# 3:exit 101
# will fail with a 101 exit code and will not reach the mv statement

define CHECK_PIPE_EXIT_CODES_WITH_VERBOSE
for r in $${PIPESTATUS[*]}; do     if (($$r != 0)); then        echo pipe failed;        exit $$r;    else        
echo pipe success;    fi;done; 
endef

define CHECK_PIPE_EXIT_CODES
for r in $${PIPESTATUS[*]}; do if (($$r != 0)); then exit $$r; fi;done; 
endef

endif #DEFINITIONS_MK_INCLUDED
