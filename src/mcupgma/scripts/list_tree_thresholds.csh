#!/bin/csh
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


#foreach d (iteration_? iteration_?? iteration_???)

tail -n1 iteration_?/tree* iteration_??/tree* iteration_???/tree* | removeBlankLines.pl \
	| awk '{if (/=/) { printf "%s\tlast merge:", $2;} else {print "E = "$3,"\tID = "$4;}}'

