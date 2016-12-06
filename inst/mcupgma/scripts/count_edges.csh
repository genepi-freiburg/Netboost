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

# Usage: records how many edges remain in *.edges.gz files per iteration_X folder into iteration_X.count

foreach dir (iteration_? iteration_?? iteration_???)
	test -s $dir.count || scat $dir/*.edges.gz | wc -l >>! $dir.count
end
