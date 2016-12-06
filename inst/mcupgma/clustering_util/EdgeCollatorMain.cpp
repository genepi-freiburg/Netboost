/***************************************************************************
 *                                                                         *
 * MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA  *
 *             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9. *
 *									   *
 * Copyright (C) Yaniv Loewenstein				 	   *
 *               School of Computer Science And Engineering		   *
 *               Hebrew University of Jerusalem				   *
 *									   *
 *     All Rights Reserved						   *
 *									   *
 *     This source code is distributed under the terms of the		   *
 *     GNU General Public License. See the file LICENSE			   *
 *     for details.                                                        *
 *									   *
 ***************************************************************************/

/**********************************************************************
 * file: EdgeCollatorMain.cpp
 *
 * Yaniv Loewenstein, 2007
 *
 * A main for an executable that reads edges which are components in "thicker"
 * average edges connecting larger up-to-date clusters. The output is the up-to-date
 * average edges between connected clusters. Input contains edges
 * with cluster indices modified to up-to-date cluster IDs, and the weight of this
 * component in the sum.
 *
 * It is necessary that all relevant components to a pair of clusters
 * are processed by the same instance of this program, to allow correct average computation,
 * since unread components are assumed to be missing in the input graph, and are substituted
 * by the missing edge value (psi).
 *
 ***********************************************************************/


#include "ClusteringUtil.hpp"
#include "EdgeCollatorMain.cmdline.h"
#include "Tree.hpp"


using namespace std;
int main(int argc , char ** argv){
   gengetopt_args_info args_info;
   if (cmdline_parser (argc, argv, &args_info) != 0){
      die("illegal options usage");      
   }
   Tree tree(args_info.missing_val_arg);
   if (args_info.sizes_file_given){
      if (args_info.tree_file_given){
         die("Please supply either tree or sizes file argument (use -help) not both");
      }
      tree.init_sizes(args_info.sizes_file_arg);
   } else {
      if (!args_info.tree_file_given){
         die("Please supply either tree or sizes file argument (use -help)");
      }
      tree.init(args_info.tree_file_arg);
   }
   WeightedEdgesSet weSet;
   //   if (VERBOSE) { clog << "c-tor WeightedEdgeReader start" <<endl;}
   WeightedEdgeReader wer(cin);
   //   if (VERBOSE) { clog << "c-tor WeightedEdgeReader end" <<endl; }
   //   if (VERBOSE) { clog << "calling collate_edges(..)" <<endl; }
   
   /****************************************
    * Read actual weighted edges input, where all components of a specific thick edge are summed into
    * a single thick edge element in a set indexed by current cluster IDs
    ****************************************/
try{
   collate_edges(/*WeightedEdgeReader(cin)*/wer,weSet);

   //   if (VERBOSE) {   clog << "fin collate_edges(..)" <<endl; }
   //   if (VERBOSE) {   clog << "we have a total of " << weSet.size() << " thick valid edges" << endl;}
   //   if (VERBOSE) {   clog << "calling thicken_edges(..)" <<endl;   }

   /****************************************
    * Add missing edges to allow correct averaging for
    * sparse inputs.
    ****************************************/
   thicken_edges(tree.sizes,weSet,args_info.missing_val_arg);

   //   if (VERBOSE) {clog << "fin thicken_edges(..)" <<endl; }
} catch (std::bad_alloc e) {
	 cerr << "caught std::bad_alloc" << endl;
	 cerr << e.what() << endl;
	 die_exit("FATAL ERROR: Seems like we're out of memory, try to split the input into smaller chunks.",201);	
}
   return 0;
}	
