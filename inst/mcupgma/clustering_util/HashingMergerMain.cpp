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
 * file: HashingMergerMain.cpp
 *
 * Yaniv Loewenstein, 2007
 *
 * Reads edges from specified input, and transforms cluster
 * indices from previous round, into valid cluster indices.
 * The input edges may, or may not, be singleton edges, but
 * they are regarded as thin edges, relative to the new clusters
 * edges which are "thicker".
 *
 * The output is broken into multiple output files based on a hash
 * function on the new cluster indices. This assures that all input
 * edges that are components of the same new "thick" edge are sent
 * to the  * the same output file, corresponding to the same merging
 * process, in parallel edge merging.
 *
 * All edges might be sent to a single output file if the command
 * line arguemnts do not specify otherwise
 *
 **********************************************************************/
#include "Tree.hpp"
#include "HashingMergerMain.cmdline.h"
#include "ClusteringUtil.hpp"
#include <cstdlib> //gengetopt

using namespace std;


/* ouput filenames go from A.1 .. A.N instead of A.0 .. A.N-1 */
#define HASH_FN_INDEX_FROM_1


/**
 * alloc_ofile_vector
 *
 * allocates a vector of file output streams, opens them, and validates the ostream validity.
 */
void alloc_ofile_vector(string const& base_fn, ofstream_ptr_vec& v){
   char buf[10];
   for (unsigned int k = 0 ; k < v.size() ; k++){
#ifdef HASH_FN_INDEX_FROM_1
      sprintf(buf,".%d",k+1);
#else
      
#endif
      string fn = base_fn + string(buf);
      v[k] = new ofstream(fn.c_str());
      if (!v[k]->is_open())
	 die(string("couldn't open:")+fn);
   }
}

/**
 * dealloc_ofile_vector
 *
 * deallocates the vector of file output streams, closes the streams, and explicitly calls destructors.
 */
void dealloc_ofile_vector(ofstream_ptr_vec& v){
   for (unsigned int k = 0 ; k < v.size() ; k++){
      v[k]->close();
      delete v[k];
      v[k] = NULL;
   }
}
/**
 * main
 *
 * See general description in header above
 */
int main(int argc , char ** argv){

   gengetopt_args_info args_info;
   if (cmdline_parser (argc, argv, &args_info) != 0){
      cerr << "illegal options usage" <<endl;
      exit(1) ;
   }
   
   const dist_t MISSING_EDGE = args_info.missing_val_arg;
   if (MISSING_EDGE < 0)
      die("missing edge value should be non negative");
   if (VERBOSE)
      clog << "Using " << args_info.missing_val_arg<< " for " << "--missing_val option." << endl ;

   string new_edges_filename = args_info.new_edges_file_arg; 
   string old_edges_filename = args_info.existing_edges_file_arg;

   
   //parse the tree into parents and sizes vectors:
   Tree m(MISSING_EDGE); 
   m.init(args_info.tree_file_arg);

   EdgePtrReader er(cin);

   /** there's only a single output file for "old edges", i.e. edges that have not been updated
       after this round of clustering */
   ofstream old_edges_ostream(old_edges_filename.c_str());
   if (!old_edges_ostream.is_open()) die(string("couldn't open for write:")+ old_edges_filename);   

   /* the actual input reading of "thin" edges is done by nodes2parents */
   if (args_info.num_splits_given){
      /***********************************************************
       * Case 1 : Output is broken into multiple output files,
       *   using the hashing version of nodes2parents() and a
       *   vector of ofstreams (used for parallel edge merging processes)
       ***********************************************************/
      const unsigned int N_OFILES = args_info.num_splits_arg;
      if (N_OFILES  <=1 || N_OFILES > 200)
	 die("number of output files to split to - N, should be 1<N<200");
      ofstream_ptr_vec v(N_OFILES);
      alloc_ofile_vector(new_edges_filename,v);
      nodes2parents(m.parents,m.sizes,er,old_edges_ostream,v);
      dealloc_ofile_vector(v);
   } else {
      /************************************************************
       * Case 2 : Output streamed into a single file, and there is no hashing of edges
       *   Into specific output files
       ************************************************************/
      ofstream new_edges_ostream;
      new_edges_ostream.open(new_edges_filename.c_str());
      if (!new_edges_ostream.is_open())
	 die(string("couldn't open for write:")+ new_edges_filename);
      nodes2parents(m.parents,m.sizes,er,old_edges_ostream,new_edges_ostream);
      new_edges_ostream.close();
   }
   old_edges_ostream.close();

   /************************************************
    * dump clusters sizes into given argument file
    ************************************************/
   if (args_info.sizes_file_given){
      ofstream osizes(args_info.sizes_file_arg);
      if (!osizes.is_open())
	 die(string("can't open for write file:")+args_info.sizes_file_arg);      
      dump_container<vector<uint>::iterator, true>(m.sizes.begin(),m.sizes.end(),osizes,m.sizes.size());
      osizes.close();
   }
   return 0;

}
