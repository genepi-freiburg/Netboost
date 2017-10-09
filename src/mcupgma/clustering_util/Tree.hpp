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

/************************************************************************
 *
 * file: Tree.hpp
 *
 * Yaniv Loewenstein, 2007
 *
 * Binary clustering tree data-structure, with pointers to parents
 *  and cluster sizes. 
 *
 ***********************************************************************/

#include <ext/hash_map>
#include <vector>
#include <iostream>
#include "Edge.hpp"

#ifndef _TREE_HPP
#define _TREE_HPP
using namespace std;


class Tree {
private:
   //static const dist_t MISSING_EDGE_WEIGHT = 99999;
   typedef __gnu_cxx::hash_map<Edge *, unsigned int, IJ_hasher, IJ_equals> e2uint;

   const dist_t MISSING_EDGE_WEIGHT;
public:
   vector<index_t>       parents;
   vector<unsigned int>  sizes;
private:

   bool was_init;
   
   void parse_tree(istream&  tree_stream );
   void parse_sizes(istream& sizes_stream);
   //  void new_edges(istream & in,vector<Edge *> & res){   }
public:
   Tree(dist_t MISSING_EDGE_WEIGHT, index_t N = 0) : MISSING_EDGE_WEIGHT(MISSING_EDGE_WEIGHT) , parents(N), sizes(N) , was_init(false){ ; }
      //   vector<index_t> parents = vector<index_t>(N,0);
      //vector<uint > sizes     = vector<uint>(N,0);
   void init(string tree_file);
   void init_sizes(string sizes_file);
   void valid_edges(istream& edge_stream, ostream& live_edges,ostream& new_edges);
   void verbose();
};
#endif //_TREE_HPP
