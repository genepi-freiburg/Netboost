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
 * file: Tree.cpp
 *
 * Yaniv Loewenstein, 2007
 *
 * See Tree.hpp
 *
  ***********************************************************************/



//#include "EdgeReader.hpp"
#include "Tree.hpp"
#include "ClusteringUtil.hpp"
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
using namespace std;
//#define MAX(x,y) {  (((x) > (y)) ? (x) : (y)) }
#define edge_count second
//#define VERBOSE 2 //uses ClusteringUtil::Tree

/**
 * inits members: parents, sizes
 *   from tree
 */


#define MAX_LINE_LENGTH 1000
void
Tree::parse_tree(istream& tree_stream){
   assert(parents.size() == sizes.size());
   index_t i,j,k;
   unsigned int size;
   dist_t d;
   char buf[MAX_LINE_LENGTH];
   while(tree_stream.getline(buf,MAX_LINE_LENGTH)){
      assert(strlen(buf) > 0);
      istringstream line_stream(buf,istringstream::out);
#ifdef TREE_HAS_SIZES
      line_stream  >> i >> j >> d >> k >> size;
      if (!line_stream)
	 die(string("seems like the input doesn't have 5 records per line:\"")+buf+string("\""));
#else
      line_stream  >> i >> j >> d >> k;
      if (!line_stream)
	 die(string("seems like the input doesn't have 4 records per line:\"")+buf+string("\""));	    
#endif

      //while(tree_stream  >> i >> j >> d >> k >> size){
      if (VERBOSE >1) clog << "i="<<i<<", j="<<j<<", d="<<d<<", k="<<k<<", size="<<size<<endl;
      assert(k > i && k > j);
      if (k >= parents.size()){//if (i > max || j > max);
	 parents.resize(k+1);

	 sizes.resize(k+1);

      }
      parents[i] = parents[j] = k;

      if (!sizes[i])
	 sizes[i] = 1;
      if (!sizes[j])
	 sizes[j] = 1;
      sizes[k] = sizes[i] + sizes[j];
#ifdef TREE_HAS_SIZES
      assert(sizes[k] == size);
#endif
      
   }
   /*     #ifdef TREE_HAS_SIZES */
   for (unsigned int t=0; t<sizes.size(); t++) 
      if (sizes[t] == 0)
	 sizes[t] = 1;
   /*     #endif */
   
}
/**
 * init members from current clustering tree
 */
void
Tree::init(string tree_file){
   ifstream tree_stream(tree_file.c_str());
   if (!tree_stream.is_open())
      die(string("couldn't open tree file:")+tree_file);
   parse_tree(tree_stream);
   tree_stream.close();
   was_init = true;
   if (VERBOSE) verbose();

}
void
Tree::verbose(){
   clog << "Tree: was_init =?" << was_init << endl;
   clog << "Tree: |parents| = "<<parents.size()<<endl;
   clog << "Tree: |sizes| = "<<sizes.size()<<endl;
   if (!parents.empty())
      clog << "Tree: parents["<<parents.size()-1<<"] = "<<parents[parents.size()-1]<<endl;
   if (!sizes.empty())
	 clog << "Tree: sizes["<<sizes.size()-1<<"] = "<<flush<<sizes[sizes.size()-1]<<endl;
}
void
Tree::parse_sizes(istream& sizes_stream){
   cerr << "parse_sizes(..)" <<endl;
   read_vector(sizes_stream,sizes);
}

void
Tree::init_sizes(string sizes_file){
   ifstream sizes_stream(sizes_file.c_str());
   if (!sizes_stream.is_open())
      die(string("couldn't open sizes file:")+sizes_file);
   parse_sizes(sizes_stream);
   sizes_stream.close();
   was_init = true;
   if (VERBOSE) verbose();
}  

//never tested yet:
void
Tree::valid_edges(istream& edge_stream, ostream& live_edges,ostream& new_edges){
   if (!was_init) 
      die("illegal: must call init(..) before this");
   /******************************************************
    *** Read edges - live (thin) edges to stream
    ***              dead edges accumalte in hash
    ******************************************************/
   EdgePtrReader er = EdgePtrReader(edge_stream);
   e2uint e2count;
   while(er.has_next()){
      Edge * e = er.next();
      assert(parents.size() > MAX(e->i(),e->j()));
      if (parents[e->i()] == 0 && parents[e->j()] == 0){
	 /** Case A - both vertices are valid **/
	 live_edges << e;
      } else {
	 /** Case B - at least one of the vertices was merged */
	 //climb up code here..	 
	 if (e->i() == e->j()){
	    /** Case 0 - an intra-cluster edge - irrelevant as it was already merged */	    	 	 
	 } else {
	    e2uint::iterator itr;
	    if ((itr = e2count.find(e)) == e2count.end()){
	       /** Case 1 - we have'nt seen an edge between these clusters yet **/
	       e2count[e] = 1;
	       continue; //this the only case where we don't want to delete e
	    } else {
	       /** Case 2 - we have seen an edge between these clusters  **/
	       //we cannot modify the key so we erase the element from the hash
	       // call copy c-tor to assign (new) comulative dist and increment count
	       Edge * was = itr->first;
	       itr->edge_count += 1;
	       assert(was->i() == e->i() && was->j() == e->j());
	       //e2count.erase(itr++);
	       *was = Edge(was->i(),was->j(),was->dist()+e->dist());
	       //itr--;
	       //e2count.insert(itr,e2uint::value_type(was,count));	    
	    }	
	 }
      }
      delete e;      
   } //end while(has_next())  
     /*******************************************************
      ***  Calculate thick edges
      *******************************************************/
   vector<Edge *> res = vector<Edge *>();
   res.reserve(e2count.size());
   for (e2uint::iterator itr = e2count.begin(); itr != e2count.end() ; itr++){
      Edge * e   = itr->first;
      uint count = itr->second;
      uint expected_count = sizes[e->i()] * sizes[e->j()];
      *e = Edge(e->i(), e->j(), (e->dist() + (expected_count - count) * MISSING_EDGE_WEIGHT) / expected_count);
      res.push_back(e);
   }
   e2count.clear();
   sort(res.begin(),res.end(),LT_Edge_dist());
   
   for (vector<Edge *>::iterator itr = res.begin() ; itr != res.end() ; itr++){
      new_edges << *itr;
      delete *itr;
   }
   res.clear();
}
