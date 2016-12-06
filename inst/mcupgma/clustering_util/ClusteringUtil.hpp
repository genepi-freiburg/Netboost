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
 * file: ClusteringUtil.hpp
 *
 * Yaniv Loewenstein, 2007
 *
 * Static utility functions for  external merging of edges, debugging, parsing etc.
 ***********************************************************************/

/**
 * Revisions
 *
 * Aug 2008 - lonshy added PATH_COMPRESSION option which did not seem to improve the running time signficantly.
 * 
 */


#ifndef _CLUSTERINGUTIL_HPP
#define _CLUSTERINGUTIL_HPP



/**
 * Use the path compression heuristic to accelerate current cluster ID lookup from merged clsuter IDs.
 * finding in large clusters, will otherwise require climbing multiple levels, often repeating climb paths multiple times.
 * The path compression heuristic updates the parent of all nodes to in the path to be the current valid cluster, thus
 * next search along this path will reach the top in O(1). See disjoint-set data-structures or union-find operetion
 * in Cormen et al., Introduction to Algorithms, 2nd edition.
 *
 * The tree is updated when this macro is toggled, and thus it is passed by value, not by reference.
 */
//#define PATH_COMPRESSION



//#define ALIVE_LOOKUP_SPEEDUP



#include <cassert>
#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include "Definitions.hpp"

#include <map>
//#include <utility>
#include <iostream>


static const dist_t default_dist = 100.0;

#define die(x)      { std::cerr << x << " - in " << __FILE__ << ":" << __LINE__ << std::endl; exit(1);}
#define die2(x,y)   { std::cerr << x << y << " - in " << __FILE__ << ":" << __LINE__ << std::endl; exit(1);}
#define die3(x,y,z) { std::cerr << x << y << z << " - in " << __FILE__ << ":" << __LINE__ << std::endl; exit(1);}
#define die_exit(x,ex)      { std::cerr << x << " - in " << __FILE__ << ":" << __LINE__ << std::endl; exit(ex);}
#define die2_exit(x,y,ex)   { std::cerr << x << y << " - in " << __FILE__ << ":" << __LINE__ << std::endl; exit(ex);}
#define die3_exit(x,y,z,ex) { std::cerr << x << y << z << " - in " << __FILE__ << ":" << __LINE__ << std::endl; exit(ex);}


#define PRINT_SPARSITY_INFO

using std::clog;
using std::endl;

typedef std::vector<std::ofstream *>       ofstream_ptr_vec;
typedef std::set<WeightedEdge, IJ_LT_Edge> WeightedEdgesSet;
//typedef __gnu_cxx ::hash_set<WeightedEdge, const IJ_Hasher, IJ_equals> WeightedEdgesSet;



/*******************************************************************************************
 ********** Some general utilities *********************************************************
 *******************************************************************************************/



/**
 * run_or_die()
 *
 * Try to run a command, and die if it fails
 */
void
run_or_die(std::string const & command);


/**
 * MAX(x,y)
 */
template<typename T>
inline
T const &
MAX(T const& x,T const& y){
   return (x > y) ? x : y;
}



/**
 * read_vector()
 *
 * useful for debugging, and loading sparse data
 *
 * reads integer indexed input lines
 * resize by one to get good memory performance, on expense of run time
 */

template <typename T>
void
read_vector(istream& in,std::vector<T>& v){
   //   if (VERBOSE) clog << "read_vector(..) before: v.size() == " << v.size() << endl;
   //#if VERBOSE > 1
   //   uint linecount = 0;
   //#endif
   unsigned int ind;
   T t;
   while (in >> ind >> t){      
      //#if VERBOSE > 1
      //      if (!(++linecount % 10000)) clog  << "read_vector(..) read " << linecount << "lines" <<endl;
      //#endif
      if (v.size() <= ind){
         v.resize(ind+1);
      } 
      v[ind] = t;
   }   
   //   if (VERBOSE) clog << "read_vector(..) after : v.size() == " << v.size() << endl;
}




/**
 * dump_container
 *
 * Print first N entries of an STL container (or array) into a given stream
 */

template<typename ConstInputItr,const bool IS_PRINT_INDEX >
inline
void
dump_container(ConstInputItr start, ConstInputItr end,  std::ostream & os, uint n/* = t.size()*/)
{
   for (uint i=0; i<n  && start!=end ; i++,++start){
      if (IS_PRINT_INDEX) { os << i << '\t'; }
      os << *start << std::endl;
   }
}



/*******************************************************************************************
 ********* External merging related utilities **********************************************
 *******************************************************************************************/


/***************************************
 *** class ClusteringUtil **************
 ***************************************/

class ClusteringUtil {
public:

   /**
    * ClusteringUtil::nodes2parents
    *
    * replaces i,j in the input edges with their parents if they were previously merged
    * this version is non static, and dumps all output into a single stream in contrast
    * to the other two nodes2parents methods.
    */
   static void nodes2parents(index_t * parents,EdgePtrReader& er,std::ostream& os);
};





/**
 * nodes2parents - single ostream version
 *
 * Everything is dumped into one ostream, no hashing.
 */


template<typename T,typename T2 >
inline
void 
nodes2parents(
#ifdef PATH_COMPRESSION  //make a local copy of the tree which is modified by path compression 
	      T parents,
#else
	      T const& parents,
#endif	           
	      T2 const& sizes,EdgePtrReader& er,std::ostream& old_edges_stream,std::ostream& new_edges_stream){

// #ifdef PATH_COMPRESSION
//    clog << "Using path compression in " << __FILE__ << ":" << __LINE__ << endl;
// #else
//    clog << "Not using path compression in " << __FILE__ << ":" << __LINE__ << endl;
// #endif


#ifdef ALIVE_LOOKUP_SPEEDUP //assumes T (the parents container type) is not an array, but has the size() method.
   vector<bool> alive(parents.size());
   for (unsigned int cur = 0; cur < parents.size() ; cur++)
      alive[cur] = (parents[cur] == 0);
#endif
   
   IJ_equals ij_comparator;
   Edge * e = NULL;
   index_t i,j;
   while(er.has_next()){	 
      e = er.next();
      assert(e!=NULL);

#ifdef ALIVE_LOOKUP_SPEEDUP
      if (alive[e->i()] && alive[e->j()]){
	 old_edges_stream << *e << std::endl;
	 delete e;
	 continue;
      }
#endif
      
      //      if (VERBOSE)  clog << ">>\t from:" << *e ;
      i = e->i();
      j = e->j();
      while(parents[i] > 0){
	 i = parents[i];
      }
      while(parents[j] > 0){
	 j = parents[j];
      }
#ifdef PATH_COMPRESSION
      index_t k,next_k;
      k = e->i();
      while(k != i) {
	 next_k = parents[k];
	 parents[k] = i;
	 k = next_k;
      }
      k = e->j();
      while(k != j) {
	 next_k = parents[k];
	 parents[k] = j;
	 k = next_k;
      }
#endif
      if (i!=j){
	 //Edge e_out(i,j,e->dist());
	 WeightedEdge e_out(Edge(i,j,e->dist()),sizes[e->i()]*sizes[e->j()]);
	 if (ij_comparator(e,&e_out)){ //equals 
	    old_edges_stream << e_out << std::endl;
	    //	    if (VERBOSE)  clog << "\tto: "<<e_out << " OLD" <<endl;
	 } else {
	    new_edges_stream << e_out << std::endl;
	    //	    if (VERBOSE)  clog << "\tto: "<<e_out << " NEW" <<endl;
	 }

      } else {
	 //	 if (VERBOSE)  clog << "\tself" << endl;
      }
      
      delete e;
   }
}

/**
 * nodes2parents - multiple ofstreams version (takes an ofstream_ptr_vec)
 *
 * In this version, a vector of multiple ofstrems (output file streams) is taken as an argument instead
 * of one specific ostream.
 * Let N be the number of ofstreams (i.e. N = new_edges_stream.size())
 *
 * Each edge is sent to a different ofstream based on a hashing function h(i,j) --> {0..N-1}
 * Idealy, h is selected to assure uniform hashing, so that output file sizes are balanced, such
 * that the load is balanced on different processing units.
 *
 * here we use: h(i,j) = (i+j)mod(N) ,
 * This is an ad-hoc hashing function which has shown uniform hashing qualities empirically.
 * Ofcourse if i,j are adversrly design to pitfall this function, it will fail, since it does
 * not have any of the theoretical guarantees associated with uniform hashing functions.
 */


template<typename T,typename T2>
inline
void
nodes2parents(
#ifdef PATH_COMPRESSION  //make a local copy of the tree which is modified by path compression
	      T parents,
#else
	      T const& parents,
#endif
	      T2 const & sizes,EdgePtrReader& er,std::ostream& old_edges_stream,ofstream_ptr_vec& new_edges_streams){

// #ifdef PATH_COMPRESSION
//    clog << "Using path compression in " << __FILE__ << ":" << __LINE__ << endl;
// #else
//    clog << "Not using path compression in " << __FILE__ << ":" << __LINE__ << endl;
// #endif

#ifdef ALIVE_LOOKUP_SPEEDUP //assumes T (the parents container type) is not an array, but has the size() method.
   vector<bool> alive(parents.size());
   for (unsigned int cur = 0; cur < parents.size() ; cur++)
      alive[cur] = (parents[cur] == 0);
#endif
   
   const unsigned int N = new_edges_streams.size();
   IJ_equals ij_comparator;
   Edge * e = NULL;
   index_t i,j;
   while(er.has_next()){	 
      e = er.next();
      assert(e!=NULL);


#ifdef ALIVE_LOOKUP_SPEEDUP
      if (alive[e->i()] && alive[e->j()]){
	 old_edges_stream << *e << std::endl;
	 delete e;
	 continue;
      }
#endif
      
      
      //      if (VERBOSE)  clog << ">>\t from:" << *e ;
      i = e->i();
      j = e->j();
      while(parents[i] > 0){
	 //	 std::clog << "parent of i="<<i << " is " << parents[e->i]<<std::endl;
	 i = parents[i];
      }
      while(parents[j] > 0){
	 //	 std::clog << "parent of j="<<i << " is " << parents[e->j]<<std::endl;
	 j = parents[j];
      }
#ifdef PATH_COMPRESSION
      index_t k,next_k;
      k = e->i();
      while(k != i) {
	 next_k = parents[k];
	 parents[k] = i;
	 k = next_k;
      }
      k = e->j();
      while(k != j) {
	 next_k = parents[k];
	 parents[k] = j;
	 k = next_k;
      }
#endif

      if (i!=j){
	 //now i,j are the up-to-date indices, and e->i(),e->j() are the input cluster indices.
	 WeightedEdge e_out(Edge(i,j,e->dist()),sizes[e->i()]*sizes[e->j()]);
	 
	 if (ij_comparator(e,&e_out)){ //equals 
	    old_edges_stream << e_out << std::endl;
	    //	    if (VERBOSE)  clog << "\tto: "<<e_out << " OLD" <<endl;
	 } else {
	    /** we use the hash function H(i,j,N) = (i+j)mod(N), to find the matching
		output file (stream), where (the locals) i,j are the up-to-date cluster indices. */
	    *new_edges_streams[(i+j)%N] << e_out << std::endl;
	    //	    if (VERBOSE)  clog << "\tto: "<<e_out << " NEW" <<endl;
	 }

      } else {
	 //	 if (VERBOSE)  clog << "\tself" << endl;
      }
      
      delete e;
   }
}



/**
 * Replace multiple edges between two clusters created by nodes2parents
 * by a single average linkage edge, including default weight for missing edges.
 *
 * Assumes input is sorted by i,j.
 * T needs to have operator[](int), will work with vector or array
 */ 

template<typename T>
inline
void
average_edges(T counts, EdgePtrReader& er,std::ostream& os){
   if (!er.has_next()) return;
   Edge * e = NULL;
   index_t i= er.peak_next()->i();
   index_t j= er.peak_next()->j();
   dist_t d = 0;
   unsigned int n_edges = 0;
   //      unsigned int expected_n_edges;
   while(er.has_next()){	 
      e = er.next();	 
      assert(e!=NULL);
      std::clog << "read:"<<*e;
      if  (i == e->i() && j == e->j()){
	 /*** case 1 - another edge between the same clusters as previous
	      increase the count and total (comulative) dist */
	 d += e->dist();
	 n_edges++;
	    std::clog << std::endl;
      } else {
	 assert(counts[i]);
	 assert(counts[j]);
	 /*** case 2 - new cluster IDs in edge - flush the previous one
	      and initialize locals accoridingly. */
	 
	 //add the missing edges with default weight:
	 assert(counts[i] * counts[j] >= n_edges);
	 d += ((dist_t)(counts[i] * counts[j] - n_edges)) * default_dist;

	 //output one edge with average weight:
	 os << Edge(i,j,d/(dist_t)(counts[i] * counts[j])) << std::endl;
	    
	 //std::clog << "\tcount["<<i<<"]="<<counts[i]<<", counts["<<j<<"]="<<counts[j]<<", d="<<d;
	 //std::clog << "\toutputing:" << Edge(i,j,d/(dist_t)(counts[i] * counts[j])) << std::endl;
	 
	 //init:
	 j = e->j();
	 i = e->i();
	 d = e->dist();
	 n_edges = 1;
      }
      //	 std::clog << "\tdelete:"<<*e << std::endl;
      delete e;
   }      
}










/**
 * takes a stream of edges, some with duplicate I,J coordinates
 * and fills a set of WeightedEdge's (indexed by I,J) holding the comulative distance of all edges with I, J
 * to the number of such edges, based on WeightedEdge '+' operator.
 */

void
collate_edges(WeightedEdgeReader& wer/*,T const & sizes*/ ,WeightedEdgesSet & weSet) ;


/**
 * thicken_edges
 *
 * takes the set of edges with counts (i.e. WeightedEdges objects) and completes
 * the missing edges with MISSING_VAL to compute average-linkage edges connecting
 * cluster pairs, which are possibly only partially connected. 
 */
template <typename T>
void
thicken_edges(T const & sizes ,WeightedEdgesSet & weSet,dist_t MISSING_VAL) {
   //   if (VERBOSE) { clog << "thicken_edges - set starts with " <<weSet.size() << " edges" <<endl; }
   std::vector<WeightedEdgesSet::value_type> to;
   WeightedEdgesSet::iterator iend  = weSet.end();
   
   for (WeightedEdgesSet::iterator i= weSet.begin(); i!= iend; ){
      to.push_back(*i);
      //i++;
      weSet.erase(i++);
   }
   //   if (VERBOSE) { clog << "thicken_edges - set ends with " <<weSet.size() << " edges" <<endl;}
   //   if (VERBOSE) { clog << "thicken_edges - vector ends with " <<to.size() << " edges" <<endl;}
   count_t expected_n_edges = 0;
   for (std::vector<WeightedEdgesSet::value_type>::iterator to_itr = to.begin(); to_itr != to.end(); ++to_itr){
      expected_n_edges = sizes[to_itr->i()]*sizes[to_itr->j()]; 
      //we print the 'averge edge' which is adding N_missing_edges * MISSING_VAL to the comulative distance and
      // then normalize it all 
      cout << Edge(to_itr->i(),to_itr->j(),
		   (
		    (expected_n_edges - to_itr->count()) * MISSING_VAL
		    +
		    to_itr->dist()
		    )    / expected_n_edges)
	 //the following no longer holds since count holds missing values as well 
	 //#ifdef PRINT_SPARSITY_INFO
	 //	   << '\t' << "had " << to_itr->count() << " edges out of expected " << expected_n_edges << " (with comulative distance:"<< to_itr->dist() << ")"
	 //#endif
	   << endl;

   }    
}






#endif // _CLUSTERINGUTIL_HPP
   
   
   
