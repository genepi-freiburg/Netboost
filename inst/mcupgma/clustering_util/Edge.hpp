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
 * file: Edge.hpp
 *
 * Yaniv Loewenstein, 2007
 *
 * An Edge data-structure (for a sparse graph), I/O operators, and
 * WeightedEdge data-structure to facilitate edge arithmetic
 *
 ***********************************************************************/


#ifndef _EDGE_HPP
#define _EDGE_HPP
#include <iostream>
#include <cfloat> //LDBL_MIN
#include <climits> //UINT_MAX

#include <iterator> //istream_iterator
#include <iomanip>

#include <algorithm> //swap(a,b)
#include <cassert>

#include <ext/hash_set> //for hash<Edge>

/*
  macros:
    FORMATTED_OUTPUT - nicer printing of edges using IO manipulator, possible slower
    INPUT_I_LT_J     - when reading a line from input using>>, will swap(i,j) if j<i
*/

#define FORMATTED_OUTPUT
#define INPUT_I_LT_J


#include <set> //mypos
#define CLUSTER_ID_IOMANIP std::setw(8)
#ifdef FORMATTED_OUTPUT
#define EDGE_HEADER "CLUST_ID_1_(I)\tCLUST_ID_2_(J)\tDISTANCE"
#else
#define EDGE_HEADER "I\tJ\tDist"
#endif

//#define DIST_OUTPUT_N_DIGITS 8

typedef unsigned int  index_t;
typedef unsigned int  count_t;
typedef double        dist_t; //distance , cin fails when this is a float
typedef unsigned char uchar;

const dist_t  INVALID_DIST  = LDBL_MIN;
const index_t INVALID_INDEX = UINT_MAX;
const count_t INVALID_COUNT = UINT_MAX;

using std::ostream;
using std::istream;

//#define SAFE_READ_EDGE //won't compile for now

#ifdef SAFE_READ_EDGE //check for format errors
#include <sstream>
#define MAX_LINE_LENGTH 1000
char linebuf[MAX_LINE_LENGTH];
std::istringstream edge_sstream;
#endif


/**********************************************************************************
 ****************** class Edge   **************************************************
 **********************************************************************************/

class Edge{
   /*** friends for building from an istream or printing to ostream */
   
   friend std::istream& operator>>(std::istream& is,Edge  & e);      //assign to existing edge
   friend std::istream& operator>>(std::istream& is,Edge * & e);     //e assigned new Edge(..) or  NULL 
   //friend class istream_iterator; //hide default c-tor from non friends
   //friend std::ostream& operator<<(std::ostream& os,const Edge& e);  //ouput
protected:
   /****   private members  *********/
   index_t _i,_j;
   dist_t  _dist;
#ifdef FLAGS
   //adds an additional member to the Edge data structure, and the appropriate enum
   uchar _flags;
   enum states { i_approx, j_approx,both_approx,exact};
#endif

   
public:
   /******* C-tors *****/   
   /* Default c-tor - needed for stl istream_iterator<Edge> */
   Edge(): _i(INVALID_INDEX), _j(INVALID_INDEX) , _dist(INVALID_DIST)
#ifdef FLAGS
	 ,_flags(exact)
#endif
   {  ; } 
   Edge(index_t i,index_t j, dist_t dist) : _i(i), _j(j), _dist(dist)
#ifdef FLAGS
					  ,_flags(exact)
#endif
   {
#ifdef INPUT_I_LT_J
      if (_j<_i) std::swap(_i,_j);
#endif   
      ; }
   Edge(Edge const& e)  : _i(e._i), _j(e._j), _dist(e._dist)
#ifdef FLAGS
			,_flags(e._flags)
#endif
      //      std::clog << "Edge copy c-tor" << std::endl;
   { ; }; 
   /***** Member - getters  ******/
   inline index_t i()    const {return _i;}
   inline index_t j()    const {return _j;}
   inline dist_t  dist() const {return _dist;}
   inline bool is_valid() const{return (i() != INVALID_DIST && j() != INVALID_DIST || dist() != INVALID_DIST);}

   /*** Output formatting ***/
   inline ostream & print_edge(ostream& os) const{
#ifdef FORMATTED_OUTPUT
      return os << std::left
		<< CLUSTER_ID_IOMANIP << i() << '\t'
		<< CLUSTER_ID_IOMANIP << j() << '\t'
		<< std::scientific << dist();
#else
      return os << i() << '\t' << j() << '\t' << dist();      
#endif
   }

      
   
   inline void addDist(dist_t distance){
      _dist += distance;
   }
#ifdef FLAGS
   inline void set_state(states s){
      _flags = s;
   }
   inline bool is_approx(){
      return (_flags == exact);
   }
#endif
   /** Other (operators etc) ***/
   bool operator<(Edge const& rhs) const; //compare by distance, than i,j

   
};

/**********************************************************************************
 *********************   WEIGHTED EDGE (i.e. with count) **************************
 **********************************************************************************/
// this class also adds a '+' operator for edge arithmetic  
class WeightedEdge : public Edge {
   friend std::istream& operator>>(std::istream& is,WeightedEdge   & e);      //assign to existing edge
   friend std::istream& operator>>(std::istream& is,WeightedEdge * & e);     //e assigned new Edge(..) or  NULL 
protected:
   count_t _count;
public:
   WeightedEdge() : Edge() , _count(INVALID_COUNT) { ; }
   WeightedEdge(Edge const& e,count_t count = 1) : Edge(e) , _count(count) {
      _dist *= static_cast<dist_t>(_count);
      ;
   }
   WeightedEdge(WeightedEdge const& w) : Edge(w), _count(w.count()) { ; }
   
   
   inline count_t count() const{
      return _count;
   }
   ostream & print_edge(ostream& os) const{
      //return Edge::print_edge(os) << '\t' << count();
      Edge e(i(),j(),dist()/static_cast<dist_t>(count()));
      return e.print_edge(os) << '\t' << count();

   }
   WeightedEdge& operator+=(WeightedEdge const& w){
      assert(i() == w.i());
      assert(j() == w.j());
      
      // fix#1: we add the distance weighted since we assume it was normalized in w
      // so we de-normalize it, to prevent normalizing twice in the end
      //   moved to operator>>
      //_dist +=  (w._count * w._dist); 
      _dist += w._dist;
      _count += w._count;
      return *this;
   }
};


/**********************************************************************************
 ************************   EDGE I/O         **************************************
 **********************************************************************************/

/**
 * formatted Edge print to stream
 */



inline std::ostream& operator<<(std::ostream& os,Edge const & e){
   return e.print_edge(os);
}


inline std::ostream& operator<<(std::ostream& os, WeightedEdge const& e){
   return e.print_edge(os);// << '\t' << e.count();
}	


/**
 * formatted Edge* print to stream, NULL pointer proof
 */

inline std::ostream& operator<<(std::ostream& os,const Edge * e){
   if (e == NULL)    return os << "Edge * -> [NULL]";
   else              return os << "Edge * -> [" <<*e<<"]";
}


/**
 * Build an Edge* from istream - allocates edge pointer on heap, or sets it to null
 */
inline std::istream& operator>>(std::istream& is,Edge *& e){
//    index_t i,j; dist_t d;
   
//    if (is>>i>>j>>d){
// #ifdef INPUT_I_LT_J
//       if (j<i) std::swap(i,j);
// #endif
//       e = new Edge(i,j,d);
//    } else {
//       e = NULL;
//    }
//    return is;
// }

    e = new Edge();   
    if (is>>e->_i>>e->_j>>e->_dist){
#ifdef INPUT_I_LT_J
      if (e->_j<e->_i) std::swap(e->_i,e->_j);
#endif
   } else {
      delete e;
      e = NULL;
   }
   return is;

}


/**
 * Build an Edge from istream
 */
inline std::istream& operator>>(std::istream& is,Edge & e){   
#ifdef SAFE_READ_EDGE
   is.getline(linebuf,MAX_LINE_LENGTH);
   assert(strlen(linebuf) > 4); //"i j e"
   edge_sstream.str(linebuf);   //inits the string stream
   
   if (!(edge_sstream >> e._i >> e._j >> e._dist))
#else 
   if (!(is >> e._i >> e._j >> e._dist))
#endif //SAFE_READ_EDGE
      {
         e._i = INVALID_INDEX;
         e._j = INVALID_INDEX;
         e._dist = INVALID_DIST;   
      }
#ifdef INPUT_I_LT_J
   if (e._j<e._i) std::swap(e._i,e._j);
#endif

   return is;

}

/**
 * Build an Edge from istream
 */
inline std::istream& operator>>(std::istream& is,WeightedEdge & e){   
   if (!(is >> e._i >> e._j >> e._dist >> e._count)){      
      e._i = INVALID_INDEX;
      e._j = INVALID_INDEX;
      e._dist = INVALID_DIST;
      e._count = INVALID_COUNT;
   } else {
      //fix#1: we de-normalize the dist by multiplying in the count, see explanation in operat
      //   we add the distance weighted since we assume it was normalized in w
      //   so we de-normalize it, to prevent normalizing twice in the end      
      e._dist *= e._count;
   }
#ifdef INPUT_I_LT_J
   if (e._j<e._i) std::swap(e._i,e._j);
#endif
   return is;
}




/**********************************************************************************
 *********************   EDGE COMPARATORS   ***************************************
 **********************************************************************************/


struct LT_Edge_dist {
   bool operator()(const Edge* s1, const Edge* s2) const
   {
      if (s1->dist() == s2->dist()){
	 if (s1->i() == s2->i()){
	    return  (s1->j() < s2->j());
	 }
	 return  (s1->i() < s2->i());
      }	
      return (s1->dist() < s2->dist());
   }
};
struct GT_Edge_dist {
   bool operator()(const Edge* s1, const Edge* s2) const
   {
      if (s1->dist() == s2->dist()){
	 if (s1->i() == s2->i()){
	    return  (s1->j() > s2->j());
	 }
	 return  (s1->i() > s2->i());
      }	
      return (s1->dist() > s2->dist());      
   }
};
// is left < right using i and then j ?
struct IJ_LT_Edge {
   bool operator()(const Edge* s1, const Edge* s2) const {
      return ( //at least 1 comparison, at most 3 comparison
              (s1->i() < s2->i())
              ||
              ((s1->i() == s2->i()) && s1->j() < s2->j())
              );
   }
   bool operator()(const Edge& s1, const Edge& s2) const {
      return ( //at least 1 comparison, at most 3 comparison
              (s1.i() < s2.i())
              ||
              ((s1.i() == s2.i()) && s1.j() < s2.j())
              );
   }
   bool operator()(const WeightedEdge& s1, const WeightedEdge& s2) const {
      return ( //at least 1 comparison, at most 3 comparison
              (s1.i() < s2.i())
              ||
              ((s1.i() == s2.i()) && s1.j() < s2.j())
              );
   }

   
         
};

struct IJ_Hasher : public __gnu_cxx::hash<Edge> {
   //the return hashcode is dependent on I,J but not on the distance (or weight)
   size_t operator()(Edge const& edge) const {
      return edge.i() + edge.j();
   }
   size_t operator()(WeightedEdge const& edge) const {
      return edge.i() + edge.j();
   }
   size_t operator()(Edge * const& edge) const{
      return edge->i() + edge->j();
   }
   
};

struct IJ_equals {
   bool operator()(const Edge* s1, const Edge* s2) const {
      return (s1->i() == s2->i() && s1->j() == s2->j());
   }
   bool operator()(const Edge s1, const Edge s2) const {
      return (s1.i() == s2.i() && s1.j() == s2.j());
   }
};
#include <ext/hash_set>
//template<size_t N>
struct IJ_hasher : public __gnu_cxx::hash<const Edge *> {
   static const size_t N = 10 * 1000 * 1000;
   size_t operator()(const Edge * e) const {
      return N*(size_t)e->i() + (size_t)e->j();
   }	
};

static struct LT_Edge_dist s;
inline
bool
Edge::operator<(Edge const& rhs) const {
   //std::clog << "\n\t\tcomparing:"<< *this << "<" << rhs << "==" << s(this,&rhs)<<std::endl;
   return s(this,&rhs);
};





// //static friend:
// inline bool re_approx(Edge * e,vector<Edge *> tree,Heap heap){
//    assert(e.is_aprrox());
//    assert(_flags != both_approx); // otherwise ? ? ?
//    // ==> either i or j (not both) is the approximation
//    index_t aprrox_node,exact_node;
//    if (_flags == i_approx) {
//       approx_node = i; exact_node = j;
//    }  else {
//       assert(e._flags = j_approx);
//       approx_node = j; exact_node = i;
//    }
//    index_t I,J;
//    I = tree[approx_node - N]->i();
//    J = tree[approx_node - N]->j();
//    Edge * tmp;
//    dist_t d_I_exact_node,d_J_exact_node;
//    bool is_approx = false;
//    //Check if we have the missing edge by now, if not take max from the heap
//    // and hope that it grew since last approximation
//    if ((tmp = heap(I,exact_node)) != NULL){
//       d_I_exact_node = tmp->dist();
//    } else {
//       assert(alive[I] && alive[exact_node]); //? ? ?
//       is_approx = true;
//       d_I_exact_node = heap.max_dist();
//    } 

//    if ((tmp = heap(J,exact_node)) != NULL){
//       d_J_exact_node = tmp->dist();
//    } else {
//       assert(alive[I] && alive[exact_node]); //? ? ?
//       is_approx = true;
//       d_J_exact_node = heap.max_dist();
//    }
//    if (!is_approx)
//       e.set_state(exact);
   
// }
   
  
   
//    if ((tmp = heap(J,j)) == NULL){
// 	    assert(alive[J] && alive[j]);
// 	    is_approx = true;
// 	    d_Jj = default_weight;
// 	 } else {
// 	    d_Jj = tmp->dist();
// 	 }
#endif //_EDGE_HPP
