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
 * file: ClusteringUtil.cpp
 *
 * Yaniv Loewenstein, 2007
 *
 * See ClusteringUtil.hpp
 * 
 ***********************************************************************/

#include "ClusteringUtil.hpp"
#include <cassert>
#include <iostream>


void 
ClusteringUtil::nodes2parents(index_t * parents,EdgePtrReader& er,std::ostream& os){
#ifdef PATH_COMPRESSION
   //we don't want to modify parents directly, and we're not using this function ATM anyways, so
   //   won't waste time on modifying it.
   cerr << "warning: PATH_COMPRESSION macro not implmeneted for nodes2parents(..) in single output stream version." << endl;
#endif
   Edge * e = NULL;
   index_t i,j;
   while(er.has_next()){	 
      e = er.next();
      assert(e!=NULL);
      std::clog << ">>\t from:" << *e ;
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

      
      if (i!=j){
	 Edge e_out(i,j,e->dist());
	 os << e_out << std::endl;
	 //std::clog << "\tto: "<<e_out<<std::endl;
      } //else { std::clog << "\tself" << std::endl;}
      delete e;
   }
}


void
run_or_die(std::string const & command){
   std::clog << "will run system(`" << command <<"`)" << std::endl;
   if (system(command.c_str()) != 0){
      die(("error running:`"+command+"`"));
   }
}
void
collate_edges(WeightedEdgeReader& wer/*,T const & sizes*/ ,WeightedEdgesSet & weSet) {
//a map from edge holding comulative distance to count

   assert(weSet.empty());
   WeightedEdgesSet::value_type w,w_prev;   //WeightedEdge w;
   WeightedEdgesSet::iterator   itr = weSet.end();

   while(wer.has_next()){
      w = wer.next();
      if (VERBOSE) clog << w << " - ";
      
      if ((itr = weSet.find(w)) == weSet.end()){      

	 weSet.insert(w);
         if (VERBOSE) clog << "first time" <<endl;
      }

      else {
         //can't modify the key, so we remove it and put it back after
         // modifying in the same position (efficient)
         // *itr += w;  //can't for set, even though I,J stay the same

         w_prev = *itr;
         if (VERBOSE) clog << "was inserted previously:" << w_prev << endl;
         assert(w.i() == w_prev.i());
         assert(w.j() == w_prev.j());
         /* === for set: ===*/
         weSet.erase(itr++); 
         if (itr != weSet.begin() ){
            weSet.insert(--itr,w_prev += w); //insert in the correct place
         } else {
            weSet.insert(itr,w_prev += w); //can't decrement iterator past the beginning 
         }
         /* === for hash_set === */
         //weSet.erase(itr);
         //weSet.insert(w_prev += w); //insert in the correct place

      }
      /*      */
   }
      
  
}
