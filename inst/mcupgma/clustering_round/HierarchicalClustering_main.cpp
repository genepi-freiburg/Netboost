/*****************************************************************************
 * MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA    *
 *             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.   *
 *                                                                           *
 *                                                                           *
 * Copyright (C), 2007  Elon Portugaly, Yaniv Loewenstein                    *
 *                School of Computer Science And Engineering                 *
 *                Hebrew University of Jerusalem                             *
 *                                                                           *
 *      All Rights Reserved                                                  *
 *                                                                           *
 *      This source code is distributed under the terms of the               *
 *      GNU General Public License. See the file LICENSE                     *
 *      for details.                                                         *
 *                                                                           *
 *****************************************************************************/

/**
 * LOG: 
 * - lonshy (13/08/2007) - added to event a merge counter
 * - Documentation added by lonshy 2/8/2007
 * - More documentation added by lonshy 5/8/2008
 * - added macro CHECK_INPUT_VALIDITY
 *        - checks for 0 clusterID
 *        - check for self edges
 *        - warn on edges >= max_distance argument (or just > when ALLOW_PSI_EDGES)
 * - 03/09/08 remove VERBOSE_HEAP macro
 * - 03/09/08 - lonshy added ALLOW_PSI_EDGES and changed main() and removeEdge() accordingly
 * - 03/09/08 - changed verbose output when the heap is partially loaded and the max loaded edge can't be
 *              retrieved from the heap directly. moved scope of "loaded edges..." printout.
 * - 07/08/08 - added option to determine behaviour for minimal inexact intervals at runtime.
 */


/* ---------------------------------------------------------------------
 * - 03/09/08 allow merging edges == psi - log of changes
 *
 * problem: since input for require lb < heap.lb edges the size of the maximum may not enter the heap, and will not be clustered
 * solution:   It's easy to allow insertion to the heap of these edges, by changing the if from '<'
 *   to '<=' in the heap insertion for loop (in the main). However, to allow their merging we need to allow to
 *   merge edges whose ub <= lambda (rather than the present strictly <). However, these valid edges (==psi) might pop
 *   out of the heap in the clustering process after invalid edges (which are also equal psi, due to the heap initialization).
 *   To solve the latter, we can initialize the heap elements to be '2.0 * max_distance_arg' instead of 'max_distance_arg'.
 *   But, now we can't stop the clustering in time, and we'll continue to merge into invalid edges. So we change
 *   removeEdge(...) to set the distances to 2.0 * ub_unknown as well.
 *
 * Summary of changes: 
 * 1. replaced < with <= in if in heap insertion for loop (main)
 * 2. replaced in nextMerge(...):
 *    if (! compareDistances(lb_unknown, mergingEdge.distance_ub))  return false;
 *     with
 *    if (! compareDistances(lb_unknown, mergingEdge.distance_ub))
 *         if (compareDistances(mergingEdge.distance_ub,lb_unknown))
 *                 return false;
 * 3. in main replace 
 *     vector<tHierarchicalClustering::tEdge> edges(args_info.number_of_input_edges_arg, args_info.max_distance_arg);
 *      with
 *     vector<tHierarchicalClustering::tEdge> edges(args_info.number_of_input_edges_arg, 2.0 * args_info.max_distance_arg);
 * 4. in removeEdge(...) replaced 
 *    edges[edgeId].distance_lb = ub_unknown;
 *    edges[edgeId].distance_ub = ub_unknown;
 *     with:
 *    edges[edgeId].distance_lb = 2.0 * ub_unknown;
 *    edges[edgeId].distance_ub = 2.0 * ub_unknown;
 *
 *    changes 1,2 are not required, once invalid edges are equal 2.0 * psi by changes 3 and 4, so in the bottom line
 *    the macro ALLOW_PSI_EDGES applies changes 3 and 4.
 *
 *    Update: we also override the definition of max_loaded_edge with psi when the heap was not used entirely by 
 *	      input edges, otherwise it is mistakenly believed to be 2*psi and the clustering is stopped 
 *            prematurely.
 * ------------------------------------------------------------------------
 */
 




using namespace std;


/** Some optional lonshy macros **/

//some debug printing for lonshy
//#define VERBOSE 


#define COUNT_MERGES
#define CHECK_INPUT_VALIDITY //lonshy: add some simple checks for input validity
//depreceted : #define WARN_MAX_DIST_ARG //warn on edges >= max dist arg when CHECK_INPUT_VALIDITY is set.

/** end of lonshy macros **/


#include "HierarchicalClustering_main.cmdline.h"
#include "HierarchicalClustering_with_unknown_edges.hpp"
#include <myutils/StreamFromFileNameGenerator.hpp>
#include <myutils/Exception.hpp>
#include <iterator>





/*********************************************
 ** struct fMergeEventLogger
 *********************************************/

struct fMergeEventLogger {
  fMergeEventLogger(ostream& os, tHierarchicalClustering::tDistance& max_clustered_distance
#ifdef COUNT_MERGES
		    , unsigned int& _count) : count(_count),
#else
					      ) :
#endif
   os(os),
   max_clustered_distance(max_clustered_distance) {
   }
void operator()(tHierarchicalClustering::tMergeEvent const & e) const {
   os << e;
   max_clustered_distance = e.distance;
#ifdef COUNT_MERGES
   ++count;
#endif
}

private:
#ifdef COUNT_MERGES
   unsigned int& count;
#endif
   ostream& os;
   tHierarchicalClustering::tDistance& max_clustered_distance;

};




/*********************************************
 ** main()   
 *********************************************/

int main(int argc, char** argv) {
   try {
      gengetopt_args_info args_info;
      if (cmdline_parser (argc, argv, &args_info) != 0)
	 throw tException();
      //lonshy bug fix, for user entering '-number-of-input-edges= 10' instead of '--number-of-input-edges=10'
      if (args_info.number_of_input_edges_arg <= 0){
	 cerr << "number of edges must be > 0, you requested:"
	      <<args_info.number_of_input_edges_arg
	      << " (hint: make sure that there is no whitespace after '=' in command line) "
	      << endl;
	 exit(1);
      }	   
      if (args_info.allow_non_dendrogram_flag)
	 cerr << "*** Warning: allowing inexact merges *** - tree heights are meaningless." << endl <<
	    "\t(at present, the average of the lower and upper bounds are output)." << endl;
      else 
	 cerr << "Allowing only exact merges, cluster heights are given in output tree (dendrogram)" << endl;

#ifdef CHECK_INPUT_VALIDITY   
      tHierarchicalClustering::tDistance max_distance_arg = static_cast<tHierarchicalClustering::tDistance>(args_info.max_distance_arg);
#endif      
    
      RCPtr<istream> in;

      /* A vector of N edges, N from command line, initialized as max_distance */
#ifdef ALLOW_PSI_EDGES 
      vector<tHierarchicalClustering::tEdge> edges(args_info.number_of_input_edges_arg, 2.0 * args_info.max_distance_arg);
#else
      vector<tHierarchicalClustering::tEdge> edges(args_info.number_of_input_edges_arg, args_info.max_distance_arg);
#endif
      
      tHierarchicalClustering::tDistance max_loaded_distance;
    {

       /** a reverse heap -  max at top, not min (used to find the minimal M edges) [lonshy]*/
       tHeap<
       vector<tHierarchicalClustering::tEdge>::const_iterator  ,
	    tHierarchicalClustering::fCompareEdgesByLowerBound<less<tHierarchicalClustering::tDistance> >
          > reverse_heap(edges.begin(), edges.size());
       in = InputStream(args_info.input_cluster_edges_file_arg);
       
       tHierarchicalClustering::tEdge e;
       //lonshy added for monitoring:
       unsigned int count_insertions = 0; //how many times did we actually insert to the heap ( != heap.size(), if input isn't sorted)
      
       /** read edges loop: *******/
       for(*in >> e; 	  *in;	  *in >> e) {



#ifdef CHECK_INPUT_VALIDITY
          if (e.lowCluster == e.highCluster) {
             cerr << "Input contains an illegal self edge:" << e << endl;
             exit(1);
          }
          else if (e.lowCluster == 0) {
             cerr << "Input contains an illegal 0 clusterID:" << e << endl;
             exit(1);
          }
#ifndef ALLOW_PSI_EDGES 
	  else if (e.distance_ub >= max_distance_arg){
	     cerr << "Warning: Input contains an edge >= max_distance argument (psi) which will not be clustered:"<< e << endl;
	     continue; //make sure edges between psi and 2*psi won't enter the heap..
	  }
#else
	  else if (e.distance_ub > max_distance_arg){
	     cerr << "Warning: Input contains an edge > max_distance argument (psi) which will not be clustered:"<< e << endl;
	     //it can't enter the heap, since all heap edges <= psi by construction.
	  }
#endif
#endif
#ifdef ALLOW_PSI_EDGES
	  //keep track of the max, in case the heap doesn't become full and we can't find the max of only valid edges.
	  if ((count_insertions < args_info.number_of_input_edges_arg && e.distance_lb > max_loaded_distance) || !count_insertions)
	     max_loaded_distance = e.distance_lb;
#endif
	  tHierarchicalClustering::tEdgeId i = reverse_heap.top();
       //lonshy - replace equal edges by the most recently read edge
       //         this allows insertion of edges == psi
       // if (e.distance_lb <= edges[i].distance_lb) {  //this edge is better than the worst at hand ==> replace
	  if (e.distance_lb < edges[i].distance_lb) {  //this edge is better than the worst at hand ==> replace 
	    count_insertions++;
	    edges[i] = e;
	    reverse_heap.changed(i);
	  }
	  
       } //end for 
       /** finished reading input **/
       
       

       cerr << "finished reading input edges." << endl;
      //cerr << "finished reading input edges, has top()=" << reverse_heap.top() +1<< " edges" << endl; //lonshy added
      cerr << "number of heap insertions:"<<count_insertions<<endl;
#ifdef ALLOW_PSI_EDGES //the maximum heap element may be erroniously 2*psi, if the heap isn't full
	// In which case we have put 2*psi edges in the heap, so we need to override the original setting
	// that checks the heap, and provide the value of psi specifically to the clustering(..) method - not the 
	// actual maximum loaded edge which might be < psi

      if (count_insertions >= args_info.number_of_input_edges_arg)
	 max_loaded_distance = edges[reverse_heap.top()].distance_lb;
      else 
	 max_loaded_distance = static_cast<tHierarchicalClustering::tDistance>(args_info.max_distance_arg);
#else //the maximum heap element will be psi in this case
      max_loaded_distance = edges[reverse_heap.top()].distance_lb;
#endif
      cerr << "loaded edges. maximum loaded distance: " << max_loaded_distance << endl;
    }
    
    
    tHierarchicalClustering clustering(edges,
				       static_cast<tHierarchicalClustering::tClusterId>(args_info.max_cluster_index_arg + 1),
				       max_loaded_distance,  //lambda
				       static_cast<tHierarchicalClustering::tDistance>(args_info.max_distance_arg), //psi
				       tHierarchicalClustering::fAverager::generate(args_info.average_type_arg),
				       args_info.allow_non_dendrogram_flag
				       );

    cerr << "arranged all edges" << endl; //neighbor lists constructed
    
    if (args_info.input_cluster_sizes_file_name_given) {
       in = InputStream(args_info.input_cluster_sizes_file_name_arg);
       istream_iterator<tHierarchicalClustering::tIdAndSize> b(*in);
       istream_iterator<tHierarchicalClustering::tIdAndSize> e;
       clustering.load_cluster_sizes(b,e);
    }
    cerr << "loaded cluster sizes" << endl;
    
    RCPtr<ostream> out = OutputStream(args_info.output_merges_file_name_arg);
    tHierarchicalClustering::tDistance max_clustered_distance;

    
    /**  actual clustering --> call cluster() */

#ifdef COUNT_MERGES
    unsigned int merge_count = 0;
    clustering.cluster(fMergeEventLogger(*out, max_clustered_distance,merge_count));
    cerr << "merged " << merge_count << " clusters" << endl;
    if (merge_count) {
       cerr << "done clustering up to distance " << max_clustered_distance << endl;
    } else {
       cerr << "no clustering done."  << endl;
    }
#else
    clustering.cluster(fMergeEventLogger(*out, max_clustered_distance));
    cerr << "done clustering up to distance " << max_clustered_distance << endl;
#endif

    /**  termination - print out cluster sizes if requested */
    
    if (args_info.output_cluster_sizes_file_name_given) {
      out = OutputStream(args_info.output_cluster_sizes_file_name_arg);
      ostream_iterator<tHierarchicalClustering::tIdAndSize> o(*out);
      clustering.save_cluster_sizes(o);
    }
    
//     if (args_info.output_edges_file_name_given) {
//       out = OutputStream(args_info.output_edges_file_name_arg);
//       ostream_iterator<tHierarchicalClustering::tEdge> o(*out);
//       clustering.save_edges(o);
//     }
    out = 0;
      
  } catch (tException const & e) {
    cerr << e << endl;
    return 11;
  } catch (std::bad_alloc & e) {
    cerr << "caught std::bad_alloc" << endl;
    cerr << e.what() << endl;
    cerr << "Seems like we're out of memory, try to use a smaller value for the number of held edges (the memory constraint)." <<endl;
    return 202;

  } catch (std::exception const & e) {
    cerr << "caught std::exception" << endl;
    cerr << e.what() << endl;
    return 12;
  } catch (...) {
    cerr << "exception thrown" << endl;
    return 13;
  }
  return 0;

}


