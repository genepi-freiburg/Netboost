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
 * Major Revisions
 *
 * 26/08/08 - lonshy added the ALLOW_INEXACT_MERGES macro and related code.
 *
 */

using namespace std;

#include "HierarchicalClustering_with_unknown_edges.hpp"
#include <math.h>

RCPtr<tHierarchicalClustering::fAverager> tHierarchicalClustering::fAverager::generate(string const & type) {
  if (type == "arithmetic")
    return new fArithmeticAverage;
  if (type == "geometric")
    return new fGeometricAverage;
  if (type == "harmonic")
    return new fHarmonicAverager;
  if (type == "minimum")
    return new fMinimum;
  if (type == "maximum")
    return new fMaximum;
  throw (tException() << "average can only be one of: {arithmetic, geometric, minimum, maximum}");
}

tHierarchicalClustering::fAverager::~fAverager() {
}

tHierarchicalClustering::tDistance tHierarchicalClustering::fArithmeticAverage::operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2) {
  return (v1 * n1 + v2 * n2) / (n1 + n2);
}

tHierarchicalClustering::tDistance tHierarchicalClustering::fGeometricAverage::operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2) {
  return pow((pow(v1,n1) * pow(v2,n2)), 1.0 / (n1 + n2));
}

tHierarchicalClustering::tDistance tHierarchicalClustering::fMinimum::operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2) {
  return (v1<v2?v1:v2);
}

tHierarchicalClustering::tDistance tHierarchicalClustering::fMaximum::operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2) {
  return (v1>v2?v1:v2);
}

tHierarchicalClustering::tDistance tHierarchicalClustering::fHarmonicAverager::operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2) {
  return (n1 + n2) / (n1 / v1 + n2 / v2);
}


/**
 * C-tor
 */

tHierarchicalClustering::tHierarchicalClustering(tHierarchicalClustering::tEdges& edges, tClusterId singletonIdUB, tDistance lb_unknown, tDistance ub_unknown, RCPtr<tHierarchicalClustering::fAverager> averager, bool is_allow_inexact_merges ) :
   lb_unknown(lb_unknown), // lambda
   ub_unknown(ub_unknown), // psi
   nextClusterId(singletonIdUB), //k
  edges(edges),
  ubEdgeHeap(edges.begin(), edges.size(), CompareEdgesByUpperBound(isLeftGreater)),
  lbEdgeHeap(edges.begin(), edges.size(), CompareEdgesByLowerBound(isLeftGreater)),
  clusters(singletonIdUB * 2),
   averager(averager),
   IS_ALLOW_INEXACT_MERGES(is_allow_inexact_merges) {
  initNeighborLists();
}



/*
 * Assumptions:
 *  existingClusterId is one of the clusters of the edge a priori
 *  existingClusterId and newClusterId are both ids of valid clusters
 *  newClusterId is higher then the id of any other cluster alive
 *  existingClusterId is higher then the id of any current neighbor of newCluster
 */
void tHierarchicalClustering::resetEdge(tEdgeId edgeId, tClusterId existingClusterId, tClusterId newClusterId, tDistance distance_lb, tDistance distance_ub) {
  if (isnan(distance_lb) || isnan(distance_ub))
    throw (tException() << "distance nan " << edgeId);
  tEdge& edge = edges[edgeId];
//  cerr << "resetEdge. existing: " << existingClusterId << " new: " << newClusterId << " before: " << edge << flush;
  if (existingClusterId == edge.lowCluster)
    edge.highCluster = newClusterId;
  else {
    assert(existingClusterId == edge.highCluster);
    edge.lowCluster = edge.highCluster;
    edge.highCluster = newClusterId;
  }
  tCluster& existingCluster = clusters[existingClusterId];
  if (existingCluster.neighborsEnd == existingCluster.neighborsRegionEnd)
    clearInvalidNeighbors(existingCluster);
  assert(existingCluster.neighborsEnd < existingCluster.neighborsRegionEnd);
  neighbors[existingCluster.neighborsEnd++] = make_pair(newClusterId, edgeId);
  tCluster& newCluster = clusters[newClusterId];
  assert(newCluster.neighborsEnd < newCluster.neighborsRegionEnd);
  neighbors[newCluster.neighborsEnd++] = make_pair(existingClusterId, edgeId);
  edge.distance_lb = distance_lb;
  edge.distance_ub = distance_ub;
  ubEdgeHeap.changed(edgeId);
  lbEdgeHeap.changed(edgeId);
//  cerr << " after: " << edge << endl;
}

void tHierarchicalClustering::clearInvalidNeighbors(tCluster& cluster) {
//   cerr << "clearing invalid neighbors " << endl;
//   cerr << "before:";
//   for(tNeighbors::const_iterator i = neighbors.begin() + cluster.neighborsBegin;
//       i != neighbors.begin() + cluster.neighborsEnd;
//       ++i)
//     cerr << " (" << i->first << ": " << edges[i->second] << ')';
//   cerr << endl;
   //lonshy - STL function
  cluster.neighborsEnd = remove_if(neighbors.begin() + cluster.neighborsBegin,
				   neighbors.begin() + cluster.neighborsEnd,
				   fIsClusterInvalid(clusters)) - neighbors.begin();
//   cerr << "after:";
//   for(tNeighbors::const_iterator i = neighbors.begin() + cluster.neighborsBegin;
//       i != neighbors.begin() + cluster.neighborsEnd;
//       ++i)
//     cerr << " (" << i->first << ": " << edges[i->second] << ')';
//   cerr << endl;
}
//lonshy - somehow remove 'invalid' neighbors
void tHierarchicalClustering::compactNeighborList() {
  cerr << "compacting neighbor list. size before: " << clusters[nextClusterId - 1].neighborsRegionEnd << flush;
  unsigned int pos = 0;
  for(tClusterId i = 0; i != nextClusterId; ++i) {
    tCluster& c = clusters[i];
    c.neighborsEnd
      = remove_copy_if(neighbors.begin() + c.neighborsBegin,
		       neighbors.begin() + c.neighborsEnd,
		       neighbors.begin() + pos,
		       fIsClusterInvalid(clusters)) - neighbors.begin();
    c.neighborsBegin = pos;
    pos += static_cast<unsigned int>((c.neighborsEnd - c.neighborsBegin) * 1.1) + 2;
    assert(pos <= c.neighborsRegionEnd);
    c.neighborsRegionEnd = pos;
  }
  cerr << " size after: " << pos << endl;
}

namespace {
  struct fCompareNeighbor {
    bool operator()(pair<tHierarchicalClustering::tClusterId, tHierarchicalClustering::tEdgeId> const & l,
		    pair<tHierarchicalClustering::tClusterId, tHierarchicalClustering::tEdgeId> const & r) const {
      return (l.first < r.first);
    }
  };
}

void tHierarchicalClustering::initNeighborLists() {
  unsigned int neighborListSize = 0;
  {
     //neighbours per cluster:
    vector<unsigned int> neighborCount_tmp(nextClusterId, 0);
    for(tEdges::const_iterator i = edges.begin();
	i != edges.end();
	++i) {
      ++neighborCount_tmp[i->lowCluster];
      ++neighborCount_tmp[i->highCluster];
    } //end for

    for(tClusterId i = 0; i != nextClusterId; ++i) {
      tCluster& c = clusters[i];
      c.neighborsBegin = neighborListSize;
      c.neighborsEnd = neighborListSize;
      neighborListSize += static_cast<unsigned int>(neighborCount_tmp[i] * 1.1) + 2;
      c.neighborsRegionEnd = neighborListSize;
    } //end for
  }
  neighbors.resize(static_cast<unsigned int>(neighborListSize * 1.1));
  cerr << "allocated space for neighbor list" << endl;

  for(tEdgeId i = 0; i != edges.size(); ++i) {
    neighbors[clusters[edges[i].lowCluster].neighborsEnd++] = make_pair(edges[i].highCluster, i);
    assert(clusters[edges[i].lowCluster].neighborsEnd <= clusters[edges[i].lowCluster].neighborsRegionEnd);
    neighbors[clusters[edges[i].highCluster].neighborsEnd++] = make_pair(edges[i].lowCluster, i);
    assert(clusters[edges[i].highCluster].neighborsEnd <= clusters[edges[i].highCluster].neighborsRegionEnd);
  }
  cerr << "initial neighbor lists loaded" << endl;

  for(tClusters::iterator i = clusters.begin();
      i != clusters.end();
      ++i) {
    i->size = 1;
    sort(neighbors.begin() + i->neighborsBegin, neighbors.begin() + i->neighborsEnd, fCompareNeighbor());
  }
  cerr << "initial neighbor lists sorted" << endl;
}

/**
 * nextMerge() - the actual clustering main routine implementation
 */

   /* lonshy check that the merge is a valid edge whose merged_edge.ub <= min(lb) (including self)
      stop if
      (1) mergingEdge.ub is not strictly greater than the smallest possible on disk (lambda = lb_unknown)
      OR
      (2) mergingEdge.ub is     stritcly greater than some lb (self included)

      Condition (2) is relaxed when the cluster height (dendrogram) is not required.

      when inexact merges are allowed, the minimal edge might be an interval, i.e. lb < ub, and we allow it.
      we need to make sure that the interval does not overlap any *other* interval for which we look at second lb
      using the second_top() method of the heap. This method is O(logN) instead of the O(1) top() method, so we only use
      it when the O(1) test is not sufficient.
   */


tHierarchicalClustering::tMergeEvent tHierarchicalClustering::nextMerge() {
   tEdgeId mergingEdgeId = ubEdgeHeap.top();    
   tEdge&  mergingEdge   = edges[mergingEdgeId]; //edge of minimal lb
   if (! isLeftGreater(lb_unknown, mergingEdge.distance_ub)) //condition (1)
#ifdef ALLOW_LAMBDA_MERGES
   if (isLeftGreater(mergingEdge.distance_ub,lb_unknown)) 
      return false;
   //otherwise, ===> (mergingEdge.distance_ub == lb_unknown) ===> continue!
#else
   return false; //stop since lambda > current_edge_ub does not hold 
#endif
      
   if (isLeftGreater(mergingEdge.distance_ub, edges[lbEdgeHeap.top()].distance_lb)){ // condition (2)
      if (IS_ALLOW_INEXACT_MERGES) {
	 if (isLeftGreater(mergingEdge.distance_ub, edges[lbEdgeHeap.second_top()].distance_lb))
	    return false;
	 //else continue - this interval does not overlap the next one
      } else {
	 return false; //the minimal edge is inexact since ub > lb ==> can't continue merging
      }
   }
   tMergeEvent event(mergingEdge.lowCluster, mergingEdge.highCluster, mergingEdge.distance_ub, nextClusterId);
//  cerr << "merging: " << mergingEdge << endl;

//   {
//     unsigned int pos = 0;
//     for (tClusterId i = 0; i != nextClusterId; ++i) {
//       tCluster& c = clusters[i];
//       if (! (c.neighborsBegin == pos) &&
// 	  (c.neighborsBegin <= c.neighborsEnd) &&
// 	(c.neighborsEnd <= c.neighborsRegionEnd))
// 	cerr << "range problem: " << i << ' ' << pos << ' ' << c.neighborsBegin << ' ' << c.neighborsEnd << ' ' << c.neighborsRegionEnd << endl;
//       pos = c.neighborsRegionEnd;
//     }
//     if (!(pos <= neighbors.size()))
//       cerr << "range problem: " << pos << ' ' << neighbors.size() << endl;
//   }


//lonshy added for more informative error messages:
  if (! (
	 (IS_ALLOW_INEXACT_MERGES || mergingEdge.exact()) &&

	 mergingEdge.valid() &&  (mergingEdge.lowCluster < mergingEdge.highCluster) && 
         (clusters[mergingEdge.lowCluster].valid() && clusters[mergingEdge.highCluster].valid()))){
     
     
     cerr << "an error occured in " << __FILE__ << ":" << __LINE__ <<  ", when merging edge:" << mergingEdge; //endl contained
     
     
        if (!mergingEdge.valid()) 
           cerr << "Hint: the input seems to contain an illegal \"0\" cluster-ID"  <<endl;
  }     //end of lonshy add

 
  assert(mergingEdge.valid());
  assert(IS_ALLOW_INEXACT_MERGES || mergingEdge.exact());
  assert(mergingEdge.lowCluster < mergingEdge.highCluster);
  assert(clusters[mergingEdge.lowCluster].valid());
  assert(clusters[mergingEdge.highCluster].valid());
  
  tCluster& lowCluster = clusters[mergingEdge.lowCluster];
  tCluster& highCluster = clusters[mergingEdge.highCluster];
  tCluster& newCluster = clusters[nextClusterId];

  newCluster.neighborsBegin = clusters[nextClusterId - 1].neighborsRegionEnd;
  newCluster.neighborsEnd = newCluster.neighborsBegin;
  newCluster.neighborsRegionEnd = newCluster.neighborsBegin + static_cast<unsigned int>(((lowCluster.neighborsEnd - lowCluster.neighborsBegin) + (highCluster.neighborsEnd - highCluster.neighborsBegin)) * 1.1);
  if (newCluster.neighborsRegionEnd > neighbors.size()) {
    compactNeighborList();
    newCluster.neighborsBegin = clusters[nextClusterId - 1].neighborsRegionEnd;
    newCluster.neighborsEnd = newCluster.neighborsBegin;
    newCluster.neighborsRegionEnd = newCluster.neighborsBegin + static_cast<unsigned int>((lowCluster.neighborsEnd - lowCluster.neighborsBegin + highCluster.neighborsEnd - highCluster.neighborsBegin) * 1.1);
    if (newCluster.neighborsRegionEnd > neighbors.size()) {
      cerr << "resizing neighbors list ... " << flush;
      neighbors.resize(newCluster.neighborsRegionEnd);
      cerr << "done" << endl;
    }
  }
  //lonshy: lowClusterIterator, lowClusterIteratorEnd, highClusterIterator, highClusterIteratorEnd
  // used to iterate neighbours
  vector<pair<tClusterId, tEdgeId> >::const_iterator li = neighbors.begin() + lowCluster.neighborsBegin; 
  vector<pair<tClusterId, tEdgeId> >::const_iterator lie = neighbors.begin() + lowCluster.neighborsEnd;
  vector<pair<tClusterId, tEdgeId> >::const_iterator hi = neighbors.begin() + highCluster.neighborsBegin;
  vector<pair<tClusterId, tEdgeId> >::const_iterator hie = neighbors.begin() + highCluster.neighborsEnd;
//   cerr << "lowCluster neighbor range: " << lowCluster.neighborsBegin << ' ' << lowCluster.neighborsEnd << ' ' << lowCluster.neighborsRegionEnd << endl;
//   cerr << "lowCluster neighbors:";
//   for(tNeighbors::const_iterator i = neighbors.begin() + lowCluster.neighborsBegin;
//       i != neighbors.begin() + lowCluster.neighborsEnd;
//       ++i)
//     cerr << " (" << i->first << ": " << edges[i->second] << ')';
//   cerr << endl;
//   cerr << "highCluster neighbor range: " << highCluster.neighborsBegin << ' ' << highCluster.neighborsEnd << ' ' << highCluster.neighborsRegionEnd << endl;
//   cerr << "highCluster neighbors:";
//   for(tNeighbors::const_iterator i = neighbors.begin() + highCluster.neighborsBegin;
//       i != neighbors.begin() + highCluster.neighborsEnd;
//       ++i)
//     cerr << " (" << i->first << ": " << edges[i->second] << ')';
//   cerr << endl;

  //lonshy this looks like some sort of linear search for the two edges li lj that are averaged into one edge lk between new cluster k
  // and i U j neighbours
  while((li != lie) && (hi != hie)) {
     while(li->first < hi->first) { //inner while #1
	if ((li->first != mergingEdge.highCluster) && clusters[li->first].valid())
	   //lonshy: averaging known with missing
	   resetEdge(li->second,    li->first,     nextClusterId,
		     (*averager)(lowCluster.size, edges[li->second].distance_lb,highCluster.size, lb_unknown),
		     (*averager)(lowCluster.size, edges[li->second].distance_ub,highCluster.size, ub_unknown));
	++li;
	if (li == lie)
	   break;
     } //end of inner while #1, li changed
     if (li == lie)  break;
     while(hi->first < li->first) {//inner while #2
	if ((hi->first != mergingEdge.lowCluster) && clusters[hi->first].valid())
	   //lonshy: averageing missing with known
	   resetEdge(hi->second,     hi->first,    nextClusterId,
		     (*averager)(lowCluster.size, lb_unknown, highCluster.size, edges[hi->second].distance_lb),
		     (*averager)(lowCluster.size, ub_unknown, highCluster.size, edges[hi->second].distance_ub));
	++hi;
	if (hi == hie)  break;
     } // end of inner while #2, hi changed
     if (hi == hie)     break;
     if (hi->first == li->first) {
	//lonshy averaging two 'known' edges - into one, the 2nd is removed
	if (clusters[hi->first].valid()) {
	   resetEdge(hi->second,     hi->first,    nextClusterId,
		     (*averager)(lowCluster.size, edges[li->second].distance_lb,  highCluster.size, edges[hi->second].distance_lb),
		     (*averager)(lowCluster.size, edges[li->second].distance_ub,  highCluster.size, edges[hi->second].distance_ub));
	   removeEdge(li->second);
	}
	++li;
	++hi;
     }
  } //end of outer while
  for(; li != lie; ++li)
     if ((li->first != mergingEdge.highCluster) && clusters[li->first].valid())
	resetEdge(li->second,	  li->first,	  nextClusterId,
		  (*averager)(lowCluster.size, edges[li->second].distance_lb,   highCluster.size, lb_unknown),
		  (*averager)(lowCluster.size, edges[li->second].distance_ub,   highCluster.size, ub_unknown));
  for(; hi != hie; ++hi)
     if ((hi->first != mergingEdge.lowCluster) && clusters[hi->first].valid())
	resetEdge(hi->second,	  hi->first, 	nextClusterId,
		(*averager)(lowCluster.size, lb_unknown,  highCluster.size, edges[hi->second].distance_lb),
		(*averager)(lowCluster.size, ub_unknown,  highCluster.size, edges[hi->second].distance_ub));
  
  newCluster.size = lowCluster.size + highCluster.size;
  lowCluster.invalidate();
  highCluster.invalidate();
  
  ++nextClusterId;
  removeEdge(mergingEdgeId);
//   if (newCluster.size >= maxClusterSize) {
//     for (li = neighbors.begin() + newCluster.neighborsBegin,
// 	   lie = neighbors.begin() + newCluster.neighborsEnd;
// 	 li != lie;
// 	 ++li)
//       removeEdge(li->second);
//     newCluster.invalidate();
//   }
  return event;
}

