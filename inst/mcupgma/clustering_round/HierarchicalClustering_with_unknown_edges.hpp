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

#ifndef __HierarchicalClustering_with_unknown_edges_hpp
#define __HierarchicalClustering_with_unknown_edges_hpp

/**
 * Major Revisions
 *
 * 26/08/08 - lonshy added the ALLOW_INEXACT_MERGES macro and related code.
 * 03/09/08 - lonshy added ALLOW_PSI_EDGES macro and changed removeEdge(...) and main(...)
 * 07/09/08 - ALLOW_INEXACT_MERGES macro replaced with class member IS_ALLOW_INEXACT_MERGES which allows
 *            runtime determination of behaviour for inexact minimal intervals.
 * 14/08/08 - ALLOW_LAMBDA_MERGES - added macro that allows clustering to proceed when all loaded edges are ties
 *                                  by allowing merger of edges whos ub equals the lb on missing edges
 */


#define ALLOW_PSI_EDGES
#define ALLOW_LAMBDA_MERGES //allow merging of edges  whos ub equals the lb on missing edges

/**
 * All comments added by lonshy to the best of his understanding
 */

#include <vector>
#include "Heap.hpp"
#include <algorithm>
#include <myutils/Exception.hpp>
#include <myutils/Writable.hpp>
#include <myutils/Readable.hpp>

class tHierarchicalClustering {
public:
   typedef unsigned int tClusterId;
   typedef unsigned int tEdgeId;
   typedef double tDistance;

    /*********************************************************
     * class tHierarchicalClustering::fAverager
     * An abstraction of a function for averaging distances **
     *********************************************************/
    class fAverager {
    public:
	static RCPtr<fAverager> generate(string const & type);
	virtual tDistance operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2) = 0;
	virtual ~fAverager();
    };
    /** Concrete fAverager types **/
    class fArithmeticAverage : public fAverager {
    public:
	virtual tDistance operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2);
    };
    
    class fGeometricAverage : public fAverager {
    public:
	virtual tDistance operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2);
    };
    
    class fHarmonicAverager : public fAverager {
    public:
	virtual tDistance operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2);
    };
  
    class fMinimum : public fAverager {
    public:
	virtual tDistance operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2);
    };
    
    class fMaximum : public fAverager {
    public:
	virtual tDistance operator()(unsigned int n1, tDistance v1, unsigned int n2, tDistance v2);
    };
   /** End of fAverager **/
   


    
    /************************************************************************
     * struct tHierarchicalClustering::tEdge - data held for each edge
     *
     * lowCluster has a lower id then highCluster.
     * the edge might be invalid - marked by setting lowCluster to 0.
     *
     * edges are sorted by distance.
     */
    
    struct tEdge : public tWritable<tEdge>, public tReadable<tEdge> {
	tEdge() :
	    lowCluster(0),
	    highCluster(0),
	    distance_lb(0.0/0.0),
	    distance_ub(0.0/0.0) {
	}
	
	tEdge(tDistance distance) :
	    lowCluster(0),
	    highCluster(0),
	    distance_lb(distance),
	    distance_ub(distance){
	}
	
	bool valid() const {
	    return (lowCluster != 0);
	}
	
	bool exact() const {
	    return (distance_lb == distance_ub);
	}
	/** Members **/      
	tClusterId lowCluster;
	tClusterId highCluster;
	tDistance distance_lb; //upper bound
	tDistance distance_ub; //lower bound
	
	ostream& Write(ostream& os) const {
	    return os << lowCluster << '\t' << highCluster << '\t' << distance_lb << '\n';
	}
	
	istream& Read(istream& is) {
	    is >> lowCluster >> highCluster >> distance_lb;
	    if (highCluster < lowCluster)
		swap(lowCluster, highCluster);
	    distance_ub = distance_lb;
	    return is;
	}
    }; //end of struct tEdge
      
    /********************************************
     * Comparator of edges by actual distances
     *******************************************/
    typedef greater<tDistance> fIsLeftGreater;
    /** Using fIsLeftGreater as a binary_function object **/
    template<typename fDistanceComparator>
    /** Comparator of edges by lower bound **/
    struct fCompareEdgesByLowerBound : public binary_function<tEdge, tEdge, bool> {
	fCompareEdgesByLowerBound(fDistanceComparator compareDistances = fDistanceComparator()) :
      compareDistances(compareDistances) {
	}
	
	bool operator()(tEdge const & l, tEdge const & r) const {
	    if (compareDistances(l.distance_lb, r.distance_lb))
		return true;
	    else if (compareDistances(r.distance_lb, l.distance_lb))
		return false;
	    else if (l.lowCluster < r.lowCluster)
		return true;
	    else if (l.lowCluster > r.lowCluster)
		return false;
	    else
		return (l.highCluster < r.highCluster);
	}    
    private:
	fDistanceComparator compareDistances;
    }; //end struct fCompareEdgesByLowerBound
    
  template<typename fDistanceComparator>
  fCompareEdgesByLowerBound<fDistanceComparator> CompareEdgesByLowerBound(fDistanceComparator dc) {
    return fCompareEdgesByLowerBound<fDistanceComparator>(dc);
  }
   
  template<typename fDistanceComparator>
  struct fCompareEdgesByUpperBound : public binary_function<tEdge, tEdge, bool> {
    fCompareEdgesByUpperBound(fDistanceComparator compareDistances = fDistanceComparator()) :
      compareDistances(compareDistances) {
    }
    
    bool operator()(tEdge const & l, tEdge const & r) const {
      if (compareDistances(l.distance_ub, r.distance_ub))
	return true;
      else if (compareDistances(r.distance_ub, l.distance_ub))
	return false;
      else if (l.lowCluster < r.lowCluster)
	return true;
      else if (l.lowCluster > r.lowCluster)
	return false;
      else
	return (l.highCluster < r.highCluster);
    }
    
  private:
    fDistanceComparator compareDistances;
  };

  template<typename fDistanceComparator>
  fCompareEdgesByUpperBound<fDistanceComparator> CompareEdgesByUpperBound(fDistanceComparator dc) {
    return fCompareEdgesByUpperBound<fDistanceComparator>(dc);
  }

    /**********************************
     * struct tHierarchicalClustering::tCluster 
     */
   
   struct tCluster {
    tCluster() :
      neighborsBegin(0),
      neighborsEnd(0),
      neighborsRegionEnd(0),
      size(0) {
    }
    
    bool valid() const {
      return (size > 0);
    }

    void invalidate() {
      neighborsEnd = neighborsBegin;
      size = 0;
    }
    
    unsigned int neighborsBegin;
    unsigned int neighborsEnd;
    unsigned int neighborsRegionEnd;
    
    unsigned int size;
  };

   /**********************************
    * struct tHierarchicalClustering::tIdAndSize
    **********************************/
  struct tIdAndSize : public tWritable<tIdAndSize>, public tReadable<tIdAndSize> {
    tIdAndSize() {
    }
    
    tIdAndSize(tClusterId id, unsigned int size):
      id(id),
      size(size) {
    }
    
    ostream& Write(ostream& os) const {
      return os << id << '\t' << size << '\n';
    }
    
    istream& Read(istream& is) {
      return is >> id >> size;
    }
    tClusterId id;
    unsigned int size;
  };

private:
   /*************************************
    * struct tHierarchicalClustering::fIsClusterInvalid
    */
   struct fIsClusterInvalid : unary_function<tClusterId, bool> {
      //c-tor
      fIsClusterInvalid(vector<tCluster> const & clusters) :
         clusters(clusters) {
      }
      //operator() checks that the cluster is valid
      bool operator()(pair<tClusterId, tEdgeId> const & r) const {
         return (! clusters[r.first].valid());
      }
      
   private:
      vector<tCluster> const & clusters;
   }; /* end of struct tHierarchicalClustering::fIsClusterInvalid */


  typedef vector<tEdge> tEdges;
  typedef tHeap<tEdges::const_iterator, fCompareEdgesByUpperBound<fIsLeftGreater> > tUpperBoundHeap;
  typedef tHeap<tEdges::const_iterator, fCompareEdgesByLowerBound<fIsLeftGreater> > tLowerBoundHeap;
  typedef vector<tCluster> tClusters;
  typedef vector<pair<tClusterId, tEdgeId> > tNeighbors;

public:

  /*
   * A merge event contains the ids of the merged clusters, the distance between them, which is the time of the merge event, and the id of the newly created cluster.
   *
   * Evaluates to true if contains data of a merge event, and to false if no merge was performed
   */ 
  struct tMergeEvent : public tReadable<tMergeEvent>, public tWritable<tMergeEvent> {
    tMergeEvent(bool b = false) :
      merged(0) {
      assert(! b);
    }
    
    tMergeEvent(tClusterId c1, tClusterId c2, tDistance distance, tClusterId merged) :
      c1(c1),
      c2(c2),
      distance(distance),
      merged(merged) {
      assert(merged != 0);
    }
    
    operator bool() const {
      return (merged != 0);
    }

    ostream& Write(ostream& os) const {
      return os << /*"lonshy" <<*/ c1 << '\t' << c2 << '\t' << distance << '\t' << merged << endl;
    }

    istream& Read(istream& is) {
      return is >> c1 >> c2 >> distance >> merged;
    }
     /* tMergeEvent members */
     tClusterId c1;
     tClusterId c2;
     tDistance distance;
     tClusterId merged;
     
  };
   /**
    * C-tor
    */
   tHierarchicalClustering(tEdges& edges, tClusterId singletonIdUB, tDistance lb_unknown, tDistance ub_unknown, RCPtr<fAverager> averager, bool is_allow_inexact_merges = false);

  template <typename fListener>
  void cluster(fListener listener);
   //lonshy input iter
  template <typename II>
  II load_cluster_sizes(II cluster_sizes_begin, II cluster_sizes_end);
   //lonshy output iter
  template <typename OI>
  OI save_cluster_sizes(OI cluster_sizes_oi) const;
   
private:
  void resetEdge(tEdgeId edgeId, tClusterId existingClusterId, tClusterId newClusterId, tDistance distance_lb, tDistance distance_ub);
  void removeEdge(tEdgeId edgeId);
  void clearInvalidNeighbors(tCluster& cluster);
  void compactNeighborList();
  void initNeighborLists();
   
  // merges the two nearest clusters and returns the merge event data. returns false if no two clusters with distance shorter then threshold exist.
  tMergeEvent nextMerge();

   /* -----------------------------------
      members of tHierarchicalClustering
      ----------------------------------- */
   const tDistance lb_unknown;
   const tDistance ub_unknown;
//  unsigned int maxClusterSize;
   unsigned int nextClusterId;

  tEdges& edges;
  tUpperBoundHeap ubEdgeHeap;
  tLowerBoundHeap lbEdgeHeap;
  tClusters clusters;
  tNeighbors neighbors;
  
  fIsLeftGreater isLeftGreater;
  RCPtr<fAverager> averager;
   const bool IS_ALLOW_INEXACT_MERGES;
};

/**
 * load_cluster_sizes(...)
 */
template <typename II>
II tHierarchicalClustering::load_cluster_sizes(II cluster_sizes_begin, II cluster_sizes_end) {
   for(;cluster_sizes_begin != cluster_sizes_end;++cluster_sizes_begin)
      //clusters[cluster_sizes_begin->id].size = cluster_sizes_end->size; //lonshy fix
      clusters[cluster_sizes_begin->id].size = cluster_sizes_begin->size;
  return cluster_sizes_begin;
}


/**
 * save_cluster_sizes(...)
 */
template <typename OI>
OI tHierarchicalClustering::save_cluster_sizes(OI cluster_sizes_oi) const {
  for(tClusterId i = 0;
      i != clusters.size();
      ++i)
    if ((clusters[i].valid()) && (clusters[i].size > 1))
      *(cluster_sizes_oi++) = tIdAndSize(i, clusters[i].size);
  return cluster_sizes_oi;
}


/**
 * cluster(...)
 */
template <typename fListener>
void tHierarchicalClustering::cluster(fListener listener) {   
   functiontracker ft("cluster");
   tMergeEvent event;
   while((event = nextMerge()))
      listener(event);
}



/**
 * removeEdge(...)
 */
inline
void tHierarchicalClustering::removeEdge(tEdgeId edgeId) {
//  cerr << "removing edge: " << edges[edgeId] << endl;
   
#ifdef ALLOW_PSI_EDGES
   edges[edgeId].distance_lb = edges[edgeId].distance_ub = 2.0 * ub_unknown;
#else
   edges[edgeId].distance_lb = edges[edgeId].distance_ub = ub_unknown;
#endif
   
  edges[edgeId].lowCluster = 0;
  ubEdgeHeap.changed(edgeId);
  lbEdgeHeap.changed(edgeId);
}

#endif
