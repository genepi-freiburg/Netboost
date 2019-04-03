/**
 * C++ replacement for R Dist_TOM() function (from schlosser@imbi...)
 *
 * Should be sufficient fast enough, although could be optimized further
 * (either store more data precalculated and then use C-functions instead of
 * STL (there is far too much copying involved atm)
 *
 * Uses RcppParallel for parallel execution.
 *
 * Requires: C++11
 *
 * Jun 2016, jo (jo@imbi.uni-freiburg.de)
 */

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::interfaces(r)]]

#include <R.h>
#include <Rcpp.h>
#include <RcppParallel.h>

#include <cstddef>
#include <vector>
#include <array>

// Using the Rcpp namespace, Rcpp::NumericVector can be written as NumericVector
using namespace Rcpp;
using namespace std;

// Types of R (usually int and double, but could be other)
typedef Sint Rint;
typedef Sfloat Rfloat;

/**
 * Main calculation of the distance if all partners are found.
 *
 * Function is inlined, so no overhead for calling.
 */
inline double calc_distance(const vector<Rint> &i1i2,
                            const vector<Rint> &j1j2,
                            const vector<Rint> &i_partners,
                            const vector<Rint> &j_partners,
                            const vector<Rint> &merges,
                            const NumericVector &adjacency,
                            Rint k) {
  /**
   * R:
   *  prodsum <- sum(sapply(1:length(merges),function(m){
   *     imerge <- i1i2[which(i_partners == merges[m])]
   *     jmerge <- j1j2[which(j_partners == merges[m])]
   *
   *     return(adjacency[imerge]*adjacency[jmerge])
   *  }))
   */
  double prodsum = 0;
  
  if (merges.size() > 0) {
    double adja_pos_i = .0;
    double adja_pos_j = .0;
    
    for (auto merge: merges) {
      // @todo In the given data, no duplicates of "merge" in i_partners are
      //       existing, but this may be the cause for other datasets
      //       Solve with stored indexes and loop until one of the arrays is
      //       looped through.
      for (decltype(i_partners.size()) idx = 0; idx < i_partners.size(); idx++) {
        if (i_partners[idx] == merge) {
          adja_pos_i = adjacency[i1i2[idx]];
          break;
        }
        
        if (idx == i_partners.size()) {
          Rprintf("No value in i_partners for merge m=%d (k=%d)\n", merge, k);
        }
      }
      
      for (decltype(j_partners.size()) idx = 0; idx < j_partners.size(); idx++) {
        if (j_partners[idx] == merge) {
          adja_pos_j = adjacency[j1j2[idx]];
          break;
        }
        
        if (idx == j_partners.size()) {
          Rprintf("No value in j_partners for merge m=%d (k=%d)\n", merge, k);
        }
      }
      
      prodsum += (adja_pos_i * adja_pos_j);
    }
  }
  
  /**
   * R: (new code, mail Pascal 22.06.2016):
   *   min(sum(adjacency[i1i2]),sum(adjacency[j1j2])) + 1 - adjacency[k])
   */
  double sum_i1i2 = 0;
  double sum_j1j2 = 0;
  
  for (auto i: i1i2)
    sum_i1i2 += adjacency[i];
  
  for (auto j: j1j2)
    sum_j1j2 += adjacency[j];
  
  double minsum = min(sum_i1i2, sum_j1j2) + (1.0 - adjacency[k]);
  
  // No merges: prodsum == 0
  double distance = 1.0 - ((prodsum + adjacency[k]) / minsum);
  
  return distance;
}

/**
 * Simple cache class to hold both index and partner arrays.
 * The cache is stored in two arrays later (left and right).
 * (basically only a struct with constructor)
 *
 * Data types must be vectors, as index in both arrays must be identical.
 *
 * @TODO Add sorted partners array (for search of merge)
 */
template <typename T>
class cache
{
public:
  // List of row indizes (in filter matrix)
  std::vector<T> *index;
  
  // List of node values in opposite column
  std::vector<T> *partners;
  
  cache() {
    index = new std::vector<T>();   // Indizes filter-rows
    
    // Values of opposite filter column (if indizes of left are stored, the
    // values of the right column are stored here and vice versa)
    partners = new std::vector<T>();
    
    // Reserve space instantly so copying on resizing of the vectors should
    // kept to minimum (size 50 should be sufficient for most items)
    index->reserve(50);
    partners->reserve(50);
  }
};

/**
 * Helper-structure for RcppParallel::parallelFor (inherits from Worker,
 * which does the trick).
 * The parallel task is called for a given slice and is implemented in the
 * "()" operator.
 */
template <typename T>
struct Distance_Parallel : public RcppParallel::Worker {
  // Input matrix (from R program): will be only read, so does not have to be
  // transformed to synchronized versions
  IntegerMatrix filter;
  NumericVector adjacency;
    
  // Caches for left and right values
  std::vector<cache<T>> all_i;
  std::vector<cache<T>> all_j;
  
  // Results (synchronised, as written by all threads)
  RcppParallel::RVector<double> results;
  
  // Initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Distance_Parallel(const IntegerMatrix filter, const NumericVector adjacency,
                    const vector<cache<T>> all_left, const vector<cache<T>> all_right,
                    NumericVector results) :
    filter(filter), adjacency(adjacency),
    all_i(all_left), all_j(all_right),
    results(results) {
  }
  
  // Function call operator that work for the specified range (begin/end)
  // @todo Still a bit too much copying and vector ops.
  void operator()(std::size_t begin, std::size_t end) {
    // @todo Check if user interrupts work in parallel functions at all.
    //    Rcpp::checkUserInterrupt();
    for (size_t k = begin; k < end; k++) {
      auto col_1 = filter(k, 0);
      auto col_2 = filter(k, 1);
      
      // Copy all vectors from the index (cannot be modified in-place, as
      // in use maybe several times.
      // (Note: do not use "=" operator, as this would trigger move (although
      // working with atomics -- as those are always copied)
      // Create c(i1,i2) by adding left index cache [col_1] to i1
      // Same for c(j1,j2)
      std::vector<Rint> i1i2(*(all_i[col_1].index));   // copy i1
      std::copy(all_j[col_1].index->begin(),
                all_j[col_1].index->end(),
                back_inserter(i1i2));
      
      std::vector<Rint> j1j2(*(all_i[col_2].index));   // copy j1
      std::copy(all_j[col_2].index->begin(),
                all_j[col_2].index->end(),
                back_inserter(j1j2));
      
      // Create c(filter[i1,2],filter[i2,1]) (values of filters must not be taken
      // from the matrix, but are stored in the partner arrays of the cache)
      std::vector<Rint> i_partners(*(all_i[col_1].partners));
      std::copy(all_j[col_1].partners->begin(),
                all_j[col_1].partners->end(),
                back_inserter(i_partners));
      
      std::vector<Rint> j_partners(*(all_i[col_2].partners));
      std::copy(all_j[col_2].partners->begin(),
                all_j[col_2].partners->end(),
                back_inserter(j_partners));
      
      // // Merging with the unprepared sets is quite unconfortable: i_partners and
      // // j_partners have to be sorted first (and therefore copied)
      // std::vector<int> merge_i = i_partners;             // Copy
      // std::vector<int> merge_j = j_partners;
      
      // // Sorting is required for fast merge (and esp. to use STL function for
      // // intersection at all)
      // // @TODO Sorting could be done upfront on building the cache (maybe not
      // // faster, as this increases the sequential part)
      // std::sort(merge_i.begin(), merge_i.end());    // Sort merge_i
      // std::sort(merge_j.begin(), merge_j.end());
      
      // // Result vector needs size for both vectors max (whyever twice)
      // vector<int> merges(merge_i.size() + merge_j.size());
      
      // // Intersection on the sorted vectors
      // std::vector<int>::iterator it = std::set_intersection(merge_i.begin(),
      // 							    merge_i.end(),
      // 							    merge_j.begin(),
      // 							    merge_j.end(),
      // 							    merges.begin());
      
      // // Resize result (skip trailing zeros)
      // merges.resize(it - merges.begin());
      
      // Manually find merges without sorting, but brute force.
      // This is faster than the above version (with sorting and intersect) if
      // source vectors are not too big (which is the case for test data). For
      // bigger data sets, the correct merging is most likely faster than this
      // (empirical: one or two merges)
      // Note: this brute search only wipes duplicates in j_partners, not
      // i_partners (c(1,2,3),c(2,2,2) => c(2), but c(2,2,2),c(1,2,3)=>c(2))
      std::vector<Rint> merges;
      //      merges.reserve(10);
      
      for (auto &i: i_partners) {
        for (auto &j: j_partners) {
          if (i == j) {
            merges.push_back(i);
            break;
          }
        }
      }
      
      // Calc distance value.
      results[k] = calc_distance(i1i2, j1j2, i_partners, j_partners,
                                 merges, adjacency, k);
    }
  }
};

//  @export
// TODO Preparation of cache could also be done in parallel
//' @title Function to calcutate distance
//' @details
//' Steps:
//'   1. - Sequential preparation of index and partner caches per value in filter
//'   2. - Parallel calculation of the distances with cached vectors
//'
//' @param filter Filter matrix
//' @param adjacency Vector
//' @return numeric vector
// [[Rcpp::export(name = "cpp_dist_tom")]]
NumericVector dist_tom(const IntegerMatrix &filter,
                       const NumericVector &adjacency) {
  // A full loop to fetch maximum of values (max may be far smaller than
  // nrow(filter) or even bigger)
  // Note, this reserves one spot for each [0:max[, so if only values 1 and
  // 888 are used, still 889 spots are reserved.
  // (max also delivers the type of integer values in filter).
  auto max = *(std::max_element(filter.begin(), filter.end())) + 1;
  
  // Reserve cache for all left values and all right values. Elements of the
  // cache are initialised with default constructor (so vectors inside the 
  // cache entries are already available).
  std::vector<cache<Rint>> all_i(max);
  std::vector<cache<Rint>> all_j(max);
  
  // Iterate filter once to create index of values on the left and right side
  // and also save the accordings partners.
  // With this cache later the required information for calculating can be
  // restored much faster.
  for (decltype(filter.nrow()) r = 0; r < filter.nrow(); r++) {
    auto col_1 = filter(r, 0);
    auto col_2 = filter(r, 1);
    
    all_i[col_1].index->push_back(r);
    all_i[col_1].partners->push_back(col_2);
    
    all_j[col_2].index->push_back(r);
    all_j[col_2].partners->push_back(col_1);
  }
  
  // Result vector
  NumericVector results(filter.nrow());
  
  // Init structure for parallel processing
  Distance_Parallel<Rint> distPara(filter, adjacency, all_i, all_j, results);
  
  try {
    // Parallel calculation.
    // Do not use less than 10000 rows to feed single thread (so thread handling
    // overhead does not outweight pure calculation)
    // Interval [0, nrow()[
    RcppParallel::parallelFor(0, filter.nrow(), distPara, 10000);
  }
  // STL exceptions
  catch (const std::exception &e) {
    forward_exception_to_r(e);
  }
  // Others.
  catch (...) {
    //     auto eptr = std::current_exception();
    Rcout << endl << endl << "C++ exception. " << endl << endl;
    ::Rf_error("(unknown reason)");
  }
  
  return results;
}

/**
 Note after bugfix (mail schlosser 22.06.2016):
 (min(sum(adjacency[i_partners]),sum(adjacency[j_partners]))
 is now (in C++):
 (min(sum(adjacency[c(i1,i2)]),sum(adjacency[c(j1,j2)]))
 
 Original R-code was (Mail schlosser 17.06.2016):
 
 library(parallel)
 load("jochen.Rdata")
 
#### Calculation of network adjacencies for filter ####
 tst <- Sys.time()
 adjacency <- sapply(1:nrow(filter),function(i){abs(cor(data[,filter[i,1]],data[,filter[i,2]]))^beta})
 Sys.time()-tst
 
 tst <- Sys.time()           
 DistTOM <- unlist(mclapply(1:nrow(filter),function(k){
 i1 <- which(filter[,1]==filter[k,1])
 i2 <- which(filter[,2]==filter[k,1])
 j1 <- which(filter[,1]==filter[k,2])
 j2 <- which(filter[,2]==filter[k,2])
 i_partners <- c(filter[i1,2],filter[i2,1])
 j_partners <- c(filter[j1,2],filter[j2,1])
 merges <- intersect(i_partners,j_partners)
 if(length(merges)>0){
 prodsum <- sum(sapply(1:length(merges),function(m){
 imerge <- c(i1,i2)[which(i_partners == merges[m])]
 jmerge <- c(j1,j2)[which(j_partners == merges[m])]
 return(adjacency[imerge]*adjacency[jmerge])
 }))
 return((1 -(
 (prodsum + adjacency[k])
 /(min(sum(adjacency[i_partners]),sum(adjacency[j_partners])) + 1 - adjacency[k])
 )
 )
 )
 } else {
 return((1-(adjacency[k]/
 (min(sum(adjacency[i_partners]),sum(adjacency[j_partners])) + 1 - adjacency[k]))))
 }
 }))
 Sys.time()-tst
 
 */
