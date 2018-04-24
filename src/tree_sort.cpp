/**
 * Treesorting implementation according to the R program by Pascal Schlosser.
 *
 * Standard: C++11.
 *
 * Sept/Oct 2016, jo@imbi.uni-freiburg.de
 */

// [[Rcpp::interfaces(r)]]

#include <R.h>
#include <Rcpp.h>

#include <cstddef>
#include <vector>
#include <array>
#include <limits>

#include <unordered_map>
#include <algorithm>
#include <memory>

#include "./rdtsc.h"             // Simple cycle counter

// Using the Rcpp namespace, Rcpp::NumericVector can be written as NumericVector
using namespace Rcpp;

constexpr bool SECURE = true;    // Range check on cache-set (minimal overhead)
constexpr bool DEBUG = false;    // Debugging output

template <typename T> class Tree;

using STORAGE_INT = int;         // Type of indexing integers. May exceed 32-bit.

/**
 * Storage type for IDs: either set or vector. Cache class has specialisation
 * for both containers, so those can be exchanged smoothly.
 *
 * Set:    - always unique and return sorted
 *         - smaller memory footprint (no duplicates)
 *         - very fast copy at the end
 * Vector: - very fast insert (as simply appended)
 *         - higher memory footprint (duplicates stored)
 *         - slower read at the end (as sorting/unique required)
 */
// using STORAGE_TYPE = std::set<STORAGE_INT>;
using STORAGE_TYPE = std::vector<STORAGE_INT>;

// Convenience shortcut.
using TREE_ST = Tree<STORAGE_TYPE>;

/**
 * Small helper class, basically wrapping access to a static cache array.
 * (here implemented as singleton).
 */
template<typename STORAGETYPE>
class Node_Cache {
private:
  // Cache-Map to determine tree directly from any ID.
  // Static accross all tree-objects.
  std::vector<Tree<STORAGETYPE>*> cache_idx_tree;

public:
  inline static Node_Cache *instance() {
    static Node_Cache *instance = new Node_Cache();
    return instance;
  }

  /**
   * Initialise cache: fill whole container with NULL.
   *
   * @param[in] size Amount of elements
   */ 
  void init(const size_t size) {
    cache_idx_tree.clear();
    cache_idx_tree.assign(size + 1, nullptr);
  }

  /**
   * Set cache element in given tree.
   * Due to constraint, maybe used without bounding check (if init() was called)
   * 
   * @param[in] val Cache-key
   * @param[in] tree (Raw) pointer to tree
   */
  inline void set(const size_t val, Tree<STORAGETYPE> *tree) {
    if (SECURE) {
      if (val >= 0 && (val < cache_idx_tree.size()))
        cache_idx_tree[val] = tree;
      else {
        Rcout << "SET: Accessing outside of cache. Elem: " << val
              << " SIZE: " << cache_idx_tree.size()
              << " CAP: " << cache_idx_tree.capacity()
              << std::endl;
        Rf_error("Cache set: outside access");
      }
    }
    else
      cache_idx_tree[val] = tree;
  }

  /**
   * Get cache value. Due to constraint, maybe used without bounding check.
   * 
   * @param[in] val Cache-key
   * @return Cache element (nullptr if not used atm)
   */ 
  inline Tree<STORAGETYPE>* get(const size_t val) {
    if (SECURE) {
      if (val >= 0 && (val < cache_idx_tree.size()))
        return cache_idx_tree[val];

      Rcout << "GET: Accessing outside of cache. Elem: " << val << std::endl;

      Rf_error("Cache get: outside access");

      return nullptr;
    }
    else
      return cache_idx_tree[val];
  }

  /**
   * Debugging: get list of unique Tree-pointers.
   *
   * @todo Currently also nullptr will be delivered once.
   * @return Vector with unique tree-pointer (+nullptr)
   */
  std::vector<Tree<STORAGETYPE>*> unique() {
    // Copy container (as trashed by inplace-sort).
    std::vector<Tree<STORAGETYPE>*> v(cache_idx_tree);

    // std::unique requires sorted container.
    std::sort(v.begin(), v.end());

    auto last = std::unique(v.begin(), v.end());  // Unique to front container
    v.erase(last, v.end());                       // Delete duplicates at the end

    return v;
  }
};

/**
 * Class to store a single tree.
 * This contains the IDs of the contained nodes and the row-indizes of the
 * input matrix.
 * 
 * @tparam STORAGETYPE Container-type to store the IDs and row-index.
 */
template<typename STORAGETYPE>
class Tree {
  // Originally as template parameter, but removed for global type.
  using T = STORAGE_INT;

private:
  unsigned int id;                  //!< (Internal) ID of the tree
  static unsigned int id_running;   //!< ID of last inserted tree
  static unsigned int trees;        //!< (Debugging) Amount of trees at the moment

public:
  /**
   * As boths lists (IDs and row numbers) are required to be unique, they are
   * stored as set and therefore are only inserted if not currently in the set
   * -- so making them unique afterwards is not required.
   * (as this is a binary tree, much faster search).
   */
  STORAGETYPE ids;
  // std::set<T> rows;    // Mind: 0/1 based. || @todo As each rowno is unique, could be vector
  std::vector<T> rows;    // Mind: 0/1 based. || @todo As each rowno is unique, could be vector

  /*!
   * Constructor: MUST be unique, so no merge.
   *
   * @param[in] id1 Node 1
   * @param[in] id2 Node 2
   * @param[in] id3 Node 3
   * @param[in] row Row number in original data matrix
   */
  Tree(T id1, T id2, T id3, T row) :
    ids({id1, id2, id3}),       // Add all IDs to set
    rows({row + 1}),            // Directly to 1-based
    id(++Tree::id_running)      // Set ID (and increment globally)
  {
    auto cache = Node_Cache<STORAGE_TYPE>::instance();

    cache->set(id1, this);      // Set cache tovalues to new tree.
    cache->set(id2, this);
    cache->set(id3, this);

    ++Tree::trees;              // Debugging: increase amount of trees
  }

  /*!
   * Depending on storage type for the IDs, this method will be specialized
   * later.
   * Default implementation just throws error. Must be specialized for each
   * possible storage type (set or vector).
   *
   * @param[in] id ID to add to this tree 
   */
  inline void insert_id(T id) {
    Rf_error("Unspecialized call to insert_id()");
  }
  
  // ~Tree() {
  //   Rcout << "Destroy " << id << " (curr max: " << Tree::id_running << ")\n";
  // }

  void add(T id1, T id2, T id3, T row) {
    auto cache = Node_Cache<STORAGE_TYPE>::instance();

//    ids.insert({id1, id2, id3});
    insert_id(id1);             // Single calls, as STORAGE_TYPE can be vector
    insert_id(id2);             //  or set.
    insert_id(id3);

    cache->set(id1, this);
    cache->set(id2, this);
    cache->set(id3, this);

    rows.push_back(row + 1);    // Add row number
  }
  
  /*!
   * @brief Merge two trees (add argument to this).
   *
   * set_union is not inplace.
   *
   * If there are too many trees deleted, memory could be freed here.
   * (Currently this is done via insertion in the linked list in the
   * main method).
   *
   * @param[in] add_tree Tree to merge (destroyed!)
   */
  void merge(const Tree &add_tree) {
    auto cache = Node_Cache<STORAGE_TYPE>::instance();

    for (const auto &elem: add_tree.ids) {
      // Insert ID in tree
      insert_id(elem);

      // Update element cache to new merged object
      cache->set(elem, this);
    }

    // Merge rows (simply insert/append, as both lists are already unique
    // (due to program logic). Result is not sorted of course.
    const auto at_rows = add_tree.rows;

    rows.reserve(rows.size() + at_rows.size());
    rows.insert(rows.end(), at_rows.begin(), at_rows.end());

    // Debugging: remove amount of existing trees (Note: deleting of
    // merged tree in caller!)
    --Tree::trees;
  }

  /**
   * Check if given ID is present in the tree presentation.
   *
   * @param[in] val ID
   * @return In current IDs
   */
  inline bool check_id(const STORAGE_INT val) {
    return ids.count(val) == 1;
  }

  inline int get_id() {
    return id;
  }

  /**
   * Return vector with row numbers for this tree.
   * @todo Sorting may be unrequired.
   *
   * @return Vector with row numbers (R notation, so 1-based)
   */
  inline const Rcpp::IntegerVector get_rows() {
    std::sort(rows.begin(), rows.end());
    return Rcpp::IntegerVector(rows.begin(), rows.end());
  }

  /**
   * Get IDs of this tree (return is: sorted, unique).
   * Must be specialized depending on used storage type (set or vector).
   *
   * @return Vector with sorted/unique IDs.
   */
  Rcpp::IntegerVector get_ids() {
    Rf_error("Unspecialized call for get_ids()");
  }

  // Debugging.
  void print() {
    Rcout << this << "  IDs: ";

    for (const auto &e: ids)
      Rcout << e << ", ";

    Rcout << std::endl << "ROWs: ";

    for (const auto &r: rows)
      Rcout << r << ", ";

    Rcout << std::endl;
  }

  // Debugging: current number of trees.
  static int no_trees() {
    return Tree::trees;
  }
};

/**
 * Insert into ID storage: specialization for Set-storage.
 */
template <> inline void Tree<std::set<STORAGE_INT>>::insert_id(STORAGE_INT id) {
  ids.insert(id);
}

/**
 * Insert into ID storage: specialization for Vector-storage.
 */
template <> inline void Tree<std::vector<STORAGE_INT>>::insert_id(STORAGE_INT id) {
  ids.push_back(id);
}

/**
 * Get IDs of this tree: specialization for Set-storage.
 * As set is already unique and returns to sorted, this is quite trivial.
 */
template <> Rcpp::IntegerVector Tree<std::set<STORAGE_INT>>::get_ids() {
  return Rcpp::IntegerVector(ids.begin(), ids.end());
}

/**
 * Get IDs of this tree: specialization for Vector-storage.
 * As duplicates are in the vector, sorting and making unique is done first.
 * Done inplace, so object is altered (attribute "ids").
 */
template <> Rcpp::IntegerVector Tree<std::vector<STORAGE_INT>>::get_ids() {
  // Unique values (must be sorted first). Inplace.
  std::sort(ids.begin(), ids.end());
  auto last = std::unique(ids.begin(), ids.end());
  ids.erase(last, ids.end());

  return Rcpp::IntegerVector(ids.begin(), ids.end());
}

// Init static id/counter.
template<typename T> unsigned int Tree<T>::id_running = 0;
template<typename T> unsigned int Tree<T>::trees = 0;

//  export  (NO export here, wrapper in R required)
//' @title Tree search.
//' @name cpp_tree_search
//' 
//' @description
//' Constraint: IDs 0 <= x (Integer)
//' 
//' @backref src/tree_sort.cpp 
//'
//' @param netboost_forest Input-matrix (4 columns, ids in colum 0,1,3)
//'
// [[Rcpp::export(name = "cpp_tree_search")]]
Rcpp::List tree_search(const IntegerMatrix &netboost_forest) {
  auto start = cycles();
  
  // Get max. value inside the three number colums.
  auto max1 = std::max_element(netboost_forest.column(0).begin(),
			       netboost_forest.column(0).end());
  auto max2 = std::max_element(netboost_forest.column(1).begin(),
			       netboost_forest.column(1).end());
  auto max3 = std::max_element(netboost_forest.column(3).begin(),
			       netboost_forest.column(3).end());

  size_t max = std::max(*max1, std::max(*max2, *max3));

  if (DEBUG)
    Rcout << "MAX: " << max << " - " << max << " COLMAX: " <<
      *max1 << ", " << *max2 << ", " << *max3 << std::endl;

  // "Storage" for generated trees. Only fullfills one purpose: keep
  // the stack generated elements undeleted inside the function...
//  std::list<std::unique_ptr<Tree>> trees;
  std::list<TREE_ST*> trees;
  
  const auto type_in = netboost_forest(0,0);

  unsigned int merge_delete = 0;
  
  // Init Node cache (basically an array with [0,max] elements).
  auto cache = Node_Cache<STORAGE_TYPE>::instance();
  cache->init(max);

  for (auto rowno = 0; rowno < netboost_forest.nrow(); rowno++) {
    auto row = netboost_forest.row(rowno);

    // Automatically unique set (if two nodes share two IDs).
    // Using a set to remove the duplicates of just three values is overkill
    // (but used for convenience here).
    std::set<TREE_ST*> hits;

    // Check all columns: if any value was already used (in any tree, mark
    // this subtree).
    for (auto &col: {0, 1, 3}) {
      auto node = cache->get(static_cast<size_t>(row[col]));

      if (node != nullptr) {
        hits.insert(node);
      }
    }

    // No hits: create new tree
    switch (hits.size()) {
    // Create new tree
    case 0:
      // Preserve data structure. Use unique_ptr for implicit destroy.
      //	trees.push_front(std::unique_ptr<Tree>(new Tree(row[0], row[1], row[3], rowno)));
      trees.push_front(new TREE_ST(row[0], row[1], row[3], rowno));
      
      if (DEBUG)
        Rcout << "CPP: NEW TREE: " << trees.front()->get_id() << "\n";
      
      break;
      
      // Add to the only matching tree.
    case 1: 
        {
          auto tree_add_it = hits.cbegin();
          
          (*tree_add_it)->add(row[0], row[1], row[3], rowno);
          
          if (DEBUG)
            Rcout << "CPP: ADD TO TREE: " << (*tree_add_it)->get_id() << "\n";
          
          break;
        }
      // Merge all participating trees (simply into the first tree found)
    case 2:
    case 3:
        {
          auto merge = hits.begin();
          
          // First tree is selected as merge base (at random).
          auto base_tree = *merge;
          
          ++merge;
          
          // Merge all other trees with hits for this row.
          for (; merge != hits.end(); ++merge) {
            Tree<STORAGE_TYPE> *merge_tree = *merge;
            
            // Inside the merge, all node-to-tree-pointers are changed so no
            // pointer is leftover to merge afterwards (in the node-cache).
            base_tree->merge(*merge_tree);
            
            merge_delete++;
            
            if (DEBUG)
              Rcout << "CPP: MERGE TREES: " << base_tree->get_id() << " < " << merge_tree->get_id() << "\n";
            
            // Merged tree (in "merge") could be removed/destroyed here, but
            // as there can only be max(node_element[,c(1,2,4)]) subtrees,
            // memory is not a problem and therefore we keep them wiping
            // automatically at the end of the function.
            trees.remove(merge_tree);
            
            // Rcout << "REMOVE: ";
            // merge_tree->print();
            
            delete merge_tree;
          }
          
          // Current row also added.
          base_tree->add(row[0], row[1], row[3], rowno);
          
          break;
        }
      
    default:
      Rcout << "wtf...\n";
      break;
    }
    
    if (DEBUG)
      Rcout << "CPP: ROW: " << rowno << " trees: " << TREE_ST::no_trees() << std::endl;
  }

  // As we reading from cache (each tree is stored up to three times there), get
  // unique trees.
  auto unique_trees = cache->unique();
  
  if (DEBUG) {
    Rcout << "TREES : " << trees.size() << " (merged: " << merge_delete << ")" << std::endl;
    Rcout << "UNIQUE: " << unique_trees.size() << std::endl;
  }

  // Create list of lists (with COPIES(!)) of vectors.
  std::vector<Rcpp::List> ret;

  // List::create is slow, but usually there should be only a very small number of
  // trees in the end.
  for (const auto &tree: unique_trees) {
    if (tree != nullptr) {
      auto cpy = Rcpp::List::create(Rcpp::Named("ids") = tree->get_ids(),
                                    Rcpp::Named("rows") = tree->get_rows());
      ret.push_back(cpy);
    }
  }

  // Free memory completely (would be automatically if unique_ptrs would be used).
  for (auto rem_tree_cache: trees)
    delete rem_tree_cache;

  return wrap(ret);
}
