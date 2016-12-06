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


#ifndef __Heap_hpp
#define __Heap_hpp

#include <vector>
#include <assert.h>
#include <myutils/functiontracker.hpp>
#include <myutils/Exception.hpp>

//lonshy : RAI - random access iterator

template <typename tElementRAI, typename tElementComparator>
class tHeap {
  typedef typename iterator_traits<tElementRAI>::difference_type tIndex;
  tElementRAI elements;
  tIndex numOfElements;
  tIndex offset;
  tElementComparator comp;
  vector<tIndex> nodes;

  static tIndex pow2roundup(tIndex n);
  static tIndex parent(tIndex n);
  static tIndex sibling(tIndex n);

  void build();
public:
  tHeap(tElementRAI elements, tIndex numOfElements, tElementComparator comp = tElementComparator());

  void reset(tElementRAI elements, tIndex numOfElements);

  tIndex top() const;

  tIndex second_top() const;

  void changed(tIndex i);

  ostream& trace(ostream& os, tIndex i) const;

  void verify() const;

  unsigned int memorysize_in_bytes() const;
};

/**
 * C-tor
 */

template <typename tElementRAI, typename tElementComparator>
tHeap<tElementRAI, tElementComparator>::tHeap(tElementRAI elements, tIndex numOfElements, tElementComparator comp) :
  elements(elements),
  numOfElements(numOfElements),
  offset(pow2roundup(numOfElements) - 1),
  comp(comp),
  nodes(parent(numOfElements + offset) + 2, numOfElements) {
  build();
}
/**
 * reset(elements,numOfElements)
 */

template <typename tElementRAI, typename tElementComparator>
void tHeap<tElementRAI, tElementComparator>::reset(tElementRAI elements, tIndex numOfElements) {
  this->elements = elements;
  this->numOfElements = numOfElements;
  offset = pow2roundup(numOfElements) - 1;
  vector<tIndex> tmp(parent(numOfElements + offset) + 2, numOfElements);
  nodes.swap(tmp);
  build();
}
/**
 * build()
 */
template <typename tElementRAI, typename tElementComparator>
void tHeap<tElementRAI, tElementComparator>::build() {
  nodes[0] = 0;
  tIndex p;
  if (offset > 0)
     for(tIndex i = 0; i != numOfElements; ++i) {
	for(p = parent(i + offset);
	    (nodes[p] == numOfElements) || ((p != 0) && comp(elements[nodes[p]], elements[i]));
	    p = parent(p))
	   nodes[p] = i;
	if ((p == 0) && comp(elements[nodes[p]], elements[i]))
	nodes[p] = i;
     }
  //end outer if
}

/**
 * top()
 */

// lonshy - returns index of the maximum
template <typename tElementRAI, typename tElementComparator>
inline
typename tHeap<tElementRAI, tElementComparator>::tIndex tHeap<tElementRAI, tElementComparator>::top() const {
  return nodes[0];
}

/**
 * changed(i)
 */
template <typename tElementRAI, typename tElementComparator>
void tHeap<tElementRAI, tElementComparator>::changed(tIndex i) {
  if (offset != 0) {
    tIndex mi = i;
    tIndex si = sibling(i + offset) - offset;
    if ((si != numOfElements) && (comp(elements[mi], elements[si])))
      mi = si;
    tIndex p = parent(i + offset);
    while(((nodes[p] != mi) || (nodes[p] == i)) && (p != 0)) {
      nodes[p] = mi;
      si = nodes[sibling(p)];
      if ((si != numOfElements) && (comp(elements[mi], elements[si])))
	mi = si;
      p = parent(p);
    }
    if (p == 0)
      nodes[p] = mi;
  }
}

/**
 * second_top
 */
template <typename tElementRAI, typename tElementComparator>
typename tHeap<tElementRAI, tElementComparator>::tIndex tHeap<tElementRAI, tElementComparator>::second_top() const {
  assert(numOfElements >= 2);
  bool second_element_so_far_initiated = false;
  tIndex second_element_so_far;
  tIndex max_node = nodes[0];
  tIndex sibling_element = sibling(max_node + offset) - offset;
  if (sibling_element != numOfElements) {
    second_element_so_far_initiated = true;
    second_element_so_far = sibling_element;
  }
  for(max_node = parent(max_node + offset); max_node != 0; max_node = parent(max_node)) {
    sibling_element = nodes[sibling(max_node)];
    if ((sibling_element != numOfElements) && ((! second_element_so_far_initiated) || (comp(elements[second_element_so_far], elements[sibling_element])))) {
      second_element_so_far_initiated = true;
      second_element_so_far = sibling_element;
    }
  }
  return second_element_so_far;
}

/**
 * trace(os,i)
 */
template <typename tElementRAI, typename tElementComparator>
ostream& tHeap<tElementRAI, tElementComparator>::trace(ostream& os, tIndex i) const {
  tIndex si = sibling(i + offset) - offset;
  os << i << '\t' << elements[i] << '\t' << si << '\t' << elements[si] <<endl;
  tIndex p = parent(i + offset);
  while(p != 0) {
    si = sibling(p);
    os << "( " << p << '\t' << nodes[p] << '\t' << elements[nodes[p]] << " )( " << si << '\t' << nodes[si] << '\t' << elements[nodes[si]] << " )" << endl;
    p = parent(p);
  }
  os << p << '\t' << nodes[p] << '\t' << elements[nodes[p]] << endl;
  return os;
}
/**
 * verify()
 */
template <typename tElementRAI, typename tElementComparator>
void tHeap<tElementRAI, tElementComparator>::verify() const {
  if (offset != 0) {
    for (tIndex i = 0; i != numOfElements; ++i)
      if (comp(elements[nodes[parent(i + offset)]], elements[i]))
	throw (tException() << "heap verification failed. element: " << i << " value: " << elements[i] << " parent: " << parent(i + offset) << " parent element: " << nodes[parent(i + offset)] << " parent value: " << elements[nodes[parent(i + offset)]]);
    for (tIndex i = 1; i != nodes.size(); ++i)
      if (nodes[i] != numOfElements) {
	tIndex e;
	for (e = nodes[i] + offset; e > i; e = parent(e));
	if (e != i)
	  throw (tException() << "heap verification failed. node: " << i << " element: " << nodes[i] << " not decendent of it. e: " << e);
      }
    for (tIndex i = 1; i != nodes.size(); ++i)
      if ((nodes[i] != numOfElements) && (comp(elements[nodes[parent(i)]], elements[nodes[i]])))
	throw (tException() << "heap verification failed. node: " << i << " element: " << nodes[i] << " value: " << elements[nodes[i]] << " parent: " << parent(i) << " parent element: " << nodes[parent(i)] << " parent value: " << elements[nodes[parent(i)]]);

  }
}
/**
 * memorysize_in_bytes()
 */

template <typename tElementRAI, typename tElementComparator>
unsigned int tHeap<tElementRAI, tElementComparator>::memorysize_in_bytes() const {
  return (sizeof(tIndex) * nodes.size());
}
/**
 * pow2roundup(n)
 */
template <typename tElementRAI, typename tElementComparator>
typename tHeap<tElementRAI, tElementComparator>::tIndex tHeap<tElementRAI, tElementComparator>::pow2roundup(tIndex n) {
  tIndex pi = 0;
  tIndex i;
  for(i = 1;
      i < n;
      i *= 2) {
    assert(pi < i);
    pi = i;
  }
  return i;
}    
/**
 * parent
 */
template <typename tElementRAI, typename tElementComparator>
typename tHeap<tElementRAI, tElementComparator>::tIndex tHeap<tElementRAI, tElementComparator>::parent(tIndex n) {
  return (n - 1) / 2;
}
/**
 * sibling
 */
template <typename tElementRAI, typename tElementComparator>
typename tHeap<tElementRAI, tElementComparator>::tIndex tHeap<tElementRAI, tElementComparator>::sibling(tIndex n) {
  if (n % 2)
    return n + 1;
  else
    return n - 1;
}

  
#endif
