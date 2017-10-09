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
 * file: IteratorAdaptor.hpp 
 *
 * Yaniv Loewenstein, 2007
 *
 * A different interface for an STL iterator.
 *
 ***********************************************************************/


#ifndef _ITERATOR_ADAPTOR_HPP
#define _ITERATOR_ADAPTOR_HPP
#include "Edge.hpp"
//#include <iterator>
using namespace std;
template <typename T>

/**
 * Supply a next() and has_next() interface to an STL istream_iterator to replace
 * an older implementation (without STL).
 */
class IteratorAdaptor {
   typedef istream_iterator<T> InputIteratorT;
protected:
   InputIteratorT _begin,_end;
public:
   IteratorAdaptor(istream & in) : _begin(InputIteratorT(in)) 
   { ; }
   IteratorAdaptor(InputIteratorT begin,InputIteratorT end) : _begin(begin), _end(end)
   { ; }  
   inline bool has_next() const {
      return _begin != _end;
   }	
   T next() {
      assert(has_next());
      return *(_begin++);
   }
   T peak_next() {
      return *_begin;
   }
   InputIteratorT begin() {
      return _begin;
   }
   InputIteratorT end() {
      return _end;
   }
};
#endif // _ITERATOR_ADAPTOR_HPP
