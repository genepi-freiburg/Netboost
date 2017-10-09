/*****************************************************************************
 * MC-UPGMA  - Accurate huge scale clustering by Memory Constrained UPGMA    *
 *             Loewenstein et al. Bioinformatics. 2008 Jul 1;24(13):i41-9.   *
 *                                                                           *
 * test_heap.cpp - Heap driver and tester                                    *
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


using namespace std;

#include "test_heap.cmdline.h"
#include "Heap.hpp"
#include "myutils/StreamFromFileNameGenerator.hpp"
#include "myutils/Exception.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <stdlib.h>

int main(int argc, char** argv) {
  try {
    gengetopt_args_info args_info;
    if (cmdline_parser (argc, argv, &args_info) != 0)
      throw tException();

    if (args_info.rand_seed_given)
      srand(args_info.rand_seed_arg);
    else
      srand(time(0));

    vector<long> v;
    multiset<long> s;
    for (int i = 0; i < args_info.init_num_numbers_arg; ++i) {
      long k = static_cast<long>((static_cast<float>(rand()) / RAND_MAX) * args_info.range_arg);
      v.push_back(k);
      s.insert(k);
    }
    tHeap<vector<long>::iterator, greater<long> > h(v.begin(), v.size());
    int iteration = 0;
    try {
      while(v[h.top()] < args_info.range_arg) {
	h.verify();
	if (s.empty())
	  throw (tException() << "multiset empty before heap at iteration " << iteration);
	if (v[h.top()] != *(s.begin()))
	  throw (tException() << "mismatch as iteration " << iteration << " " << v[h.top()] << ' ' << *(s.begin()));
	if ((s.size() > 1) && (v[h.second_top()] != *(++s.begin())))
	  throw (tException() << "second_top mismatch as iteration " << iteration << " " << v[h.second_top()] << ' ' << *(++s.begin()));
	if (args_info.break_at_given && (iteration == args_info.break_at_arg))
	  throw(tException() << "requested break at iteration " << iteration);
	s.erase(s.begin());

	bool reentry = (static_cast<float>(rand()) / RAND_MAX < args_info.reentry_arg);
	long k;
	if (reentry)
	  k = static_cast<long>((static_cast<float>(rand()) / RAND_MAX) * args_info.range_arg);
	else
	  k = args_info.range_arg;
	v[h.top()] = k;
	h.changed(h.top());
	if (reentry)
	  s.insert(k);
	++iteration;
      }
      cout << "heap is empty" << endl;
      h.verify();
      if (! s.empty())
	throw (tException() << "heap empty before multiset at iteration " << iteration);
    } catch (tException e) {
      cerr << "v contents:" << endl;
      copy(v.begin(), v.end(), ostream_iterator<long>(cerr, " "));
      cerr << endl;
      cerr << "s contents:" << endl;
      copy(s.begin(), s.end(), ostream_iterator<long>(cerr, " "));
      cerr << endl;
      cerr << "heap trace:" << endl;
      for (int i = 0; i != v.size(); ++i) {
	cerr << "heap trace " << i << "\n";
	h.trace(cerr, i);
	cerr << "\n";
      }
      throw;
    }
  } catch (tException const & e) {
    cerr << e << endl;
    return 1;
  } catch (std::exception const & e) {
    cerr << "caught std::exception" << endl;
    cerr << e.what() << endl;
    return 1;
  } catch (...) {
    cerr << "exception thrown" << endl;
    return 1;
  }
  return 0;
}
