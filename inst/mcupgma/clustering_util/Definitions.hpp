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
 * file: Definitions.hpp 
 *
 * Yaniv Loewenstein, 2007
 *
 * Global defitions to be shared across other header files 
 *
 ***********************************************************************/

#ifndef _EDGEREADER2_HPP
#define _EDGEREADER2_HPP

#include "Edge.hpp"
#include "IteratorAdaptor.hpp"

#define VERBOSE 0

typedef IteratorAdaptor<Edge *>        EdgePtrReader;
typedef IteratorAdaptor<WeightedEdge > WeightedEdgeReader;

typedef unsigned int uint;
#endif //_EDGEREADER2_HPP
 
