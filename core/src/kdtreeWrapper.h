/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

#ifndef _GIMLI_KDTREEWRAPPER__H
#define _GIMLI_KDTREEWRAPPER__H

#include "gimli.h"

#include <iostream>

#ifdef _MSC_VER
	// warning from kdtree
	#pragma warning( disable: 4396)
#endif

//#define KDTREE_DEFINE_OSTREAM_OPERATORS
#include <kdtree++/kdtree.hpp>

struct _Node_Bracket_accessor { 
    typedef double result_type;
    result_type operator()(GIMLI::Node * n, size_t i) const;
};

typedef KDTree::KDTree< 3, GIMLI::Node *, _Node_Bracket_accessor > NodeKDTree;

namespace GIMLI{

//! Interface class for a kd-search tree. We use it for fast nearest neighbor point search in three dimensions.
/*! Interface class for a kd-search tree. Currently we use libkdtree++ from: http://libkdtreeplus-pplus-p.sourcearchive.com/
We use it for fast nearest neighbor point search in three dimensions. The tree is designed to cooperate with \ref Mesh thus it has to be feeded by pointers of \ref Node. */
class KDTreeWrapper{
public:
    /*! Standard constructor */
    KDTreeWrapper();

    /*! Standard destructor */
    ~KDTreeWrapper();

    /*! Insert new node to the tree */
    void insert( Node * node );

    /*! Find the nearest \ref Node to the coordinates pos. */
    Node * nearest( const RVector3 & pos );

    /*! Return the amount of nodes inside the tree. */
    uint size( ) const;

    /*! Return a pointer to the base libkdetree++  */
    NodeKDTree * tree();

protected:
    NodeKDTree * tree_;
};

} // namespace GIMLI

#endif // _GIMLI_KDTREEWRAPPER__H
