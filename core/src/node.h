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

#ifndef _GIMLI_NODE__H
#define _GIMLI_NODE__H

#include "gimli.h"

#include "baseentity.h"
#include "pos.h"

#include <set>

namespace GIMLI{

enum NodeState{No, Original, Secondary, Connected};

//! 3D Node
/*!
 * Node is a basic entity of a mesh at a 3D position x/y/z (a vertex),
 * which owns a marker, an id, and information about connected boundarys and cells.
 * For a 2D node ignore $z$ or let $z$-coord = 0.0.
 */
class DLLEXPORT Node : public BaseEntity {

public:
    /*! Construct node with non valid position and marker = 0 */
    Node();

    /*! Construct node from coordinates $x,y,[z]$ with marker = 0 */
    Node(double x, double y, double z=0.0);

    /*! Construct node from RVector3 */
    Node(const RVector3 & pos);

    /*! Construct node from RVector3 with marker and optional id */
    Node(const RVector3 & pos, int marker, int id=-1);

    /*! Copy constructor */
    Node(const Node & node);

    /*! Assignement operator */
    Node & operator = (const Node & node);

    /*! Destruct the node and all containing informations */
    ~Node();

    /*! Unchecked index operator to pos */
    inline double & operator [] (const Index i) { return this->at(i); }

    /*! Unchecked index operator to pos */
    inline const double & operator [] (const Index i) const { return this->at(i); }

    inline uint rtti() const { return MESH_NODE_RTTI; }

    inline void setPos(const RVector3 & pos) { changed_(); pos_ = pos; }

    inline const RVector3 & pos() const { return pos_; }

    /*! Don't use this if you plan to modify the nodes posisiton. The mesh won't get any notifications of the change and this not remove geometry caches. */
    inline RVector3 & pos() { return pos_; }

    inline void insertBoundary(Boundary * bound){ boundSet_.insert(bound); }

    inline void eraseBoundary(Boundary * bound){ boundSet_.erase(bound); }

    inline void insertCell(Cell * cell){ cellSet_.insert(cell); }

    inline void eraseCell(Cell * cell){ cellSet_.erase(cell); }

    inline const std::set < Boundary * > & boundSet() const { return boundSet_; }

    inline std::set < Boundary * > & boundSet() { return boundSet_; }

    inline const std::set < Cell * > & cellSet() const { return cellSet_; }

    inline std::set < Cell * > & cellSet() { return cellSet_; }

    void transform(const RMatrix & mat);

    inline void scale(const RVector3 & s) { changed_(); pos_.scale(s); }

    inline void translate(const RVector3 & t) { changed_(); pos_.translate(t); }

    inline void translate(double x, double y=0.0, double z=0.0) {
        changed_(); pos_.translate(x, y, z); }

    inline void rotate(const RVector3 & r) { changed_(); pos_.rotate(r); }

    inline void swap(Index i, Index j) { changed_(); pos_.swap(i, j); }

    inline const double & x() const { return pos_[0]; }

    inline const double & y() const { return pos_[1]; }

    inline const double & z() const { return pos_[2]; }

    inline double & at(Index i) { return pos_[i]; }

    inline const double & at(Index i) const { return pos_[i]; }

    inline double dist(const Node & n) const { return pos_.dist(n.pos()); }

    /*!*/
    void smooth(uint function);

    /*!Little helper to identify the state of this node after some merging.*/
    void setState(NodeState s) { this->_state = s; }

    /*!Return the state of this node.*/
    const NodeState state() const { return this->_state; }

    void setSecondaryParent(MeshEntity * e) { this->_secondaryParent = e; }
    /*!Return the state of this node.*/
    MeshEntity * secondaryParent() { return this->_secondaryParent; }

protected:

    void copy_(const Node & node);

    void changed_();

    void init_();

    RVector3 pos_;

    std::set < Boundary * > boundSet_;
    std::set < Cell * >     cellSet_;

    NodeState _state;
    MeshEntity *_secondaryParent;
}; // class Node

DLLEXPORT std::ostream & operator << (std::ostream & str, const GIMLI::Node & node);

DLLEXPORT std::ostream & operator << (std::ostream & str, const std::vector < GIMLI::Node * > & nodes);

inline bool operator == (const Node & n1, const Node & n2) { return n1.pos() == n2.pos(); }

} // namespace GIMLI

#endif // _GIMLI_NODE__H
