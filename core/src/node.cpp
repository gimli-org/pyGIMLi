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

#include "node.h"

#include "matrix.h"
#include "meshentities.h"
#include "shape.h"

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const GIMLI::Node & n){
    str << "Node: "<< &n << " id: " << n.id() << "\t" << n.pos();
    str << " marker: " << n.marker();
    return str;
}

std::ostream & operator << (std::ostream & str, 
                            const std::vector < GIMLI::Node * > & nodes){
    for (auto *n: nodes) {
        str << n->id() << " "; 
    }
    str << std::endl;
    return str;
}

Node::Node(){
    init_();
    marker_ = 0;
}

Node::Node(double x, double y, double z)
: pos_(RVector3(x, y, z)) {
    init_();
}

Node::Node(const RVector3 & pos)
: pos_(pos) {
    init_();
}

Node::Node(const RVector3 & pos, int marker, int id)
: pos_(pos) {
    init_();
    setMarker(marker);
    setId(id);
}

Node::Node(const Node & node){
    init_();
    copy_(node);
}

Node & Node::operator = (const Node & node){
    if (this != &node){
        copy_(node);
    }
    return *this;
}

Node::~Node(){
    //std::cout << " delete Node " << pos_ << " " << id_ << " at " << this << std::endl;
}
void Node::transform(const RMatrix & mat) {
    this->changed_(); pos_.transform(mat);
}
void Node::changed_(){
    for (auto &b : boundSet_){
        b->changed();
    }
    for (auto &c : cellSet_){
        c->changed();
    }
}

void Node::copy_(const Node & node){

    // std::cout << "copy node from "  << & node << " into " << this << std::endl;
    init_();
    pos_    = node.pos();
    marker_ = node.marker();
    setId(node.id());
}

void Node::init_(){
    setId(-1);
    _state = No;
    _secondaryParent = 0;
}

void Node::smooth(uint function){
    std::set< Node * > common(commonNodes(this->boundSet()));
    //** Achtung konkave gebiete koennen entstehen wenn zu festen knoten benachbarte gesmooth werden
    //** aufzeichen -> pruefen -> fixen.
    //** was passiert bei node at interface or boundary
    RVector3 c(0.0, 0.0, 0.0);
    for (std::set< Node * >::iterator it = common.begin(); it != common.end(); it ++){
        c += (*it)->pos();
    }
    this->setPos(c / (double)common.size());
}

} // namespace GIMLI{

