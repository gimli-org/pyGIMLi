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

#include "meshentities.h"

#include "line.h"
#include "mesh.h"
#include "node.h"
#include "shape.h"

#include <map>
#include <algorithm>

namespace GIMLI{

Boundary * findBoundary_(const std::set < Boundary *> & common){
    if (common.size() == 1) {
        return *common.begin();
    } else {
        if (common.size() > 1){
            std::cerr << WHERE_AM_I << " pls. check, this should not happen. "
                    " There is more than one boundary defined." <<
                    common.size() << std::endl;
            std::for_each(common.begin(), common.end(), cerrPtrObject());
        }
    }
    return NULL;
}

Boundary * findBoundary(const Node & n1) {
    if (n1.boundSet().size()) return *n1.boundSet().begin();
    return NULL;;
}

Boundary * findBoundary(const Node & n1, const Node & n2) {
    std::set < Boundary * > common;
    intersectionSet(common, n1.boundSet(), n2.boundSet());
    return findBoundary_(common);
}

Boundary * findBoundary(const Node & n1, const Node & n2, const Node & n3) {
    std::set < Boundary * > common;
    intersectionSet(common, n1.boundSet(), n2.boundSet(), n3.boundSet());
    return findBoundary_(common);
}

Boundary * findBoundary(const Node & n1, const Node & n2, const Node & n3,
                        const Node & n4){
    std::set < Boundary * > common;
    intersectionSet(common, n1.boundSet(), n2.boundSet(),
                    n3.boundSet(), n4.boundSet());
    return findBoundary_(common);
}

Boundary * findBoundary(const std::vector < Node * > & n) {
    if (n.size() == 1) return findBoundary(*n[0]);
    else if (n.size() == 2) return findBoundary(*n[0], *n[1]);
    else if (n.size() == 3) return findBoundary(*n[0], *n[1], *n[2]);
    else if (n.size() == 4) return findBoundary(*n[0], *n[1], *n[2], *n[3]);

    return findBoundary_(findBoundaries(n));
}

std::set < Boundary * > findBoundaries(const std::vector < Node * > & n){
    std::vector < std::set< Boundary * > > bs(n.size());

    for (uint i = 0; i < n.size(); i ++) bs[i] = n[i]->boundSet();

    std::set < Boundary * > common;
    intersectionSet(common, bs);
    return common;
}

Boundary * findCommonBoundary(const Cell & c1, const Cell & c2){
    for (Index i = 0; i < c1.boundaryCount(); i ++){
        Boundary * b = findBoundary(c1.boundaryNodes(i));
        if ((b->leftCell() == &c1 && b->rightCell() == &c2) ||
            (b->leftCell() == &c2 && b->rightCell() == &c1)){
            return b;
        }
    }
    return NULL;
}

Cell * findCommonCell(const std::vector < Node * > & n, bool warn) {
  //** search the cell[s] which is[are] the intersection of all cells in n
    std::vector < std::set< Cell * > > bs;
    for (uint i = 0; i < n.size(); i ++) bs.push_back(n[i]->cellSet());

    std::set < Cell * > common;
    intersectionSet(common, bs);

    if (common.size() == 1) {
        return *common.begin();
    } else {
        if (common.size() > 1){
            if (warn){
                for (uint i = 0; i < n.size(); i ++) std::cout << n[i]->id()<< " " ;
                std::cout <<std::endl;
                std::cerr << WHERE_AM_I << " pls. check, this should not happen. there is more than one cell defined for the given nodes." <<
                common.size() << std::endl;
            }
            return *common.begin();
        }
    }
    return NULL;
}

std::ostream & operator << (std::ostream & str, const MeshEntity & e){
    str << "MeshEntity " << &e << " rtti: " << e.rtti() << " id: " << e.id() << " rtti: " << e.rtti() << "\tN: " ;
    for (uint i = 0; i < e.nodeCount(); i ++) str << e.node(i).id() << " " ;
    return str;
}

std::ostream & operator << (std::ostream & str, const Boundary & e){
    str << "Boundary " << &e << " rtti: " << e.rtti() << " id: " << e.id() << "\tN: " ;
    for (uint i = 0; i < e.nodeCount(); i ++) str << e.node(i).id() << " " ;
    str << " marker: " << e.marker();
    return str;
}

std::ostream & operator << (std::ostream & str, const Edge & e){
    str << "Edge " << &e << " id: " << e.id() << "\t"
      << e.node(0).id() << " " << e.node(1).id()
      << " marker: " << e.marker();
  return str;
}

std::ostream & operator << (std::ostream & str, const TriangleFace & t){
  str << "TriangleFace " << &t << " id: " << t.id() << "\t"
      << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()
      << " attribute: " << t.marker();
  return str;
}

std::ostream & operator << (std::ostream & str, const EdgeCell & c){
  str << "EdgeCell " << &c << " id: " << c.id() << "\tN: ";
  for (uint i = 0; i < c.nodeCount(); i ++) str << c.node(i).id() << " " ;
  str << " attribute: " << c.attribute();
  return str;
}

std::ostream & operator << (std::ostream & str, const Cell & c){
  str << "Cell " << &c << " id: " << c.id() << "\tN: ";
  for (uint i = 0; i < c.nodeCount(); i ++) str << c.node(i).id() << " " ;
  str << " attribute: " << c.attribute();
  return str;
}

std::ostream & operator << (std::ostream & str, const Triangle & t){
  str << "Triangle " << &t << " id: " << t.id() << "\t"
      << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()
      << " attribute: " << t.attribute();
  return str;
}

std::ostream & operator << (std::ostream & str, const Quadrangle & t){
  str << "Quadrangle " << &t << " id: " << t.id() << "\t"
      << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()  << " " << t.node(3).id()
      << " attribute: " << t.attribute();
  return str;
}

std::ostream & operator << (std::ostream & str, const Tetrahedron & t){
    str << "Tetrahedron " << &t << " id: " << t.id() << "\t"
        << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id()  << " " << t.node(3).id()
        << " attribute: " << t.attribute();
  return str;
}

std::ostream & operator << (std::ostream & str, const Hexahedron & t){
    str << "Hexahedron " << &t << " id: " << t.id() << "\t"
        << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id() << " "
        << t.node(3).id() << " " << t.node(4).id() << " " << t.node(5).id() << " "
        << t.node(6).id() << " " << t.node(7).id()
        << " attribute: " << t.attribute();
    return str;
}

std::ostream & operator << (std::ostream & str, const TriPrism & t){
    str << "TrianglePrism" << &t << " id: " << t.id() << "\t"
        << t.node(0).id() << " " << t.node(1).id() << " " << t.node(2).id() << " "
        << t.node(3).id() << " " << t.node(4).id() << " " << t.node(5).id() << " "
        << " attribute: " << t.attribute();
    return str;
}

MeshEntity::MeshEntity()
    : BaseEntity(), shape_(0){
}

MeshEntity::~MeshEntity(){
    for (auto *n: secondaryNodes_) this->deRegisterSecNode_(n);
}

RVector3 MeshEntity::center() const {
    if (!shape_){
        log(Error, "no shape defined");
        return 0;
    }
    return shape_->center();
}

double MeshEntity::size() const {
    if (!shape_){
        log(Error, "no shape defined");
        return 0;
    }
    return shape_->domainSize();
}

void MeshEntity::fillShape_(){
    if (shape_){
        shape_->setNodesPtr(this->nodeVector_);
        //* create Shape function and cache them to avoid multithreading problems ...
        shape_->changed();
        ShapeFunctionCache::instance().shapeFunctions(*shape_);
    }
}

void MeshEntity::setNodes(const std::vector < Node * > & nodes){
    if (nodes.size() > 0){
        deRegisterNodes_();
        if (nodeVector_.size() != nodes.size()) {
            nodeVector_.resize(nodes.size());
        }
        std::copy(nodes.begin(), nodes.end(), &nodeVector_[0]);
        registerNodes_();
        fillShape_();
    } else {
        std::cerr << WHERE_AM_I << " not enough nodes to fill meshEntity " << std::endl;
    }
}

void MeshEntity::registerSecNode_(Node *n){
}
void MeshEntity::deRegisterSecNode_(Node *n){

}

void MeshEntity::addSecondaryNode(Node * n) {
    secondaryNodes_.push_back(n);
    registerSecNode_(n);
}

void MeshEntity::delSecondaryNode(Node * n) {
    deRegisterSecNode_(n);
    secondaryNodes_.erase(std::remove(secondaryNodes_.begin(),
                                      secondaryNodes_.end(), n),
                          secondaryNodes_.end());
}
const std::vector < Node * > & MeshEntity::secondaryNodes() const {
    return secondaryNodes_;
}

const std::vector < Node * > MeshEntity::allNodes() const {
    std::vector < Node * > ns;
    for (Index i = 0; i < this->nodeVector_.size(); i++ ){
        ns.push_back(this->nodeVector_[i]);
    }
    for (Index i = 0; i < secondaryNodes_.size(); i++ ){
        ns.push_back(secondaryNodes_[i]);
    }
    return ns;
}

Index MeshEntity::allNodeCount() const{
    return nodeCount() + secondaryNodes_.size();
}

void MeshEntity::registerNodes_(){
}

void MeshEntity::deRegisterNodes_(){
}

IndexArray MeshEntity::ids() const {
    IndexArray idVec(nodeCount());

    for (uint i = 0; i < nodeCount(); i ++) {
        idVec[i] = node(i).id();
    }
    return idVec;
}

RVector3 MeshEntity::rst(uint i) const {
    return shape_->rst(node(i).pos()).round(TOLERANCE);
}

std::vector < PolynomialFunction < double > > MeshEntity::createShapeFunctions() const{
    std::cerr << "need shape function implementation for meshEntity " << rtti() << std::endl;
    THROW_TO_IMPL
    return std::vector < PolynomialFunction < double > > ();
}

RVector MeshEntity::N(const RVector3 & rst) const {
    RVector n(nodeCount(), 0.0);
    this->N(rst, n);
    return n;
}

void MeshEntity::N(const RVector3 & rst, RVector & n) const {
    const std::vector< PolynomialFunction < double > > &N = ShapeFunctionCache::instance().shapeFunctions(*this);

    for (Index i = 0; i < N.size(); i ++) {
        n[i] = N[i](rst);
    }
}

RVector MeshEntity::dNdL(const RVector3 & rst, uint i) const {

    const std::vector< PolynomialFunction < double > > &dNL =
        ShapeFunctionCache::instance().deriveShapeFunctions(*this, i);

    RVector ret(dNL.size());

    for (Index i = 0; i < dNL.size(); i ++) {
         ret[i] = dNL[i](rst);
    }

    return ret;
}

RMatrix MeshEntity::dNdL(const RVector3 & rst) const {
    RMatrix ret;
    ret.push_back(this->dNdL(rst, 0));
    ret.push_back(this->dNdL(rst, 1));
    ret.push_back(this->dNdL(rst, 2));
    return ret;
}

double MeshEntity::pot(const RVector3 & xyz, const RVector & u) const {
    return sum(u(this->ids()) * this->N(shape().rst(xyz)));
}

RVector3 MeshEntity::vec(const RVector3 & xyz,
                         const R3Vector & v) const{
    return RVector3(pot(xyz, x(v)), pot(xyz, y(v)), pot(xyz, z(v)));
}

RVector3 MeshEntity::grad(const RVector3 & xyz, const RVector & u) const {

    RVector3 rst(shape_->rst(xyz));

    RMatrix MdNdL;
    MdNdL.push_back(dNdL(rst, 0));
    MdNdL.push_back(dNdL(rst, 1));
    MdNdL.push_back(dNdL(rst, 2));

    RVector up(u(this->ids()));
    RVector3 gr;
    gr[0] = sum(up * MdNdL.transMult(shape_->invJacobian().col(0)));
    gr[1] = sum(up * MdNdL.transMult(shape_->invJacobian().col(1)));
    gr[2] = sum(up * MdNdL.transMult(shape_->invJacobian().col(2)));
    return gr;
}

void MeshEntity::changed(){
    this->shape_->changed();
    uCache_.setValid(false);
    gradUCache_.setValid(false);
}

bool MeshEntity::enforcePositiveDirection(){
    // __MS(this->shape_->createJacobian().det())
    // if (this->shape_->createJacobian().det() < 0){
    //     std::reverse(nodeVector_.begin(), nodeVector_.end());
    //     this->changed();

        // __MS(this->shape_->createJacobian().det())
    //     return true;
    // }
    __M
    return false;
}

//############### Cell ##################

Cell::Cell() : MeshEntity(), attribute_(){
}

Cell::Cell(const std::vector < Node * > & nodes)
    : MeshEntity(), attribute_(0.0){
    this->setNodes(nodes);
}

Cell::~Cell(){
    for (auto *n: secondaryNodes_) this->deRegisterSecNode_(n);
    deRegisterNodes_();
}

void Cell::registerSecNode_(Node *n){
    n->insertCell(this);
}

void Cell::deRegisterSecNode_(Node *n){
    n->eraseCell(this);
}


Node * Cell::oppositeTo(const Boundary & bound){
    THROW_TO_IMPL
    //** maybee obsolete pls check

//     std::set < Node * > nodes;
//     for (int i = 0; i < 3; i++) nodes.insert(&node(i));
//
//     for (uint i = 0; i < bound.nodeCount(); i ++) nodes.erase(& bound.node(i));
//
//     if (nodes.size() == 1){
//         return (*nodes.begin());
//     } else {
//         std::cerr << bound << std::endl << *this << std::endl;
//         std::cerr << WHERE_AM_I << " this triangle does not contain the edge" << std::endl;
//     }
     return NULL;
}

void Cell::cleanNeighborInfos(){
    for (uint i = 0; i < this->neighborCellCount(); i ++){
        neighborCells_[i] = NULL;
    }
}

Boundary * Cell::boundaryTo(const RVector & sf){

    double maxSF = max(sf);
    double minSF = min(sf);

    IndexArray maxIdx(find(sf == maxSF));
    IndexArray minIdx(find(sf == minSF));

    std::set < Boundary * > common;

    if (maxIdx.size() > 1){
        std::vector < std::set< Boundary * > > bs;
        for (Index i = 0; i < maxIdx.size(); i ++) {
            bs.push_back(node(maxIdx[i]).boundSet());
        }

        intersectionSet(common, bs);
    } else {
        common = node(maxIdx[0]).boundSet();
    }

    if (common.size() == 0) {
        // this happens when the cell has no boundaries.
        return NULL;
    }
    if (common.size() == 1) return (*common.begin());

    // remove bounds from foreign cells
//     std::cout << WHERE << "common.size() remove foreign: " << common.size()<< std::endl;
    for (std::set < Boundary * >::iterator it = common.begin(); it != common.end();){

        if ((*it)->leftCell() != this && (*it)->rightCell() != this){
//             std::cout << WHERE << *(*it)<< std::endl;
//             std::cout << WHERE << this << " " << (*it)->leftCell() << " " << (*it)->rightCell() << std::endl;
            common.erase(it++);
        } else {
            ++it;
        }
    }

//     std::cout << WHERE << "common.size() remove foreign: " << common.size()<< std::endl;

    std::set < Boundary * > commonSub;

    if (minIdx.size() > 1){
        std::vector < std::set< Boundary * > > bs;
        for (Index i = 0; i < minIdx.size(); i ++) bs.push_back(node(minIdx[i]).boundSet());

        intersectionSet(commonSub, bs);
    } else {
        commonSub = node(minIdx[0]).boundSet();
    }

    for (std::set < Boundary * >::iterator it = commonSub.begin();
         it != commonSub.end(); it++){
        common.erase(*(it));
    }

//     std::cout << WHERE << "common.size() remove foreign: "
//      << common.size()<< std::endl;

    if (common.size() == 0) {
        std::cerr << " this.should not happen" << std::endl;
        std::cout << rtti() << " " << *this  << std::endl;
        std::cout << sf  << std::endl;
        THROW_TO_IMPL
        return NULL;
    }

//     std::cout << "common.size() remove sub: " << common.size()<< std::endl;

    return (*common.begin());
}

Boundary * Cell::boundary(Index i){
    ASSERT_RANGE(i, 0, boundaryCount())
    return findBoundary(boundaryNodes(i));
}

Cell * Cell::neighborCell(const RVector & sf){
    if (haveInfNaN(sf)){
        __MS("fixme " << sf)
        return NULL;
    }

    //** hack for edge, triangle and tetrahedron. pls refactor
    if (((sf.size() == 2) && (shape_->dim() == 1)) ||
        ((sf.size() == 3) && (shape_->dim() == 2))
        //|| ((sf.size() == 4) && (shape_->dim() == 3))
    ){
        Index minId = find(sf == min(sf))[0];

//         Boundary * b = boundaryTo(sf);
//        __MS(b)
//        __MS(b->leftCell())
//        __MS(b->rightCell())

        return neighborCells_[minId];
    }

    Boundary * b = boundaryTo(sf);

//     std::cout << *b << " " << b->leftCell() << " "<< b->rightCell() << std::endl;
//     this->setMarker(-33);

    if (b){
        if (b->rightCell() == this) return b->leftCell();
        if (b->leftCell() == this) return b->rightCell();
    }

    return NULL;
}

void Cell::findNeighborCell(uint i){
    if (!neighborCells_[i]){
        std::vector < Node * > n(boundaryNodes(i));

        std::set < Cell *> common;
        std::set < Cell *> commonTmp;

        if (n.size() > 1) {
            intersectionSet(common, n[0]->cellSet(), n[1]->cellSet());

            for (uint j = 2; j < n.size(); j ++){
                commonTmp = common;
                intersectionSet(common, commonTmp, n[j]->cellSet());
            }
        } else {
            common = n[0]->cellSet();
        }

        common.erase(this);
        if (common.size() == 1) neighborCells_[i] = *common.begin(); else neighborCells_[i] = NULL;
    }
}

void Cell::registerNodes_(){
    for (auto n: nodeVector_) n->insertCell(this);
}

void Cell::deRegisterNodes_(){
    for (auto n: nodeVector_) n->eraseCell(this);
}

Boundary::Boundary()
    : MeshEntity(), leftCell_(NULL), rightCell_(NULL) {
}

Boundary::Boundary(const std::vector < Node * > & nodes)
    : MeshEntity(), leftCell_(NULL), rightCell_(NULL) {
    this->setNodes(nodes);
}

Boundary::~Boundary(){
    deRegisterNodes_();
}

void Boundary::registerNodes_(){
    for (auto n: nodeVector_) n->insertBoundary(this);
}
void Boundary::deRegisterNodes_(){
    for (auto n: nodeVector_) n->eraseBoundary(this);
}
void Boundary::registerSecNode_(Node *n){
    n->insertBoundary(this);
}
void Boundary::deRegisterSecNode_(Node *n){
    n->eraseBoundary(this);
}


RVector3 Boundary::rst(uint i) const {
    if (this->nodeCount() == shape_->nodeCount()) return shape_->rst(i);
    std::cerr << "need local coordinate function implementation for meshEntity " << rtti() << std::endl;
    THROW_TO_IMPL
    return RVector3(0.0, 0.0, 0.0);
}

RVector3 Boundary::norm() const {
    return shape_->norm();
}

RVector3 Boundary::norm(const Cell & c) const {
    if (!this->normShowsOutside(c)) return -this->norm();
    return this->norm();
}

bool Boundary::normShowsOutside(const Cell & cell) const {
    RVector3 n(this->norm());
    RVector3 bc(this->center());
    RVector3 cc(cell.center());
    return (cc-(bc+n)).abs() > (cc-(bc-n)).abs();
}

void Boundary::swapNorm(bool withNeighbours){
    std::reverse(nodeVector_.begin(), nodeVector_.end());
    if (withNeighbours){
        Cell * left = this->leftCell();
        this->setLeftCell(this->rightCell());
        this->setRightCell(left);
    }
    this->changed();
}

NodeBoundary::NodeBoundary(const std::vector < Node * > & nodes)
    : Boundary(nodes){
    shape_ = new NodeShape(this);
}

NodeBoundary::NodeBoundary(Node & n1)
    : Boundary() {
    shape_ = new NodeShape(this);
    setNodes(n1);
}

NodeBoundary::~NodeBoundary(){
    delete shape_;
}

void NodeBoundary::setNodes(Node & n1){
    const std::vector < Node * > nodes{&n1};
    MeshEntity::setNodes(nodes);
}

RVector3 NodeBoundary::norm() const{
    const Cell *c = &this->leftCell();
    if (c != 0){
        return (this->center() - c->center()).norm();
    } else {
        return RVector3(1.0, 0.0, 0.0);
    }
}

Edge::Edge(const std::vector < Node * > & nodes)
    : Boundary(nodes){
    shape_ = new EdgeShape(this);
}

Edge::Edge(Node & n1, Node & n2){
    shape_ = new EdgeShape(this);
    setNodes(n1, n2);
}

Edge::~Edge(){
    delete shape_;
}

void Edge::setNodes(Node & n1, Node & n2){
    if ((&n1 == &n2)){
        throwError(WHERE + " Edge nodes not valid " +
                   str(n1) + " " + str(n2) );
    }
    const std::vector < Node * > nodes{&n1, &n2};
    MeshEntity::setNodes(nodes);
}

std::vector < PolynomialFunction < double > > Edge::createShapeFunctions() const {
     return createPolynomialShapeFunctions(*this, 2, true, false);
}

int Edge::swap(){
    if ((marker_ != 0) || (leftCell_ == NULL) || (rightCell_ == NULL)) {
    //    cout<< WHERE_AM_I << " " << marker_ << " " << leftCell_  << " " << rightCell_ << std::endl;
        return 0;
    } else if ((leftCell_->rtti() == MESH_TRIANGLE_RTTI) && (rightCell_->rtti() == MESH_TRIANGLE_RTTI)) {

        Node * oA = &node(0);
        Node * oB = &node(1);
        Triangle *left = dynamic_cast< Triangle * >(leftCell_);
        Triangle *right = dynamic_cast< Triangle * >(rightCell_);

        Node * oL = left->oppositeTo(* this);
        Node * oR = right->oppositeTo(* this);

        if (oL == NULL || oR == NULL){
            std::cout << *this << std::endl
                << left << std::endl
                << right << std::endl;
        if (oL != NULL) std::cout << "oL " << oL->id() << std::endl;
        if (oR != NULL) std::cout << "oR " << oR->id() << std::endl;
        throwError(WHERE);
    }

    //** swap only when the resulting triangles have the same sign
    //** abort swapping for concav domain, possible check by angles
    if (sign(jacobianDetXY(oL->pos(), oR->pos(), oB->pos())) !=
	       sign(jacobianDetXY(oL->pos(), oA->pos(), oR->pos()))){
        return 0;
    }

    right->setNodes(* oL, * oA, * oR);

    setNodes(* oL, * oR);
    //    cout << "\t" << *this << std::endl;

    if (leftCell_ == rightCell_){
        std::cerr << WHERE << " Edge " << id() << " wrong swapped " << std::endl;
        std::cerr << "LeftElement: " << left->id()
        << "; RightElement: " << right->id() << std::endl;
        std::cerr << "NodeA: " << oA->id() << ", NodeB: " << oB->id()
        << ", NodeL: " << oL->id() << ", NodeR: " << oR->id() << std::endl;
        return 0;
    }

    left->setNodes(* oL, * oR, * oB);
    right->setNodes(* oL, * oA, * oR);

    return 1;
  }
  return 0;
}

// RVector Edge::N(const RVector3 & L) const {
//     return shape_->N(L);
// }

Edge3::Edge3(const std::vector < Node * > & nodes) : Edge(nodes){
}

Edge3::~Edge3(){
}

RVector3 Edge3::rst(uint i) const {
    if (i == 2) return RVector3(0.5, 0.0, 0.0);
    return shape_->rst(i);
}

std::vector < PolynomialFunction < double > > Edge3::createShapeFunctions() const {
    return createPolynomialShapeFunctions(*this, 3, true, false);
}

TriangleFace::TriangleFace(const std::vector < Node * > & nodes)
    : Boundary(nodes){
    shape_ = new TriangleShape(this);
}

TriangleFace::TriangleFace(Node & n1, Node & n2, Node & n3){
    shape_ = new TriangleShape(this);
    setNodes(n1, n2, n3);
}

TriangleFace::~TriangleFace(){
    delete shape_;
}

std::vector < PolynomialFunction < double > > TriangleFace::createShapeFunctions() const {
    return createPolynomialShapeFunctions(*this, 2, true, false);
}

void TriangleFace::setNodes(Node & n1, Node & n2, Node & n3){
    if ((&n1 == &n2) || (&n1 == &n3) || (&n2 == &n3)){
        std::cerr << WHERE << " TriangleFace nodes not valid " << n1 << " " << n2 << " " <<  n3 << std::endl;
        throwError(WHERE);
    }
    const std::vector < Node * > nodes{&n1, &n2, &n3};
    MeshEntity::setNodes(nodes);
}

Triangle6Face::Triangle6Face(const std::vector < Node * > & nodes) : TriangleFace(nodes){
}

Triangle6Face::~Triangle6Face(){
}

std::vector < PolynomialFunction < double > > Triangle6Face::createShapeFunctions() const {
    return createPolynomialShapeFunctions(*this, 3, true, false);
}

RVector3 Triangle6Face::rst(uint i) const {
    if (i == 3) return RVector3((shape_->rst(0) + shape_->rst(1)) / 2);
    if (i == 4) return RVector3((shape_->rst(1) + shape_->rst(2)) / 2);
    if (i == 5) return RVector3((shape_->rst(2) + shape_->rst(0)) / 2);
    return shape_->rst(i);
}


QuadrangleFace::QuadrangleFace(const std::vector < Node * > & nodes)
    : Boundary(nodes){
    shape_ = new QuadrangleShape(this);
}

QuadrangleFace::QuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4){
    shape_ = new QuadrangleShape(this);
    setNodes(n1, n2, n3, n4);
}

QuadrangleFace::~QuadrangleFace(){
    delete shape_;
}

void QuadrangleFace::setNodes(Node & n1, Node & n2, Node & n3, Node & n4){
    if ((&n1 == &n2) || (&n1 == &n3) || (&n2 == &n3)){
        std::cerr << WHERE << " QuadrangleFace nodes not valid " << n1 << " "
                << n2 << " " <<  n3 << " " << n4<< std::endl;
        throwError(WHERE);
    }
    const std::vector < Node * > nodes{&n1, &n2, &n3, &n4};
    MeshEntity::setNodes(nodes);
}

std::vector < PolynomialFunction < double > > QuadrangleFace::createShapeFunctions() const {
    return createPolynomialShapeFunctions(*this, 2, true, true);
}

Quadrangle8Face::Quadrangle8Face(const std::vector < Node * > & nodes) : QuadrangleFace(nodes){
}

Quadrangle8Face::~Quadrangle8Face(){
}

RVector3 Quadrangle8Face::rst(uint i) const {
    if (i == 4) return RVector3((shape_->rst(0) + shape_->rst(1)) / 2.);
    if (i == 5) return RVector3((shape_->rst(1) + shape_->rst(2)) / 2.);
    if (i == 6) return RVector3((shape_->rst(2) + shape_->rst(3)) / 2.);
    if (i == 7) return RVector3((shape_->rst(3) + shape_->rst(0)) / 2.);
    return shape_->rst(i);
}

std::vector < PolynomialFunction < double > > Quadrangle8Face::createShapeFunctions() const {
    return createPolynomialShapeFunctions(*this, 3, true, true);
}

PolygonFace::PolygonFace(const std::vector < Node * > & nodes)
    : Boundary(nodes){
    shape_ = new PolygonShape(this);
}

PolygonFace::~PolygonFace(){
}

void PolygonFace::insertNode(Node * n, double tol){
    for (Index i = 0; i < nodeCount(); i ++){
        if (n->id() == this->node(i).id()) return;

        if (this->node(i).pos().distance(n->pos()) < tol){
            __MS(*this)
            __MS(n)
            log(Error, "PolygonFace::insertNode. Duplicate node position found. "
                       "Node need to touch the Polygon face or its edge but not the corner nodes.");
        }

        Line segment(this->node(i).pos(),
                     this->node((i+1)%this->nodeCount()).pos());
        int pFkt;
        if (segment.touch1(n->pos(), pFkt, tol)){
            if (pFkt == 3){
                // __MS("insert edge")
                std::vector < Node * > nodes;
                for (Index j = 0; j < i + 1; j ++){
                    nodes.push_back(&this->node(j));
                }
                nodes.push_back(n);
                for (Index j = i + 1; j < this->nodeCount(); j ++){
                    nodes.push_back(&this->node(j));
                }
                n->setState(Connected);
                MeshEntity::setNodes(nodes);
                shape_->resizeNodeSize_(this->nodeCount());
                return;
            }
        }
    }
    this->addSecondaryNode(n);
    n->setState(Secondary);
    n->setSecondaryParent(this);
}

void PolygonFace::addSubface(const std::vector < Node * > & nodes, bool isHole){
    this->subfaces_.push_back(nodes);
    if (!isHole) for (auto *n: nodes) n->insertBoundary(this);
}
const std::vector < Node * >  & PolygonFace::subface(Index i) const {
    return this->subfaces_[i];
}
void PolygonFace::addHoleMarker(const RVector3 & pos){
    holeMarker_.push_back(pos);
}
void PolygonFace::delHoleMarker(const RVector3 & pos){
    holeMarker_.erase(std::remove(holeMarker_.begin(),
                                  holeMarker_.end(), pos),
                      holeMarker_.end());
}
const PolygonFace::HoleMarkerList & PolygonFace::holeMarkers() const {
    return holeMarker_;
}

EdgeCell::EdgeCell(const std::vector < Node * > & nodes)
    : Cell(nodes){
    shape_ = new EdgeShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

EdgeCell::EdgeCell(Node & n1, Node & n2)
    : Cell() {
    shape_ = new EdgeShape(this);
    setNodes(n1, n2);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

EdgeCell::~EdgeCell(){
    delete shape_;
}

void EdgeCell::setNodes(Node & n1, Node & n2){
    const std::vector < Node * > nodes{&n1, & n2};
    MeshEntity::setNodes(nodes);
}

std::vector < Node * > EdgeCell::boundaryNodes(Index i) const {
    std::vector < Node * > nodes(1);
    nodes[0] = nodeVector_[(i+1)%2];
    return nodes;
}

// void EdgeCell::findneighborCell(uint id){
//     //     case: 0 - 1
//     //     case: 1 - 0
//
//     std::set < Cell * > common(nodeVector_[(id + 1)%2]->cellSet());
//     common.erase(this);
//     if (common.size() == 1) neighborCells_[id] = *common.begin(); else neighborCells_[id] = NULL;
// }

std::vector < PolynomialFunction < double > > EdgeCell::createShapeFunctions() const{
//     __M
    return createPolynomialShapeFunctions(*this, 2, true, false);
}

// RVector EdgeCell::N(const RVector3 & L) const{
//     return shape_->N(L);
// }
//
// RVector EdgeCell::dNdL(const RVector3 & coord, uint dim) const {
//     RVector dN(nodeCount());
//     dN[0] = - 1.0;
//     dN[1] =   1.0;
//     return dN;
// }

Edge3Cell::Edge3Cell(const std::vector < Node * > & nodes) : EdgeCell(nodes){
}

Edge3Cell::~Edge3Cell(){
}

RVector3 Edge3Cell::rst(uint i) const {
    if (i == 2) return RVector3(0.5, 0.0, 0.0);
    return shape_->rst(i);
}

std::vector < PolynomialFunction < double > > Edge3Cell::createShapeFunctions() const{
    return createPolynomialShapeFunctions(*this, 3, true, false);
}

Triangle::Triangle(const std::vector < Node * > & nodes)
    : Cell(nodes){
    shape_ = new TriangleShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Triangle::Triangle(Node & n1, Node & n2, Node & n3)
    : Cell(){
    shape_ = new TriangleShape(this);
    setNodes(n1, n2, n3);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Triangle::~Triangle(){
    delete shape_;
}

void Triangle::setNodes(Node & n1, Node & n2, Node & n3){
    if ((&n1 == &n2) || (&n1 == &n3) || (&n2 == &n3)){
        std::cerr << WHERE << " Triangle nodes not valid " << n1 << " " << n2 << " " << n3 << std::endl;
        throwError(WHERE);
    }
    const std::vector < Node * > nodes{&n1, &n2, &n3};
    MeshEntity::setNodes(nodes);
}


std::vector < PolynomialFunction < double > > Triangle::createShapeFunctions() const{
    return createPolynomialShapeFunctions(*this, 2, true, false);
}

std::vector < Node * > Triangle::boundaryNodes(Index i) const{
    // 0 -> 1..2
    // 1 -> 2..0
    // 2 -> 0..1
    std::vector < Node * > nodes(2);
    nodes[0] = nodeVector_[(i+1)%3];
    nodes[1] = nodeVector_[(i+2)%3];
    return nodes;
}

Triangle6::Triangle6(const std::vector < Node * > & nodes) : Triangle(nodes){
}

Triangle6::~Triangle6(){
}

std::vector < PolynomialFunction < double > > Triangle6::createShapeFunctions() const{
    return createPolynomialShapeFunctions(*this, 3, true, false);
}

Quadrangle::Quadrangle(const std::vector < Node * > & nodes) : Cell(nodes){
    shape_ = new QuadrangleShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Quadrangle::Quadrangle(Node & n1, Node & n2, Node & n3, Node & n4): Cell(){
    shape_ = new QuadrangleShape(this);
    setNodes(n1, n2, n3, n4);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Quadrangle::~Quadrangle(){
    delete shape_;
}

void Quadrangle::setNodes(Node & n1, Node & n2, Node & n3, Node & n4){
    const std::vector < Node * > nodes{&n1, &n2, &n3, &n4};
    MeshEntity::setNodes(nodes);
}

// void Quadrangle::findneighborCell(uint id){
//     std::set < Cell * > common;
//
//     intersectionSet(common, nodeVector_[(id + 1)%4]->cellSet(), nodeVector_[(id + 2)%4]->cellSet());
//     common.erase(this);
//     if (common.size() == 1) neighborCells_[id] = *common.begin(); else neighborCells_[id] = NULL;

    //** start check;
//     if (this->id() ==){
//
//         int pIdx = 0;
//         shape_->touch1(neighborCells_[id]->center(), true, pIdx);
//         if (pIdx != id){
//             std::cout << "something goes wrong here idx:" << id << " pIdx = " << pIdx << " "
//                     <<  id - pIdx<< std::endl;
//             std::cout << "origin:" << *this << std::endl;
//             std::cout << "test:" << *neighborCells_[id] << std::endl;
//             for (uint i = 0; i <4 ; i ++){
//                 std::cout << this->node(i)<< std::endl;
//             }
//             for (uint i = 0; i <4 ; i ++){
//                 common.clear();
//                 intersectionSet(common, nodeVector_[(i)]->cellSet(),
//                                         nodeVector_[(i + 1)%4]->cellSet());
//                 common.erase(this);
//                 if (common.size() == 1) std::cout << i << " "<< **common.begin() << std::endl;
//             }
//
//             exit(1);
//         }
//     }
//}

std::vector < PolynomialFunction < double > > Quadrangle::createShapeFunctions() const{
    return createPolynomialShapeFunctions(*this, 2, true, true);
}

std::vector < Node * > Quadrangle::boundaryNodes(Index i) const {

    std::vector < Node * > nodes(2);
    for (Index j = 0; j < 2; j ++) nodes[j] = nodeVector_[(i+j)%4];
    return nodes;
}

Quadrangle8::Quadrangle8(const std::vector < Node * > & nodes) : Quadrangle(nodes){
}

Quadrangle8::~Quadrangle8(){
}

std::vector < PolynomialFunction < double > > Quadrangle8::createShapeFunctions() const{
//     return createPolynomialShapeFunctions(*this, 3, false, false);
    return createPolynomialShapeFunctions(*this, 3, true, true);
}

Tetrahedron::Tetrahedron(const std::vector < Node * > & nodes) : Cell(nodes){
    shape_ = new TetrahedronShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Tetrahedron::Tetrahedron(Node & n1, Node & n2, Node & n3, Node & n4) : Cell() {
    shape_ = new TetrahedronShape(this);
    setNodes(n1, n2, n3, n4);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Tetrahedron::~Tetrahedron(){
    delete shape_;
}

void Tetrahedron::setNodes(Node & n1, Node & n2, Node & n3, Node & n4){
    const std::vector < Node * > nodes{&n1, &n2, &n3, &n4};
    MeshEntity::setNodes(nodes);
}

// void Tetrahedron::findneighborCell(uint i){
//     std::set < Cell * > c1, common;
//     intersectionSet(c1, nodeVector_[(i + 1)%4]->cellSet(), nodeVector_[(i + 2)%4]->cellSet());
//     intersectionSet(common, c1, nodeVector_[(i + 3)%4]->cellSet());
//
//     common.erase(this);
//     if (common.size() == 1) neighborCells_[i] = *common.begin(); else neighborCells_[i] = NULL;
// }

std::vector < PolynomialFunction < double > > Tetrahedron::createShapeFunctions() const{
    return createPolynomialShapeFunctions(*this, 2, true, false);
}

std::vector < Node * > Tetrahedron::boundaryNodes(Index i) const {
    // 0 -> 1..2..3
    // 1 -> 2..0..3
    // 2 -> 0..1..3
    // 3 -> 0..2..1

    std::vector < Node * > nodes(3);
    for (Index j = 0; j < 3; j ++) {
        nodes[j] = nodeVector_[TetrahedronFacesID[i][j]];
//         std::cout << *nodes[j] << std::endl;
    }
    return nodes;
}

Tetrahedron10::Tetrahedron10(const std::vector < Node * > & nodes) : Tetrahedron(nodes){
}

Tetrahedron10::~Tetrahedron10(){
}

std::vector < PolynomialFunction < double > > Tetrahedron10::createShapeFunctions() const{
    return createPolynomialShapeFunctions(*this, 3, true, false);
}

Hexahedron::Hexahedron(const std::vector < Node * > & nodes)
    : Cell(nodes) {
    shape_ = new HexahedronShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Hexahedron::~Hexahedron(){
    delete shape_;
}

std::vector < PolynomialFunction < double > > Hexahedron::createShapeFunctions() const{
//     std::cout << "return createPolynomialShapeFunctions(*this, 2, false, false);" << std::endl;
//     return createPolynomialShapeFunctions(*this, 2, false, false);
    return createPolynomialShapeFunctions(*this, 2, true, true);
}

// void Hexahedron::findneighborCell(uint i){
//     std::vector < Node * > n(boundaryNodes(i));
//
//     std::set < Cell *> common;
//     std::set < Cell *> commonTmp;
//
//     if (n.size() > 1) intersectionSet(common, n[0]->cellSet(), n[1]->cellSet());
//
//     for (uint j = 2; j < n.size(); j ++){
//         commonTmp = common;
//         intersectionSet(common, commonTmp, n[j]->cellSet());
//     }
//
//     common.erase(this);
//     if (common.size() == 1) neighborCells_[i] = *common.begin(); else neighborCells_[i] = NULL;
// }

std::vector < Node * > Hexahedron::boundaryNodes(Index i) const {
    std::vector < Node * > nodes(4);
    nodes[0] = nodeVector_[HexahedronFacesID[i][0]];
    nodes[1] = nodeVector_[HexahedronFacesID[i][1]];
    nodes[2] = nodeVector_[HexahedronFacesID[i][2]];
    nodes[3] = nodeVector_[HexahedronFacesID[i][3]];
    return nodes;
}

Hexahedron20::Hexahedron20(const std::vector < Node * > & nodes) : Hexahedron(nodes) {
}

Hexahedron20::~Hexahedron20(){

}

std::vector < PolynomialFunction < double > > Hexahedron20::createShapeFunctions() const{
//     return createPolynomialShapeFunctions(*this, 3, false, false);
    return createPolynomialShapeFunctions(*this, 3, true, true);
}

std::vector < Node * > Hexahedron20::boundaryNodes(Index i) const {
    std::vector < Node * > nodes(8);
    for (Index j = 0; j < 8; j ++) nodes[j] = nodeVector_[Hexahedron20FacesID[i][j]];
    return nodes;
}

TriPrism::TriPrism(const std::vector < Node * > & nodes)
    : Cell(nodes){
    shape_ = new TriPrismShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

TriPrism::~TriPrism(){
    delete shape_;
}

std::vector < PolynomialFunction < double > > TriPrism::createShapeFunctions() const{
    //return createPolynomialShapeFunctions(*this, 2, true, true);
    RVector e2(2); e2[0] = 0; e2[1] =  1; // x

    RPolynomialFunction T3_2(e2, RVector(0));
    RPolynomialFunction T3_3(RVector(0), e2);
    RPolynomialFunction T3_1 = -(-1.0 + T3_2 + T3_3);

    RPolynomialFunction E2_2T(RVector(0), RVector(0), e2);
    RPolynomialFunction E2_1T = -(-1.0 + E2_2T);
    std::vector < PolynomialFunction < double > > ret;

    ret.push_back(T3_1 * E2_1T);
    ret.push_back(T3_2 * E2_1T);
    ret.push_back(T3_3 * E2_1T);
    ret.push_back(T3_1 * E2_2T);
    ret.push_back(T3_2 * E2_2T);
    ret.push_back(T3_3 * E2_2T);
    return ret;

}

std::vector < Node * > TriPrism::boundaryNodes(Index i) const {
    std::vector < Node * > nodes;
    for (uint j = 0; j < 3; j ++){
        nodes.push_back(nodeVector_[TriPrismFacesID[i][j]]);
    }
    if (TriPrismFacesID[i][3] != 255) nodes.push_back(nodeVector_[TriPrismFacesID[i][3]]);

    return nodes;
}

TriPrism15::TriPrism15(const std::vector < Node * > & nodes) : TriPrism(nodes){
}

TriPrism15::~TriPrism15(){
}

std::vector < PolynomialFunction < double > > TriPrism15::createShapeFunctions() const{
    //return createPolynomialShapeFunctions(*this, 3, false, false);

    RVector xy(9, 0.0);
    xy[0] = 1.;  // 1
    xy[1] = 1.;  // x
    xy[2] = 1.;  // x^2
    xy[3] = 1.;  // y
    xy[4] = 1.;  // yx
    xy[5] = 0.;  // yx^2
    xy[6] = 1.;  // y^2
    xy[7] = 0.;  // y^2x
    xy[8] = 0.;  // y^2x^2

    RVector start(27);
    start.setVal(xy, 0, 9);
    start.setVal(xy, 9, 18);
    start.setVal(xy, 18, 27);
    start[18 + 2] = 0.; // z^2x^2
    start[18 + 4] = 0.; // z^2xy
    start[18 + 6] = 0.; // z^2y^2

    return createPolynomialShapeFunctions(*this, 3, false, false, start);


    //#return createPolynomialShapeFunctions(*this, 3, true, true);
}

Pyramid::Pyramid(const std::vector < Node * > & nodes)
    : Cell(nodes){
    shape_ = new PyramidShape(this);
    neighborCells_.resize(this->neighborCellCount(), NULL);
}

Pyramid::~Pyramid(){
}

std::vector < PolynomialFunction < double > > Pyramid::createShapeFunctions() const{
    THROW_TO_IMPL
    return std::vector < PolynomialFunction < double > >();
}

std::vector < Node * > Pyramid::boundaryNodes(Index i) const {
    THROW_TO_IMPL
    return std::vector < Node * > ();
}

Pyramid13::Pyramid13(const std::vector < Node * > & nodes): Pyramid(nodes){
}

Pyramid13::~Pyramid13(){
}

std::vector < PolynomialFunction < double > > Pyramid13::createShapeFunctions() const{
    THROW_TO_IMPL
    return std::vector < PolynomialFunction < double > >();
}

std::ostream & operator << (std::ostream & str, const std::set < GIMLI::MeshEntity * > & ents){
    for (auto *e: ents){
        str << e->id() << " ";
    }
    return str;
}

DLLEXPORT std::ostream & operator << (std::ostream & str, const std::set < GIMLI::Boundary * > & bounds){
    for (auto *b: bounds){
        str << b->id() << " ";
    }
    return str;
}


} // namespace GIMLI{
