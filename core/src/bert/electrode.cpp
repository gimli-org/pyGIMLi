/******************************************************************************
 *   Copyright (C) 2005-2024 by the resistivity.net development team          *
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

#include "electrode.h"

#include <elementmatrix.h>
#include <meshentities.h>
#include <node.h>
#include <numericbase.h>
#include <shape.h>

namespace GIMLI{

// class Electrode -> ElectrodeNode -> ElectrodeRVector3 -> ElectrodeBoundary -> ElectrodeDomain

Electrode::Electrode(){
    this->setValid(false);
    pos_.assign(0.0, 0.0, 0.0);
}

Electrode::Electrode(const RVector3 & pos, int id) : pos_(pos){
    this->setId(id);
    this->setValid(true);
}

Electrode::Electrode(double x, double y, double z){
    this->setValid(true);
    pos_.assign(x, y, z);
}

Electrode::~Electrode(){
}

Electrode::Electrode(const Electrode & el){
  pos_ = el.pos();
  this->setId(el.id());
  this->setValid(el.valid());
}

Electrode & Electrode::operator=(const Electrode & el){
  if (this != & el){
    pos_ = el.pos();
    this->setId(el.id());
    this->setValid(el.valid());
  }
  return *this;
}

ElectrodeShape::ElectrodeShape()
    : Electrode() {
    size_ = 0.0;
//    minRadius_ = 0.0;
    mID_  = -1;
}

ElectrodeShape::ElectrodeShape(const RVector3 & pos)
    : Electrode(pos) {
    size_ = 0.0;
    //minRadius_ = 0.0;
    mID_  = -1;
}

ElectrodeShape::~ElectrodeShape(){
}

ElectrodeShapeNode::ElectrodeShapeNode(Node & node)
    : ElectrodeShape(node.pos()){
    valid_ = true;
    this->setNode(node);
}

ElectrodeShapeNode::~ElectrodeShapeNode(){
    // WARNING leads to mem leaks
    // it may happen that the node for this entity_ is not more existent
    // if the mesh is deleted but the electrodes not .. fop not
    // need reference counting here
//     if (entity_) delete entity_;
}

void ElectrodeShapeNode::setNode(Node & node){
    node_ = & node;
    entity_ = new NodeBoundary(node);
    mID_  = node.id();
//    minRadius_ = 0.0;
}

double ElectrodeShapeNode::geomMeanCellAttributes() const {
    std::set< Cell * > cellsPerNode(node_->cellSet());

    std::vector < double > tmpRhos;
//     for_each(cellsPerNode.begin(), cellsPerNode.end(),
//               std::back_inserter< double >(tmpRhos, std::mem_fun(& Cell::attribute)));

    if (cellsPerNode.size() == 0){
        std::cout << *this->node() << std::endl;
        THROW_TO_IMPL
    }

    for (std::set< Cell * >::iterator it = cellsPerNode.begin(); it != cellsPerNode.end(); it++){
        tmpRhos.push_back((*it)->attribute());
    }

    double val = geometricMean(tmpRhos);
    if (1.0 - val / GIMLI::max(tmpRhos) < 0.1) return val;
    return geometricMean(tmpRhos);

    std::cout << this->id() << " double ElectrodeShapeNode::geomMeanCellAttributes()"
            << " geom. mean: " << geometricMean(tmpRhos)
            << " max: " << GIMLI::max(tmpRhos)
            << std::endl;
    return GIMLI::max(tmpRhos);
}

double ElectrodeShapeNode::pot(const RVector & sol) const {
    return sol[node_->id()];
}

void ElectrodeShapeNode::assembleRHS(RVector & rhs, double value, uint matrixSize) const {
    if (node_ && matrixSize == rhs.size()){
        if (node_->id() > -1 && node_->id() < (int)rhs.size()){
            rhs[node_->id()] = value;
        } else {
            std::stringstream str1; str1 << WHERE_AM_I << " nodeID or rhs.size() invalid"
                                              << node_->id() << ", " << rhs.size() << std::endl;
            throwLengthError(str1.str());
        }
    } else {
        if (this->id() > -1 && (matrixSize + this->id()) < rhs.size()){
            rhs[matrixSize + this->id()] = value;
        } else {
            std::cerr << WHERE_AM_I << " don't know what to do " << std::endl;
            std::cerr << "Electrode-id() out of range: " << this->id() << " " << matrixSize << " "
                        << rhs.size() << std::endl;
        }
    }
}


void ElectrodeShapeNode::setSingValue(RVector & sol, double scale, double k) const {
    double minRadius = 0.0;

    if (minRadius < TOLERANCE){
        std::set< Cell * > cellsPerNode(node_->cellSet());
        std::set< Node * > neighNodes;
        for (std::set< Cell * >::iterator it = cellsPerNode.begin(); it != cellsPerNode.end(); it++){
            for (uint i = 0; i < (*it)->nodeCount(); i ++){
                neighNodes.insert(&(*it)->node(i));
            }
        }
        neighNodes.erase(node_);
        minRadius = MAX_DOUBLE;
        for (std::set< Node * >::iterator it = neighNodes.begin(); it != neighNodes.end(); it++){
            minRadius = min(minRadius, node_->pos().dist((*it)->pos()));
        }
       // std::cout << "minRadius_: " << minRadius_ << std::endl;
    }

    if (mID_ > -1.0){
        if (k > 0.0){
            if (::fabs(scale) < TOLERANCE){
                sol[mID_] = besselK0(minRadius / 6.0 * k) / (PI);
            } else {
                sol[mID_] = scale * besselK0(minRadius / 6.0 * k) / (PI) * geomMeanCellAttributes();
            }
            //std::cout << geomMeanCellAttributes() << " " << sol[mID_] << std::endl;
        } else {
            if (::fabs(scale) < TOLERANCE){
                sol[mID_] = 1.0 / (2.0 * PI * minRadius / 2.0);
            } else {
                sol[mID_] = scale * 1.0 / (2.0 * PI * minRadius / 2.0) * geomMeanCellAttributes();
            }
        }
    }
}

ElectrodeShapeNodesWithBypass::ElectrodeShapeNodesWithBypass(std::vector < Node * > & nodes)
    : ElectrodeShapeNode(*nodes[0]), nodes_(nodes){
}

ElectrodeShapeNodesWithBypass::~ElectrodeShapeNodesWithBypass(){
}

ElectrodeShapeEntity::ElectrodeShapeEntity(MeshEntity & entity, const RVector3 & pos)
    : ElectrodeShape(pos), entity_ (& entity){
    size_ = entity_->shape().domainSize();
    valid_ = true;
}

ElectrodeShapeEntity::~ElectrodeShapeEntity(){
}

void ElectrodeShapeEntity::setSingValue(RVector & sol, double scale, double k) const{

    //** check if pos_ hit node within tolerance

    int mID = -1;
    double minRadius = 0;
    for (uint i = 0; i < entity_->nodeCount(); i ++){
        Node * n = &entity_->node(i);

        if (pos_.distance(n->pos()) < 1e-4){
            mID = n->id();
            std::set< Cell * > cellsPerNode(n->cellSet());
            std::set< Node * > neighNodes;
            for (std::set< Cell * >::iterator it = cellsPerNode.begin(); it != cellsPerNode.end(); it++){
                for (uint i = 0; i < (*it)->nodeCount(); i ++){
                    neighNodes.insert(&(*it)->node(i));
                }
            }
            neighNodes.erase(n);
            minRadius = MAX_DOUBLE;
            for (std::set< Node * >::iterator it = neighNodes.begin(); it != neighNodes.end(); it++){
                minRadius = min(minRadius, n->pos().dist((*it)->pos()));
            }
        }
       // std::cout << "minRadius_: " << minRadius_ << std::endl;
    }

    if (mID > -1.0){
        if (k > 0.0){
            if (::fabs(scale) < TOLERANCE){
                sol[mID] = besselK0(minRadius / 6.0 * k) / (PI);
            } else {
                sol[mID] = scale * besselK0(minRadius / 6.0 * k) / (PI) * geomMeanCellAttributes();
            }
            //std::cout << geomMeanCellAttributes() << " " << sol[mID_] << std::endl;
        } else {
            if (::fabs(scale) < TOLERANCE){
                sol[mID] = 1.0 / (2.0 * PI * minRadius / 2.0);
            } else {
                sol[mID] = scale * 1.0 / (2.0 * PI * minRadius / 2.0) * geomMeanCellAttributes();
            }
        }
    }
}

double ElectrodeShapeEntity::geomMeanCellAttributes() const {

    //** this can be critical if the position is on a boundary or a node of the entity,
    if (entity_->parentType() == MESH_BOUNDARY_RTTI){
        Cell * le = dynamic_cast < Boundary * > (entity_)->leftCell();
        Cell * ri = dynamic_cast < Boundary * > (entity_)->rightCell();
        if (le && ri) return (le->attribute() + ri->attribute()) / 2.0;
        if (le) return le->attribute();
        if (ri) return ri->attribute();
        throwError(WHERE_AM_I + " electrode invalid. ");
    } else if (entity_->parentType() == MESH_CELL_RTTI){
        return dynamic_cast < Cell * > (entity_)->attribute();
    }

    CERR_TO_IMPL
    return 0.0;
}

double ElectrodeShapeEntity::pot(const RVector & sol) const {
    return entity_->pot(pos_, sol);
}

void ElectrodeShapeEntity::assembleRHS(RVector & rhs, double value, uint matrixSize) const {

    //** source function is \dirac(x-pos), \int f(x) dirac(x-pos)=f(pos)
    //** so right hand side entries will be shapefunctions(pos)
    if (valid_){
        if (entity_){
            rhs.setVal(entity_->N(entity_->shape().rst(pos_)), entity_->ids());
        } else {
            throwError(WHERE_AM_I + " no entity given");
        }
    } else {
        throwError(WHERE_AM_I + " electrode is not valid");
    }
}

ElectrodeShapeDomain::ElectrodeShapeDomain(const std::vector < MeshEntity * > & entities)
    : ElectrodeShape(), entities_(entities) {
    std::set < Node * > nodes;
    for (uint i = 0; i < entities_.size(); i ++) {
        size_ += entities_[i]->shape().domainSize();
        for (uint j = 0; j < entities_[i]->nodeCount(); j ++){
            nodes.insert(&entities_[i]->node(j));
        }
    }

    for (std::set < Node * >::iterator it = nodes.begin(); it != nodes.end(); it ++){
        pos_ += (*it)->pos();
    }

    pos_ /= nodes.size();

    //**! Pos muss nicht zwangl�ufig == dem center des PLC entsprechen da tetgen ja unsymmetrisch mehr Knoten reingebaut haben kann. und durch mittelung nicht zwangsl�ufig das center der einh�llenden bestimmt wird.
    valid_ = true;
//     std::cout << "ElectrodeShapeDomain::ElectrodeShapeDomain(const std::vector < MeshEnity * > & entities): "
//     << pos_ << std::endl;
}

ElectrodeShapeDomain::ElectrodeShapeDomain(const std::vector < Boundary * > & entities)
    : ElectrodeShape() {
    for (uint i = 0; i < entities.size(); i ++) {
        entities_.push_back(entities[i]);
        size_ += entities[i]->shape().domainSize();
        pos_  += entities[i]->shape().center();
    }
    pos_ /= entities.size();
    valid_ = true;
//     std::cout << "ElectrodeShapeDomain::ElectrodeShapeDomain(const std::vector < Boundary * > & entities): "
//     << pos_ << std::endl;
}

ElectrodeShapeDomain::~ElectrodeShapeDomain(){
}

double ElectrodeShapeDomain::geomMeanCellAttributes() const {
    double weighted = 0.0;
    for (uint i = 0; i < entities_.size(); i ++) {
        if (entities_[i]->parentType() == MESH_BOUNDARY_RTTI){
            Cell * le = dynamic_cast < Boundary * > (entities_[i])->leftCell();
            Cell * ri = dynamic_cast < Boundary * > (entities_[i])->rightCell();
            if (le && ri) {
                //** cem ist im moment nur f�r �ussere boundary faces definiert
                CERR_TO_IMPL
                return 0.0;
            }
            if (le) {
                weighted += le->attribute() * entities_[i]->shape().domainSize() / size_;
            } else if (ri) {
                weighted += ri->attribute() * entities_[i]->shape().domainSize() / size_;
            } else {
                std::cerr << WHERE_AM_I << " WARNING! No cell found " << std::endl;
            }
        }
    }
    return weighted;
}

void ElectrodeShapeDomain::assembleRHS(RVector & rhs, double value, uint matrixSize) const {
    if (rhs.size() > matrixSize){
      //** here we assume CEM;
        if (this->id() > -1 && (matrixSize + this->id()) < rhs.size()){
            rhs[matrixSize + this->id()] = value;
        } else {
            std::cerr << WHERE_AM_I << " don't know what to do " << std::endl;
            std::cerr << "Electrode-id() out of range: " << this->id() << " " << matrixSize << " "
                        << rhs.size() << std::endl;
        }
    } else {
        std::cerr << WHERE_AM_I << " this makes no sense, calculate complete electrode model" << std::endl;
        std::cerr << "Electrode-id() out of range: " << this->id() << " " << matrixSize << " "
                        << rhs.size() << std::endl;
    }
}

double ElectrodeShapeDomain::pot(const RVector & sol) const {
    //*** temp Hack ********
//     double pot = 0.0;
//     ElementMatrix < double > Se, Stmp;
//
//     for (uint i = 0; i < entities_.size(); i ++) {
//         Stmp = Se.u(*entities_[i]);
//
//  //           pot += Stmp * sol;
//         //std::cout << entities_[i]->shape().jacobianDeterminant() << " "  << Stmp << std::endl;
//         std::valarray < std::valarray < double > > mat(Stmp.mat());
//         std::vector < uint > idx(Stmp.idx());
//         double utmp=0.0;
//         for (uint j = 0; j < mat[0].size(); j ++){
//                     //std::cout << idx[j] << " " << mat[0][j] << " " << sol[idx[j]] << std::endl;
//             utmp += mat[0][j] * sol[idx[j]];
//         }
//         pot += utmp;
//     }
//     return pot / size_;
    //*** temp Hack ********

    throwError(" The potential can not determined directly from ElectrodeShapeDomain use the matrix extension from complete electrode model instead");
    return 0.0;
}

} //namespace BERT
