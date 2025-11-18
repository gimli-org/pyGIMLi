/******************************************************************************
 *   Copyright (C) 2006-2024 by the resistivity.net development team          *
 *   Carsten R�cker carsten@resistivity.net                                   *
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

#ifndef _BERT_ELECTRODE__H
#define _BERT_ELECTRODE__H

#include "bert.h"

#include <baseentity.h>
#include <pos.h>

namespace GIMLI{

// class Electrode -> ElectrodeNode -> ElectrodeRVector3 -> ElectrodeBoundary -> ElectrodeDomain

/*! A simple Electrode. */
class DLLEXPORT Electrode : public GIMLI::BaseEntity{
public:
  /*! Constructor, builds an empty and non valid electrode.*/
    Electrode();

    Electrode(const RVector3 & pos, int id = -1);

    Electrode(double x, double y, double z);

    virtual ~Electrode();

    Electrode(const Electrode & el);

    Electrode & operator = (const Electrode & el);

    bool operator == (const Electrode & e) const { return pos_ == e.pos(); }

    inline void setPos(const RVector3 & pos){ pos_ = pos; }

    inline const RVector3 & pos() const { return pos_; }

protected:
    RVector3 pos_;
};

//! Abstract class for an electrode with a shape, which is required for modelling
class DLLEXPORT ElectrodeShape : public Electrode{
public:
    ElectrodeShape();

    ElectrodeShape(const RVector3 & pos);

    virtual ~ElectrodeShape();

    virtual std::vector < MeshEntity * > entities() const {
        THROW_TO_IMPL; return std::vector < MeshEntity * >();
    }

    virtual double domainSize() const { return size_; }

    //virtual double minRadius() const { return minRadius_; }

    /*! Default method */
    virtual double geomMeanCellAttributes() const { THROW_TO_IMPL; return 0.0; }

    /*! Default method to collect to potential value corresponding to this electrode */
    virtual double pot(const RVector & sol) const { THROW_TO_IMPL; return 0.0; }

    /*! Default method */
    virtual void assembleRHS(RVector & rhs, double value, uint matrixSize) const {
        THROW_TO_IMPL; return;
    }

    /*! Default method */
    virtual void setSingValue(RVector & sol, double scale, double k) const {

    };

    /*! Sets the position within the system matrix, e.g. the node id for ElectrodeShapeNode.*/
    void setMID(uint id){ mID_ = id; }

    /*! Return matrix ID, i.e. the discrete index within the system of matrix or vector. maybee also dof-position, -1 for unknown.*/
    virtual int mID() const { return mID_; }

protected:

    //! size of the domain
    double size_;
    //double minRadius_;
    int mID_;
};

//! Electrode that is represented by a node
class DLLEXPORT ElectrodeShapeNode : public ElectrodeShape{
public:
    ElectrodeShapeNode(Node & node);

    virtual ~ElectrodeShapeNode();

    void setNode(Node & node);

    inline const Node * node() const { return node_; }

    virtual std::vector < MeshEntity * > entities() const {
        std::vector < MeshEntity * > ents;
        ents.push_back(entity_);
        return ents;
    }

    virtual double geomMeanCellAttributes() const;

    /*! Collect the potential value corresponding to this nodes id, this is the counterpart to
    assembleRHS. */
    virtual double pot(const RVector & sol) const;

    virtual void assembleRHS(RVector & rhs, double value, uint matrixSize) const;

    virtual void setSingValue(RVector & sol, double scale = 1.0, double k = 0.0) const;

protected:
    Node * node_;
    MeshEntity * entity_;
};


//! Electrode that is represented by a list of nodes short circuited by a bypass.
class DLLEXPORT ElectrodeShapeNodesWithBypass : public ElectrodeShapeNode {
public:
    /*! The corresponding reference Node is the first of a given vector of nodes. */
    ElectrodeShapeNodesWithBypass(std::vector < Node * > & nodes);

    virtual ~ElectrodeShapeNodesWithBypass();

protected:
    std::vector < Node * > nodes_;
};




//! Electrodeshape is singular source within a mesh entity (boundary or cell)
class DLLEXPORT ElectrodeShapeEntity : public ElectrodeShape{
public:
    ElectrodeShapeEntity(MeshEntity & entity, const RVector3 & pos);

    virtual ~ElectrodeShapeEntity();

    inline void setEntity(MeshEntity & entity) { entity_ = & entity; }
    inline MeshEntity * entity() const { return entity_; }

    virtual std::vector < MeshEntity * > entities() const {
        std::vector < MeshEntity * > ents;
        ents.push_back(entity_);
        return ents;
    }

    virtual double geomMeanCellAttributes() const;

    /*! Collect the potential value corresponding to this nodes id */
    virtual double pot(const RVector & sol) const;

    virtual void assembleRHS(RVector & rhs, double value, uint matrixSize) const;

    virtual void setSingValue(RVector & sol, double scale, double k) const;

protected:
    MeshEntity * entity_;
};


//! Electrodeshape is a domain, e.g. the boundary of a complicated geoemtry for the complete electrode model
class DLLEXPORT ElectrodeShapeDomain : public ElectrodeShape{
public:
    ElectrodeShapeDomain(const std::vector < MeshEntity * > & entities);

    /*! Temporary for pygimli
        to fit: g.ElectrodeShapeDomain(mesh.findBoundaryByMarker(g.MARKER_BOUND_ELECTRODE))
    */
    ElectrodeShapeDomain(const std::vector < Boundary * > & entities);

    virtual ~ElectrodeShapeDomain();

    inline void setEntities(std::vector < MeshEntity * > & entities) { entities_ = entities; }

    virtual std::vector < MeshEntity * > entities() const { return entities_; }

    virtual double geomMeanCellAttributes() const;

    virtual double pot(const RVector & sol) const;

    virtual void assembleRHS(RVector & rhs, double value, uint matrixSize) const;

protected:
    std::vector < MeshEntity * > entities_;
};

} //namespace BERT

#endif // _BERT_ELECTRODE__H
