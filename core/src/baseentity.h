/******************************************************************************
 *   Copyright (C) 2005-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_BASEENTITY__H
#define _GIMLI_BASEENTITY__H

#include "gimli.h"

namespace GIMLI{

  //! Base Entity
  /*! The is a base entity which holds some basic informations. An identification
  * number for the object, a validation status flag and a runtime type information (RTTI).
  * This information will be inherited in many objects i.e. the mesh and data entitys. */
class DLLEXPORT BaseEntity{
    public:
    BaseEntity()
        : id_(-1), valid_(false), marker_(0) { }

    BaseEntity(const BaseEntity & ent)
        : id_(ent.id()), valid_(ent.valid()), marker_(ent.marker()){
    }

    BaseEntity & operator = (const BaseEntity & ent){
        if (this != &ent){
            id_ = ent.id();
            valid_ = ent.valid();
            marker_ = ent.marker();
        } return * this;
    }

    virtual ~BaseEntity(){}

    /*! Return entity rtti value. */
    inline virtual uint rtti() const { return MESH_BASEENTITY_RTTI; }

    /*! Return entity valid status. */
    inline virtual bool valid() const { return valid_; }

    /*! Set entity valid status. */
    inline virtual void setValid(bool valid) { valid_ = valid ; }

    /*! Returns the entity id. */
    inline int id() const { return id_; }

    /*! Set the entity id. */
    inline void setId(int id) { id_ = id ; }

    inline void setMarker(int marker) { marker_ = marker; }

    inline int marker() const { return marker_; }

    /*! Userflag to mark the entity for something you want.
     * This will be used internal while cell search so be carefully. */
    inline void setTagged(bool tagged){ tagged_ = tagged; }

    /*! Untag the cell */
    inline void untag() { setTagged(false); }

    /*! Tag the cell */
    inline void tag() { setTagged(true); }

    /*! Return true if the cell is tagged */
    inline bool tagged() const { return tagged_; }

protected:
    int id_;

    bool valid_;

    int marker_;

    bool tagged_;

}; // class BaseEntity;

}  // namespace GIMLI{


#endif // _GIMLI_BASEENTITY__H
