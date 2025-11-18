/******************************************************************************
 *   Copyright (C) 2012-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_GRAVIMETRY__H
#define _GIMLI_GRAVIMETRY__H

#include "gimli.h"
#include "modellingbase.h"

namespace GIMLI {

//! Modelling class for gravimetry calculation using polygon integration
/*! Modelling class for gravimetry calculation using polygon integration */
class DLLEXPORT GravimetryModelling : public ModellingBase {
public:
    GravimetryModelling( Mesh & mesh, DataContainer & dataContainer, bool verbose = false );

    virtual ~GravimetryModelling() { }

    RVector createDefaultStartModel( );

    /*! Interface. Calculate response */
    virtual RVector response( const RVector & slowness );

    /*! Interface. */
    virtual void createJacobian( const RVector & slowness );

    /*! Interface. */
    virtual void initJacobian( );

protected:
};



/*! Only for a small TPOC for 2d gravimetry after WonBevis1987
 Do not use until u know what u do. */
DLLEXPORT double lineIntegraldGdz( const RVector3 & p1, const RVector3 & p2 );

/*! Do not use until u know what u do. */
DLLEXPORT RVector calcGBounds( const std::vector< RVector3 > & pos, const Mesh & mesh, const RVector & model );

/*! Do not use until u know what u do. */
DLLEXPORT RVector calcGCells( const std::vector< RVector3 > & pos, const Mesh & mesh, const RVector & model, uint nInt = 0 );

} //namespace GIMLI

#endif
