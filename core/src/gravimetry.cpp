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

#include "gravimetry.h"

#include "integration.h"
#include "mesh.h"
#include "pos.h"
#include "shape.h"

#include <cmath>

namespace GIMLI {

GravimetryModelling::GravimetryModelling( Mesh & mesh, DataContainer & dataContainer, bool verbose ){
}

RVector GravimetryModelling::createDefaultStartModel( ){
    RVector ret;
    THROW_TO_IMPL
    return ret;
}

RVector GravimetryModelling::response( const RVector & model ){
    RVector ret;
    THROW_TO_IMPL
    return ret;
}

void GravimetryModelling::createJacobian( const RVector & model ){
    THROW_TO_IMPL
}

void GravimetryModelling::initJacobian( ){
    THROW_TO_IMPL
}


double lineIntegraldGdz( const RVector3 & p1, const RVector3 & p2 ){
    double x1 = p1[ 0 ], z1 = p1[ 1 ];
    double x2 = p2[ 0 ], z2 = p2[ 1 ];

    if ( ( ::fabs( x1 ) < TOLERANCE ) && ( ::fabs( z1 ) < TOLERANCE ) ) return 0.0;
    if ( ( ::fabs( x2 ) < TOLERANCE ) && ( ::fabs( z2 ) < TOLERANCE ) ) return 0.0;

    double theta1 = ::atan2( z1, x1 );
    double theta2 = ::atan2( z2, x2 );

    double r1 = ::sqrt( x1*x1 + z1*z1 );
    double r2 = ::sqrt( x2*x2 + z2*z2 );

    // z-component of gravitational field
    double Z = 0.0;
    // x-component of gravitational field
    // double X = 0.0;

    if ( sign( z1 ) != sign( z2 ) ){

        if ( ( x1*z2 < x2*z1 ) && ( z2 >=0.0 ) ) {
            theta1 += PI2;
        } else if( ( x1*z2 > x2*z1 ) && ( z1 >=0.0 ) ){
            theta2 += PI2;
        } else if( ::fabs( x1*z2 - x2*z1 ) < TOLERANCE ){
            // X = Z = 0
            return 0.0;
        }
    }

    if ( ::fabs( x1 - x2 ) < TOLERANCE ){ // case 3
        Z = x1 * ::log( r2 / r1 );
        // X = -x1 * ( theta1 - theta2 );
    } else { // default

        double B = ( z2 - z1 ) / ( x2 - x1 );
        double A = ( ( x2 - x1 ) * ( x1 * z2 - x2 * z1 ) ) / ( ( x2 - x1 )*( x2 - x1 ) + ( z2 - z1 )*( z2 - z1 ) );

        Z = A * ( ( theta1 - theta2 ) + B * ::log( r2 / r1 ) );
        // X = A * ( -( theta1 - theta2 ) B + ::log( r2 / r1 ) );
    }

    return Z;
}

RVector calcGBounds( const std::vector< RVector3 > & pos, const Mesh & mesh, const RVector & model ){
    /*! Ensure neighborInfos() */
    RMatrix Jacobian( pos.size(), mesh.cellCount() );

    Jacobian *= 0.;

    for ( uint i = 0; i < pos.size(); i ++ ){
        for ( std::vector< Boundary * >::const_iterator it = mesh.boundaries().begin(); it != mesh.boundaries().end(); it ++ ){
            Boundary *b = *it;
            double Z = lineIntegraldGdz( b->node( 0 ).pos() - pos[ i ], b->node( 1 ).pos() - pos[ i ] );

            if ( b->leftCell() ) Jacobian[ i ][ b->leftCell()->id() ] = Jacobian[ i ][ b->leftCell()->id() ] - Z;
            if ( b->rightCell() ) Jacobian[ i ][ b->rightCell()->id() ] = Jacobian[ i ][ b->rightCell()->id() ] + Z;
        }
    }

    return Jacobian * model * 2.0 * 6.67384e-11 * 1e5;
}

double f_gz( const RVector3 & x, const RVector3 & p ){
    double r = x.dist( p );
    return (p[1]-x[1]) / (r*r);
}

RVector calcGCells( const std::vector< RVector3 > & pos, const Mesh & mesh, const RVector & model, uint nInt ){
    /*! Ensure neighborInfos() */
    RMatrix Jacobian( pos.size(), mesh.cellCount() );

    Jacobian *= 0.;

    for ( uint i = 0; i < pos.size(); i ++ ){
        for ( std::vector< Cell * >::const_iterator it = mesh.cells().begin(); it != mesh.cells().end(); it ++ ){
            Cell *c = *it;
            double Z = 0.;
            if ( nInt == 0 ){
                for ( uint j = 0; j < c->nodeCount(); j ++ ){
                    // negative Z because all cells are numbered counterclockwise
                    Z -= 2.0 * lineIntegraldGdz( c->node( j ).pos() - pos[ i ], c->node( (j+1)%c->nodeCount() ).pos() - pos[ i ] );
                }
            } else {
                for ( uint j = 0; j < IntegrationRules::instance().triAbscissa( nInt ).size(); j ++ ){
                    Z += IntegrationRules::instance().triWeights( nInt )[ j ] *
                                    f_gz( c->shape().xyz( IntegrationRules::instance().triAbscissa( nInt )[ j ] ), pos[ i ] );
                }
            }

            Jacobian[ i ][ c->id() ] = -Z;
        }
    }

    return Jacobian * model * 6.67384e-11 * 1e5;
}

} // namespace GIMLI{
