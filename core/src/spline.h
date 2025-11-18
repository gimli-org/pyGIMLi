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

#ifndef _GIMLI_SPLINE__H
#define _GIMLI_SPLINE__H

#include "gimli.h"
#include "pos.h"

#include <vector>

namespace GIMLI{

class DLLEXPORT CubicFunct{
public:
    /*! define function ax^3+bx^2+cx+d*/
    CubicFunct( const double a = 0.0, const double b = 0.0, const double c = 0.0, const double d = 0.0)
        : a_( a ), b_( b ), c_( c ), d_( d ) {}

    CubicFunct( const CubicFunct & C ) : a_( C.a_ ), b_( C.b_ ), c_( C.c_ ), d_( C.d_ ) {}

    CubicFunct operator = ( const CubicFunct & C ) {
        if ( this != & C ){
            a_ = C.a_; b_ = C.b_; c_ = C.c_; d_ = C.d_;
        } return *this;
    }

    inline double operator()( double t ) const { return t * (t * (t * a_ + b_) + c_) + d_; }

    inline double val( double t ) const { return (*this)( t ); }

    friend bool operator == ( const CubicFunct & a , const CubicFunct & b );

    double a_, b_, c_, d_;
};

inline bool operator == ( const CubicFunct & a , const CubicFunct & b ){
  if ( std::fabs( a.a_ - b.a_ ) < TOLERANCE &&
       std::fabs( a.b_ - b.b_ ) < TOLERANCE &&
       std::fabs( a.c_ - b.c_ ) < TOLERANCE &&
       std::fabs( a.d_ - b.d_ ) < TOLERANCE ) return true; else return false;
}

DLLEXPORT std::vector < CubicFunct > calcNaturalCubic( const std::vector < double > & x );
DLLEXPORT std::vector < CubicFunct > calcNaturalCubicClosed( const std::vector < double > & x );

/*! Create a vector of RealPos from input points with cubic spline interpolation. The edge between the input points is subdivided into nSegments.*/
DLLEXPORT std::vector < RVector3 > createSpline( const std::vector < RVector3 > & input, int nSegments, bool close );

/*! Create a vector of RealPos from input points with cubic spline interpolation. The edge between the input points is subdivided into 3 segments within the lokal distance localDX to the inputpoints.*/
DLLEXPORT std::vector < RVector3 > createSplineLocalDX( const std::vector < RVector3 > & input, double localDX, bool close );

} // namespace GIMLI

#endif // _GIMLI_SPLINE__H

