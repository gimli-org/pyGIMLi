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

#include "spline.h"

namespace GIMLI{

std::vector < RVector3 > createSpline( const std::vector < RVector3 > & input, int nSegments, bool close ){

  std::vector < double > inX( input.size() ), inY( input.size() );
  for ( size_t i = 0; i < input.size(); i ++ ){
    inX[ i ] = input[ i ][ 0 ];
    inY[ i ] = input[ i ][ 1 ];
  }

  std::vector < CubicFunct > X,Y;
  if ( close ){
    X = calcNaturalCubicClosed( inX );
    Y = calcNaturalCubicClosed( inY );
  } else {
    X = calcNaturalCubic( inX );
    Y = calcNaturalCubic( inY );
  }

  std::vector < RVector3 > output;

  for ( size_t i = 0; i < X.size(); i++) {
    for ( int j = 0; j < nSegments; j++ ) {
      double u = j / (double)nSegments;
      output.push_back( RVector3( X[ i ].val( u ), Y[ i ].val( u ), 0.0 ) );
    }
  }
  return output;
}

std::vector < RVector3 > createSplineLocalDX( const std::vector < RVector3 > & input, double localDX, bool close ){

  std::vector < double > inX( input.size() ), inY( input.size() );
  for ( unsigned int i = 0; i < input.size(); i ++ ){
    inX[ i ] = input[ i ][ 0 ];
    inY[ i ] = input[ i ][ 1 ];
  }

  std::vector < CubicFunct > X,Y;
  if ( close ){
    X = calcNaturalCubicClosed( inX );
    Y = calcNaturalCubicClosed( inY );
  } else {
    X = calcNaturalCubic( inX );
    Y = calcNaturalCubic( inY );
  }

  std::vector < RVector3 > output;
  output.push_back( RVector3( X[ 0 ].val( 0 ), Y[ 0 ].val( 0 ), 0.0 ) );

  for ( size_t i = 0; i < X.size(); i++) {
    output.push_back( RVector3( X[ i ].val( localDX ), Y[ i ].val( localDX ), 0.0 ) );
    output.push_back( RVector3( X[ i ].val( 1-localDX ), Y[ i ].val( 1-localDX ), 0.0 ) );
    output.push_back( RVector3( X[ i ].val( 1 ), Y[ i ].val( 1 ), 0.0 ) );
  }
  return output;
}

/*! FROM: http://www.cse.unsw.edu.au/~lambert/splines/
  calculates the closed natural cubic spline that interpolates
  x[0], x[1], ... x[n]
  The first segment is returned as
  C[0].a + C[0].b*u + C[0].c*u^2 + C[0].d*u^3 0<=u <1
  the other segments are in C[1], C[2], ...  C[n] */

std::vector < CubicFunct > calcNaturalCubicClosed( const std::vector < double > & x ){
  int n = x.size()-1;

    double * w = new double[ n + 1 ];
    double * v = new double[ n + 1 ];
    double * y = new double[ n + 1 ];
    double * D = new double[ n + 1 ];
  double z = 0.0, F = 0.0, G = 0.0, H = 0.0;
/* We solve the equation
   [4 1      1] [D[0]]   [3(x[1] - x[n])  ]
   |1 4 1     | |D[1]|   |3(x[2] - x[0])  |
   |  1 4 1   | | .  | = |      .         |
   |    ..... | | .  |   |      .         |
   |     1 4 1| | .  |   |3(x[n] - x[n-2])|
   [1      1 4] [D[n]]   [3(x[0] - x[n-1])]

   by decomposing the matrix into upper triangular and lower matrices
   and then back sustitution.  See Spath "Spline Algorithms for Curves
   and Surfaces" pp 19--21. The D[i] are the derivatives at the knots.
*/

 w[ 0 ] = v[ 0 ] =0;
 w[ 1 ] = v[ 1 ] = z = 1.0 / 4.0;
 y[ 0 ] = z * 3.0 * ( x[ 1 ] - x[ n ] );
 H = 4.0;
 F = 3.0 * ( x[ 0 ] - x[ n - 1 ] );
 G = 1.0;
 for ( int i = 1; i < n; i++) {
   v[ i + 1 ] = z = 1.0 / ( 4.0 - v[ i ] );
   w[ i + 1 ] = -z * w[ i ];
   y[ i ] = z * ( 3.0 * ( x[ i + 1 ] - x[ i - 1 ] ) - y[ i - 1 ] );
   H = H - G * w[ i ];
   F = F - G * y[ i - 1 ];
   G = -v[ i ] * G;
 }
 H = H - ( G + 1.0 ) * ( v[ n ] + w[ n ] );
 y[ n ] = F - ( G + 1 ) * y[ n - 1 ];


 D[ n ] = y[ n ] / H;
 D[ n - 1 ] = y[ n - 1 ] - ( v[ n ] + w[ n ] ) * D[ n ]; /* This equation is WRONG! in my copy of Spath */
 for ( int i = n - 2; i >= 0; i-- ) {
   D[ i ] = y[ i ] - v[ i + 1 ] * D[ i + 1 ] - w[ i + 1 ] * D[ n ];
 }

    std::vector < CubicFunct > C;
    for ( int i = 0; i < n; i++) {
        C.push_back( CubicFunct( 2.0 * ( x[ i ] - x[ i + 1 ] ) + D[ i ] + D[ i + 1 ],
                                3.0 * ( x[ i + 1 ] - x[ i ] ) - 2.0 * D[ i ] - D[ i + 1 ],
                                D [ i ], x[ i ] ) );
    }
    C.push_back( CubicFunct(   2.0 *( x[ n ] - x[ 0 ] ) + D[ n ] + D[ 0 ],
                                3.0 * ( x[ 0 ] - x[ n ] ) - 2.0 * D[ n ] - D[ 0 ],
                                D[ n ],
                                x[ n ] ) );
    delete [ ] w;
    delete [ ] v;
    delete [ ] y;
    delete [ ] D;

    return C;
}

std::vector < CubicFunct > calcNaturalCubic( const std::vector < double > & x ){
  int n = x.size()-1;

  /* We solve the equation
     [2 1       ] [D[0]]   [3(x[1] - x[0])  ]
     |1 4 1     | |D[1]|   |3(x[2] - x[0])  |
     |  1 4 1   | | .  | = |      .         |
     |    ..... | | .  |   |      .         |
     |     1 4 1| | .  |   |3(x[n] - x[n-2])|
     [       1 2] [D[n]]   [3(x[n] - x[n-1])]

     by using row operations to convert the matrix to upper triangular
     and then back sustitution.  The D[i] are the derivatives at the inots.
  */

 double * gamma = new double[ n + 1 ];
  gamma[ 0 ] = 0.5;
  for ( int i = 1; i < n; i++) {
    gamma[ i ] = 1.0 / ( 4 - gamma[ i - 1 ] );
  }
  gamma[ n ] = 1.0 / ( 2 - gamma[ n - 1 ] );

  double * delta = new double[ n + 1 ];

  delta[ 0 ] = 3.0 * ( x[ 1 ] - x[ 0 ] ) * gamma[ 0 ];
  for ( int i = 1; i < n; i++ ){
    delta[ i ] = ( 3.0 * ( x[ i + 1 ] - x[ i - 1 ] ) - delta[ i - 1 ] ) * gamma[ i ];
  }
  delta[ n ] = ( 3.0 * ( x[ n ] - x[ n - 1 ] ) - delta[ n - 1 ] ) * gamma[ n ];

  double * D = new double[ n + 1 ];
  D[ n ] = delta[ n ];
  for ( int i = n-1; i >= 0; i--) {
    D[ i ] = delta[ i ] - gamma[ i ] * D[ i + 1 ];
  }

    /* now compute the coefficients of the cubics */
    std::vector < CubicFunct > C;
    for ( int i = 0; i < n; i++ ){
        C.push_back( CubicFunct( 2.0 * ( x[ i ] - x[ i + 1 ] ) + D[ i ] + D[ i + 1 ],
                                3.0 * (x[ i + 1 ] - x[ i ] ) - 2.0 * D[ i ] - D[ i + 1 ],
                                D[ i ],
                                x[ i ] ) );
    }

    delete [ ] gamma;
    delete [ ] delta;
    delete [ ] D;
    return C;
}

}  //namespace GIMLI
