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

#ifndef _GIMLI_NUMERICBASE__H
#define _GIMLI_NUMERICBASE__H

#include "gimli.h"
#include <cmath>

namespace GIMLI{

#ifndef PI
    #define PI 3.141592653589793238462643383279502884197169399375105820974944592308
    #define PI2 6.2831853071795864769252867665590057683942387987502116419498891846
    #define PI_2 1.5707963267948967384626433682795028841971193993751058209749445923
#endif

template < class ValueType > ValueType round(const ValueType & v, ValueType tol){ return ::rint(v / tol) * tol; }

/*! DEPRECATED use logDropTol.*/
DLLEXPORT RVector logTransDropTol(const RVector & data, double logdrop=1e-6,
                                  bool normalize=true);

/*! log10 scale of data with negative values.
 * Conserve sign alternating a drop tolerance.*/
DLLEXPORT RVector logDropTol(const RVector & data, double logdrop=1e-6,
                             bool normalize=true);

/*! Converts a degree value to radian.*/
template < class ValueType > ValueType degToRad(const ValueType & deg){ return deg * (2.0 * PI) / 360.0; }

/*! Converts a radian value to degree.*/
template < class ValueType > ValueType radToDeg(const ValueType & rad){ return rad * 360.0 / (2.0 * PI); }

template < class T > T powInt(const T & a, uint dim){
    switch (dim){
        case 0 : return (T)1;
        case 1 : return a;
        case 2 : return a*a;
        case 3 : return a*a*a;
        case 4 : return a*a*a*a;
        case 5 : return a*a*a*a*a;
        case 6 : return a*a*a*a*a*a;
        default: return (T)std::pow((float)a, (float)dim);
    }
}


/*!  Given the lower and upper limits of integration x1 and x2 and given n,
  this routine returns vector x(0..n-1) and w(0..n-1) of length n,
  containing the abscissas and weights of the Gauss-Legendre
  n-point quadrature formula. */
DLLEXPORT void GaussLegendre(double x1, double x2, uint n, RVector & x, RVector & w);

/*!  Given alpha = 0.0, the parameter alpha of the Laguerre polynomials, this routine
  returns vector x[0..n-1] and w[0..n-1] containing the abscissas and weights
  of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa
  is returned in x[ 0 ], the largest in x[ n-1 ]. */
DLLEXPORT void GaussLaguerre(uint n, RVector & x, RVector & w);


//! Caluculate modified Bessel function of the first kind
/*! See Abramowitz: Handbook of math. functions; */
//DLLEXPORT double besselI0(double x);
template < class ValueType > ValueType besselI0(const ValueType & x) {
    ValueType ax = std::fabs(x), result = 0.0, y = 0.0;

    if (ax < 3.75) {
        y = x / 3.75, y = y * y ;
        result = 1.0 + y * (3.5156229 +
             y * (3.0899424 +
                   y * (1.2067492 +
                     y * (0.2659732 +
                       y * (0.360768e-1
                         + y * 0.45813e-2)))));
    } else {
        y = 3.75 / ax;
        result = (std::exp(ax) / std::sqrt(ax)) * (0.39894228 +
                        y * (0.1328592e-1 +
                          y * (0.225319e-2 +
                            y * (-0.157565e-2 +
                                  y * (0.916281e-2 +
                                    y * (-0.2057706e-1 +
                                      y * (0.2635537e-1 +
                                        y * (-0.1647633e-1 +
                                              y * 0.392377e-2))))))));
    }
    return result;
}

//! Caluculate modified Bessel function of the first kind
/*! See Abramowitz: Handbook of math. functions */
//DLLEXPORT double besselI1(double x);
template < class ValueType > ValueType besselI1(const ValueType & x) {
  ValueType ax = std::fabs(x), result = 0.0, y = 0.0;

  if (ax < 3.75) {
    y = x / 3.75, y =y * y;
    result = ax * (0.5 +
            y * (0.87890594 +
              y * (0.51498869 +
                y * (0.15084934 +
                      y * (0.2658733e-1 +
                        y * (0.301532e-2 +
                          y * 0.32411e-3))))));
  } else {
    y = 3.75 / ax;
    result = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 -y * 0.420059e-2));
    result = 0.39894228   + y * (-0.3988024e-1 + y * (-0.362018e-2 +
                            y * (0.163801e-2 + y * (-0.1031555e-1 + y * result)))) ;
    result *= (std::exp(ax) / std::sqrt(ax));
  }
  return x < 0.0 ? -result : result;
}


//! Caluculate modified Bessel function of the second kind
/*! See Abramowitz: Handbook of math. functions */
//DLLEXPORT double besselK0(double x);
template < class ValueType > ValueType besselK0(const ValueType & x){
    ValueType y = 0.0, result = 0.0;

    if (x <= 2.0) {
        y = x * x / 4.0;
        result = (-std::log(x / 2.0) * besselI0(x)) + (-0.57721566
                          + y * (0.42278420
                              + y * (0.23069756
                                  + y * (0.3488590e-1
                                      + y * (0.262698e-2
                                          + y * (0.10750e-3
                                              + y * 0.74e-5))))));
    } else {
        y = 2.0 / x;
        result = (std::exp(-x) / std::sqrt(x)) * (1.25331414
                            + y * (-0.7832358e-1 + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
                            + y * (-0.251540e-2
                            + y * 0.53208e-3))))));
    }
    return result;
}

//! Caluculate modified Bessel function of the second kind
/*! See Abramowitz: Handbook of math. functions */
//DLLEXPORT double besselK1(double x);
template < class ValueType > ValueType besselK1(const ValueType & x){
    ValueType y = 0.0, result = 0.0;

    if (x <= 2.0) {
        y = x * x / 4.0;
        result = (std::log(x / 2.0) * besselI1(x)) +
        (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 +
                                     y * (-0.1919402e-1 +
                                           y * (-0.110404e-2 +
                                             y * (-0.4686e-4))))))) ;
    } else {
        y = 2.0 / x;
        result = (std::exp(-x) / std::sqrt(x)) * (1.25331414 + y * (0.23498619 +
                                  y * (-0.3655620e-1 +
                                    y * (0.1504268e-1 +
                                      y * (-0.780353e-2 +
                                        y * (0.325614e-2 +
                                              y * (-0.68245e-3)))))));
    }
    return result;
}

//*!Spherical tangential coordinates to geocentric Cartesian coordinates
/*!
 * Convert vector field from spherical tangential coordinates
    (radial, Latitude/theta(north/south), Longitude/phi(east/west))

        \vec{V(1, lat, lon)} = V[0] * \vec{unit_r} +  V[1] * \vec{unit_theta} + V[2] * \vec{unit_phi}

        to geocentric Cartesian coordinates x/y/z:

        \vec{F(x, y, z)} = F[0] * \vec{unit_x} +  F[1] * \vec{unit_y} + F[2] * \vec{unit_z}

        Transformation via rotation matrix S
        F(x,y,z) = S * V(r,\theta,\phi)

        J = S * (1, r, r cos th)

        J (\dx, \dy, \dz) / (\dr, \d th, \d ph)

        x = r * cos ph * cos th
        y = r * sin ph * cos th
        z = r * sin th

        th = latitude degrees -pi/2 .. pi/2, 90 = north pole
        ph = longitude degrees  -pi .. pi .. west - east


    Inputs
        V B in radial direction or Magnetic field strength (B)
        lon longitude degrees (in degrees from -180 to 180)
        lat Latitude measured positive north from equator (in degrees from south pole [-90 .. 0 .. 90] north pole)
    Outputs
        F
*/
DLLEXPORT RVector3 sphTangential2Initerial(const RVector3 &V, double lat, double lon);

/*! Lineintegral after WonBevis.*/
DLLEXPORT void lineIntegralZ_WonBevis(const RVector3 &p1, const RVector3 &p2,
                                      RVector3 & dg, RVector3 & dgz);
/*! Lineintegral after WonBevis. Only return dgz.*/
DLLEXPORT double lineIntegralZ_WonBevis(const RVector3 &p1, const RVector3 &p2);

} // namespace GIMLI

#endif // _GIMLI_NUMERICBASE__H
