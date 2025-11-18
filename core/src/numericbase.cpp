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

#include "numericbase.h"
#include "vector.h"
#include "pos.h"

namespace GIMLI{

RVector logDropTol(const RVector & data, double logdrop, bool normalize){
    RVector tmp(data);
    RVector tmp2(data);
//     __MS(data)
    tmp2.setVal(logdrop, find(data == 0.0));

    for (uint i = 0; i < tmp.size(); i ++) {
        tmp[i] = max(1.0, std::fabs(tmp[i] / logdrop));
    }

    tmp = log10(tmp);
    // normalize
    if (normalize){
        double m = max(abs(tmp));
        if (m != 0.0) tmp /= m;
    }
//     __MS(tmp)
//     __MS(RVector(sign(tmp2)))
    return tmp * sign(tmp2);
}


RVector logTransDropTol(const RVector & data, double logdrop, bool normalize){
    return logDropTol(data, logdrop, normalize);
}

void GaussLaguerre(uint n, RVector & x, RVector & w){
//  taken from matlab-code (Thomas Guenther)
//  function [x, w] = gaulag(n)
//  GAULAG - Gauss-Laguerre Integration Points
//  [x,w] = gaulag(n)
//  Given alf = 0.0, the parameter alpha of the Laguerre polynomials, this routine
//  returns arrays x[1..n] and w[1..n] containing the abscissas and weights
//  of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa
//  is returned in x[1], the largest in x[n].
//  For a description of the following routines see
//  Numerical Recipes, Press et a

  if (x.size() != n) x.resize(n);
  if (w.size() != n) w.resize(n);

  double epsilon = 3.0e-11;
  int maxiter = 20;

  double z = 0.0, z1 = 0.0, p1 = 1.0, p2 = 0.0, p3 = 0.0, pp = 0.0;
  int ai = 0;

  for (size_t i = 1; i <=n; i ++){ //	Loop over desired roots
    if (i == 1){
      z = 3.0 / (1.0 + 2.4 * n);
    } else if (i == 2){
      z = z + 15.0 / (1.0 +2.5 * n);
    } else {
      ai = i - 2;
      z = z + (1.0 + 2.55 * ai) / (1.9 * ai) *(z - x[ ai - 1 ]);
    }

    for (int its = 1; its <= maxiter; its ++){
      p1 = 1.0;
      p2 = 0.0;

      for (size_t j = 1; j <= n; j ++){
	p3 = p2;
	p2 = p1;
	p1 = ((2.0 * j - 1 - z) * p2 - (j - 1) * p3) / j;
      }
      pp = n *(p1 - p2) / z;
      z1 = z;
      z  = z1 - p1 / pp;

      if (std::fabs(z - z1) <= epsilon) break;
    }
    x[ i - 1 ] = z;
    w[ i - 1 ] = -1.0 / (pp * n * p2);
  }
}

void GaussLegendre(double x1, double x2, uint n, RVector & x, RVector & w){
//  taken from matlab-code (Thomas Guenther)
//  function [x, w] = gauleg(x1, x2, n)
//  Given the lower and upper limits of integration x1 and x2 and given n,
//  this routine returns arrays x(1..n) and w(1..n) of length n,
//  containing the abscissas and weights of the Gauss-Legendre
//  n-point quadrature formula.
//  For a description of the following routines see
//  Numerical Recipes, Press et al.

  if (x.size() != n) x.resize(n);
  if (w.size() != n) w.resize(n);

  double epsilon = 3.0e-6;

  double m = (n + 1.0) / 2.0 ;
  double xm = 0.5 * (x2 + x1);
  double xl = 0.5 * (x2 - x1);

  double z = 0.0, z1 = 0.0, p1 = 0.0, p2 = 0.0, p3 = 0.0, pp = 0.0;

  for (int i = 1; i <= m; i ++){
    z = std::cos(PI * (i - 0.25) / (n + 0.5));

   // Starting with the above approximation to the ith root, we enter
   // the main loop of refinements by Newton's method
   z1 = z + 2.0 * epsilon;

   while (std::fabs(z - z1) > epsilon){
     p1 = 1.0;
     p2 = 0.0;
     for (size_t j = 1; j <= n; j ++){
       p3 = p2;
       p2 = p1;
       p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / (double)j;
     }

     // p1 is now the desired Legendre polynomial. We next compute pp,
     // its derivative, by a standard relation involving also p2, the
     // polynomial of one lower order
     pp = (double)n * (z * p1 - p2) / (z * z - 1.0);
     z1 = z;
     z = z1 - p1 / pp; // Newtons method
   }
   // Scale the root to the desired interval, and put in its
   // symmetric counterpart
   x[ i - 1 ] = xm - xl * z;
   x[ n - i ] = xm + xl * z;
   //   x[ n +1 - i ] = xm + xl * z;
   //Compute the weight and ist symmetric counterpart
   w[ i - 1 ] = 2.0 * xl / ((1.0 -z * z) *pp *pp);
   w[ n - i ] = w[ i - 1 ];
   //   w[ n + 1 - i ] = w[ i - 1 ];
  }
}

RVector3 sphTangential2Initerial(const RVector3 &V, double lat, double lon){

    double th = lat * PI/180.0;
    double ph = lon * PI/180.0;
    double ct = std::cos(th);
    double st = std::sin(th);
    double cp = std::cos(ph);
    double sp = std::sin(ph);
    return RVector3(
            (V[0]* ct + V[1] * -st) * cp + V[2] * -sp,
            (V[0]* ct + V[1] * -st) * sp + V[2] *  cp,
             V[0]* st + V[1] * ct);
}

void lineIntegralZ_WonBevis(const RVector3 &p1, const RVector3 &p2,
                            RVector3 &dg, RVector3 &dgz){

    double x1 = p1[0];
    double z1 = p1[1];
    double x2 = p2[0];
    double z2 = p2[1];

    double x21 = x2 - x1;
    double z21 = z2 - z1;
    double z21s = z21 * z21;
    double x21s = x21 * x21;

    double xz12 = x1 * z2 - x2 * z1;

    if (abs(x1) < TOLERANCE && abs(z1) < TOLERANCE){
        dg *= 0.0;
        dgz *= 0.0;
        return;
    }
    if (abs(x2) < TOLERANCE && abs(z2) < TOLERANCE){
        dg *= 0.0;
        dgz *= 0.0;
        return;
    }

    double theta1 = std::atan2(z1, x1);
    double theta2 = std::atan2(z2, x2);

    double r1s = x1 * x1 + z1 * z1;
    double r2s = x2 * x2 + z2 * z2;
    double r1 = std::sqrt(r1s);
    double r2 = std::sqrt(r2s);

    double r21s = x21s + z21s;
    double R2 = r21s;

    double rln = std::log(r2 / r1);

    double p = (xz12/r21s) * ((x1*x21 - z1*z21)/r1s - (x2*x21 - z2*z21)/r2s);
    double q = (xz12/r21s) * ((x1*z21 + z1*x21)/r1s - (x2*z21 + z2*x21)/r2s);

    double Fz = 0.0;
    double Fx = 0.0;
    double Fzx = 0.0; // dFz/dx
    double Fzz = 0.0; // dFz/dz

    if (sign(z1) != sign(z2)){
        if ((x1*z2 < x2*z1) && z2 >=0.0) theta1 = theta1 + PI2;
        if ((x1*z2 > x2*z1) && z1 >=0.0) theta2 = theta2 + PI2;
    }
    if (abs(x1*z2 - x2*z1) < TOLERANCE){
        dg *= 0.0;
        dgz *= 0.0;
        return;
    }
    double th12 = (theta1 - theta2);

    if (abs(x21) < TOLERANCE){
//         #print "case 3"
//         #case 3
        Fz = x1 * rln;
        Fx = 0.0;
        Fzz = -p;
        Fzx = q - z21s/r21s * rln;

//__MS(Fz << " " << Fx << " " <<  R2 << " " <<  x1 << " " <<  z1 << " " <<  x2 << " " <<  z2)
    } else { //: #default
        double B = z21 / x21;
        double A = (x21 * xz12) / R2;

        Fz = A * (th12 + B * rln);
        Fx = A * (-th12 * B + rln);
        double z21dx21 = z21 / x21;
        double z21x21 = z21 * x21;

        double fz = (th12 + z21dx21 * rln)/r21s;

        Fzz = -p + x21s * fz;
        Fzx = q - z21x21 * fz;

//__MS(Fz << " " << Fx << " " <<  R2 << " " <<  x1 << " " <<  z1 << " " <<  x2 << " " <<  z2)
//         #// check this
//         #fx = (th12 * z21dx21 - rln)/r21s
    }
    dg[0] = Fx;
    dg[1] = 0.0;
    dg[2] = Fz;
    dgz[0] = Fzx;
    dgz[1] = 0.0;
    dgz[2] = Fzz;
    return;
}

double lineIntegralZ_WonBevis(const RVector3 &p1, const RVector3 &p2){
    double dgz = 0.0;

    double x1 = p1[0];
    double z1 = p1[1];
    double x2 = p2[0];
    double z2 = p2[1];

    double x21 = x2 - x1;
    double z21 = z2 - z1;
    double z21s = z21 * z21;
    double x21s = x21 * x21;

    double xz12 = x1 * z2 - x2 * z1;

    if (abs(x1) < TOLERANCE && abs(z1) < TOLERANCE){
        return 0.0;
    }
    if (abs(x2) < TOLERANCE && abs(z2) < TOLERANCE){
        return 0.0;
    }

    double theta1 = std::atan2(z1, x1);
    double theta2 = std::atan2(z2, x2);

    double r1s = x1 * x1 + z1 * z1;
    double r2s = x2 * x2 + z2 * z2;
    double r1 = std::sqrt(r1s);
    double r2 = std::sqrt(r2s);

    double r21s = x21s + z21s;
    double R2 = r21s;
    double rln = std::log(r2 / r1);

    if (sign(z1) != sign(z2)){
        if ((x1*z2 < x2*z1) && z2 >=0.0) theta1 = theta1 + PI2;
        if ((x1*z2 > x2*z1) && z1 >=0.0) theta2 = theta2 + PI2;
    }
    if (abs(x1*z2 - x2*z1) < TOLERANCE){
        return 0.0;
    }

    if (abs(x21) < TOLERANCE){
        dgz = x1 * rln;
    } else { //: #default
        double B = z21 / x21;
        double A = (x21 * xz12) / R2;
        double th12 = (theta1 - theta2);
        dgz = A * (th12 + B * rln);
    }
    return dgz;
}



} //namespace GIMLI;
