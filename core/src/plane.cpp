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

#include "plane.h"
#include "line.h"
#include "pos.h"
#include "matrix.h"
//#include "numericbase.h"

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const Plane & p){
  if (p.valid()){
    str << "Plane: norm = " << p.norm() << " p = " << p.d() << " ";
  } else {
    str << "Plane: invalid ";
  }
  return str;
}

Plane::Plane()
    : valid_(false) {
}

Plane::Plane(const RVector3 & norm, double d)
    : norm_(norm), d_(d), valid_(false){
    checkValidity();
}

Plane::Plane(const RVector3 & norm, const RVector3 & x0)
    : norm_(norm), valid_(false){
    d_ = x0.abs();
    checkValidity();
}

Plane::Plane(const RVector3 & p0, const RVector3 & p1, const RVector3 & p2)
    : valid_(false){

    norm_ = p0.norm(p1, p2);
    if (std::fabs(norm_.abs() - 1.0) > TOLERANCE){
        __MS(p0)
        __MS(p1)
        __MS(p2)
        __MS(norm_)
        __MS(std::fabs(norm_.abs() - 1.0)<< " "<< TOLERANCE)
    }

    RMatrix A(3, 3);
    A[0][0] = p0[0];
    A[0][1] = p0[1];
    A[0][2] = p0[2];
    A[1][0] = p1[0];
    A[1][1] = p1[1];
    A[1][2] = p1[2];
    A[2][0] = p2[0];
    A[2][1] = p2[1];
    A[2][2] = p2[2];

    RMatrix A0(A);
    A0[0][0] = 1.0; A0[1][0] = 1.0; A0[2][0] = 1.0;
    RMatrix A1(A);
    A1[0][1] = 1.0; A1[1][1] = 1.0; A1[2][1] = 1.0;
    RMatrix A2(A);
    A2[0][2] = 1.0; A2[1][2] = 1.0; A2[2][2] = 1.0;

    double detA = det(A);
    if (std::fabs(det(A)) < TOLERANCE){
         d_ = 0.0;
    } else {
         d_ = detA / (std::sqrt(std::pow(det(A0), 2.0) + std::pow(det(A1), 2.0) + std::pow(det(A2), 2.0)));
    }

    d_ = rint(d_ / (TOLERANCE * 1e-2)) * TOLERANCE * 1e-2;
    //norm_.round();
    checkValidity();
}

Plane::Plane(double a, double b, double c, double d)
    : valid_(false){

    double abstmp = std::sqrt(a * a + b * b + c * c);
    norm_ = RVector3(a / abstmp, b / abstmp, c / abstmp);
    d_ = d / abstmp;
    checkValidity();
}

Plane::Plane(const Plane & plane){
    copy_(plane);
}

Plane::~Plane(){
}

Plane & Plane::operator = (const Plane & plane){
    if (this != & plane){
        copy_(plane);
    } return *this;
}

void Plane::copy_(const Plane & plane){
    norm_ = plane.norm();
    d_ = plane.d();
    valid_ = plane.valid();
}

bool Plane::checkValidity(double tol){
    if (std::fabs(norm_.abs() - 1.0) < tol){
        valid_ = true;
    } else {
        std::cerr << WHERE_AM_I << " WARNING! Plane NOT valid " 
                  << std::fabs(norm_.abs() - 1.0) << " / " << tol << std::endl;
        valid_ = false;
    }
    return false;
}

bool Plane::touch(const RVector3 & pos, double tol){
    if (valid_){
        return (std::fabs(this->distance(pos)) < tol);
    }
    return false;
}

RVector3 Plane::intersect(const Line & line, double tol, bool inside){
    //** http://astronomy.swin.edu.au/~pbourke/geometry/planeline/;
    if (!this->valid() || !line.valid()) return RVector3(false);

    double n = norm_.dot(line.p1() - line.p0());

    //! Line and Plane are parallel
    if (::fabs(n) < TOLERANCE) return RVector3(false);
    if (this->touch(line.p0(), tol) && this->touch(line.p1(), tol)) return RVector3(false);

    RVector3 x0(this->x0());

    double t = norm_.dot(x0 - line.p0()) / n;
   //  double t = norm.scalar(p3 - p1) / norm.scalar(p2 - p1);

    if (isnan(t) || isinf(t)) {
        std::cout << norm_.dot(x0 - line.p0()) << " " << norm_.dot(line.p1() - line.p0()) << std::endl;
        throwError(WHERE_AM_I + " t is not valid " + str(t));
        return RVector3(false);
    }

    if (inside && (t < -TOLERANCE || t > 1.0 + TOLERANCE)) return RVector3(false);

    return line(t);
}

Line Plane::intersect(const Plane & plane, double tol){
    //** both planes are identical return invalid line
    if ((*this) == plane) return Line();

    RVector3 a = norm_.cross(plane.norm());
    //** both planes are parallel return invalid line;
    if (a.abs() < tol) return Line();

    //p = c1 N1 + c2 N2 + u N1 * N2 // line definition;
    double n0n0 = norm_.dot(norm_);
    double n0n1 = norm_.dot(plane.norm());
    double n1n1 = plane.norm().dot(plane.norm());
    double det = n0n0 * n1n1 - n0n1 * n0n1;
    double d0 = this->d();
    double d1 = plane.d();

    double c0 =  (d0 * n1n1 - d1 * n0n1) / det;
    double c1 =  (d1 * n0n0 - d0 * n0n1) / det;

    RVector3 lineP0(norm_ * c0 + plane.norm() * c1);
    RVector3 lineP1(a + lineP0);

    return Line(lineP0, lineP1);
}

bool Plane::compare(const Plane & plane, double tol, bool bothDirs){
    if (bothDirs == true){
        // __MS(std::fabs(std::fabs(d_) - std::fabs(plane.d())))
        // __MS(norm_.distance(plane.norm()))
        // __MS(std::fabs(norm_.distance(plane.norm()) - 2.0))

        if (std::fabs(std::fabs(d_) - std::fabs(plane.d())) > tol) return false;
        
        if (norm_.distance(plane.norm()) < tol || 
             std::fabs(norm_.distance(plane.norm()) - 2.0) < tol) return true;

    } else {
        if (std::fabs(d_ - plane.d()) > tol) return false;
        if (norm_.distance(plane.norm()) < tol) return true;
    } 

    return false;
}

} //namespace GIMLI
