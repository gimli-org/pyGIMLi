/******************************************************************************
 *   Copyright (C) 2007-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_PLANE__H
#define _GIMLI_PLANE__H

#include "gimli.h"
#include "pos.h"

namespace GIMLI{

//! A plane
/*! A plane, defined through Hessian normal form, norm * x = -p.
    Stores Normvector norm with |norm| = 1.0 and distance d from orign. */
class DLLEXPORT Plane{
public:
    /*! Default constructor, Construct an invalid empty plane. */
    Plane();

    /*! Construct a plane based on Hessian normal form norm * x = -p.*/
    Plane(const RVector3 & norm, double d);

    /*! Construct a plane based on his unit nomal vector norm and a base position x0.*/
    Plane(const RVector3 & norm, const RVector3 & x0);

    /*! Construct a plane based on 3 real positions in R^3, respectivly parameterized style.*/
    Plane(const RVector3 & p0, const RVector3 & p1, const RVector3 & p2);

    /*! Construct a plane based on his general equation */
    Plane(double a, double b, double c, double d);

    /*! Copyconstructor.*/
    Plane(const Plane & plane);

    /*! Default destructor*/
    ~Plane();

    /*! Assignment operator*/
    Plane & operator = (const Plane & plane);

    /*! Equal_to operator */
    inline bool operator == (const Plane & plane){ 
        return this->compare(plane);
    }

    /*! Not_equal_to operator */
    inline bool operator != (const Plane & plane){ 
        return !(*this== plane);
    }

    /*! Compare two planes with a given tolerance. Check if both norms and distances are equal.
    | norm - p.norm | < tol && | d_ - p.d | < tol. 
    If bothDirection is True  | norm - p.norm | can be 2.0 to allow for coplanar check with different orientations.*/
    bool compare(const Plane & p, double tol=TOLERANCE, bool bothDirs=false);
    

    /*! Returns true if the plane is valid and pos touch this plane.
        Touch when plane.distance(pos) < tol. */
    bool touch (const RVector3 & pos, double tol = TOLERANCE);

    /*! Returns the \ref Line of intersection between 2 planes.
        Are booth planes parallel or identically the returned line is invalid. */
    Line intersect(const Plane & plane, double tol = TOLERANCE);

    /*! Returns the point of intersection between this plane and the \ref Line line. Return an
        invalid RVector3 if line and this plane are parallel. Optional inside check.
        Set inside returns invalid RVector3 if intersection point is not inside (including nodes) the line */
    RVector3 intersect(const Line & line, double tol = TOLERANCE, bool inside = false);

    /*! Return a const reference to the unit vector of this plane.
        | norm | = 1.0 */
    inline const RVector3 & norm() const { return norm_; }

    /*! Returns the orign if the unit vector. x0 = n * d. */
    inline const RVector3 x0() const { return norm_ * d_; }

    /*! Returns the distance between the Point pos and this plane.*/
    inline double distance(const RVector3 & pos) const { return norm_.dot(pos) - d_ ; }

    /*! Distance between the orign(0,0,0) and this plane.*/
    inline double d() const { return d_; }

    /*! Check if the plane spans norm vector has a length of 1.0 that means it spans a valid R3 space. */
    bool checkValidity(double tol = TOLERANCE);

    /*! Return the validity of this plane. This plane is not valid if its initialized by default constructor.*/
    inline bool valid() const { return valid_; }

protected:

    /*! internal: copy content of plane into this plane */
    void copy_(const Plane & plane);

    RVector3 norm_;
    double d_;
    bool valid_;
};

std::ostream & operator << (std::ostream & str, const Plane & p);

} //namespace GIMLI

#endif //  _GIMLI_PLANE__H
