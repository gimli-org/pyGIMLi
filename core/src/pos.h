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

#ifndef _GIMLI_POS__H
#define _GIMLI_POS__H

#include "gimli.h"
#include "vector.h"

namespace GIMLI{

std::ostream & operator << (std::ostream & str, const Pos & pos);
std::istream & operator >> (std::istream & is, Pos & pos);


DLLEXPORT RVector3 center(const R3Vector & vPos);
DLLEXPORT R3Vector normalise(const R3Vector & vPos);

DLLEXPORT double jacobianDetXY(const RVector3 & p1, const RVector3 & p2, const RVector3 & p3);
DLLEXPORT double angle(const RVector3 & p1, const RVector3 & p2, const RVector3 & p3);

DLLEXPORT bool xVari(const R3Vector & electrodeList);
DLLEXPORT bool yVari(const R3Vector & electrodeList);
DLLEXPORT bool zVari(const R3Vector & electrodeList);

/*! Return array of all x-coordinates. [:,0]*/
DLLEXPORT RVector x(const R3Vector & rv);
/*! Return array of all y-coordinates. [:,1]*/
DLLEXPORT RVector y(const R3Vector & rv);
/*! Return array of all z-coordinates. [:,2]*/
DLLEXPORT RVector z(const R3Vector & rv);

DLLEXPORT RVector absR3(const R3Vector & vPos);

DLLEXPORT void swapXY(R3Vector & rv);
DLLEXPORT void swapXZ(R3Vector & rv);
DLLEXPORT void swapYZ(R3Vector & rv);

DLLEXPORT std::vector < RVector3 > loadRVector3(const std::string & fileName);
DLLEXPORT void saveRVector3(const std::vector < RVector3 > l, const std::string & fileName);

/*! Create one dimensional array from R3Vector
 * return = [vev[0][0], vev[0][1], vev[0][2], vev[1][0] .. ] */
DLLEXPORT RVector toArray(const R3Vector & vec);

/*! Create two dimensional [n x 3] array from R3Vector */
DLLEXPORT RMatrix toMatrix(const R3Vector & vec);

/*! Temporary transformation vor R3Vector until std::vector < RVector3 > will be removed. */
DLLEXPORT R3Vector stdVectorRVector3ToR3Vector(const std::vector < RVector3 > & rv);

/*! Temporary transformation vor R3Vector until std::vector < RVector3 > will be removed. */
DLLEXPORT std::vector < RVector3 > R3VectorTostdVectorRVector3(const R3Vector & rv);

//! 3 dimensional vector
/*! 3 dimensional vector */

class DLLEXPORT Pos {
public:

  //  static const RVector3 ZERO(0.0, 0.0, 0.0);
//   static const RVector3 UNIT_X;
//   static const RVector3 UNIT_Y;
//   static const RVector3 UNIT_Z;
//   static const RVector3 NEGATIVE_UNIT_X;
//   static const RVector3 NEGATIVE_UNIT_Y;
//   static const RVector3 NEGATIVE_UNIT_Z;
//   static const RVector3 UNIT_SCALE;

    /*! Construct an empty in 3d at (0, 0, 0) */
    Pos() : valid_(true) { assign(0.0, 0.0, 0.0); }

    /*! Construct an empty in 3d at (0, 0, 0). Optional set valid flag. */
    Pos(bool valid) : valid_(valid) { assign(0.0, 0.0, 0.0); }

    Pos(double x, double y) : valid_(true) { assign(x, y, 0.0); }
    Pos(double x, double y, double z) : valid_(true) { assign(x, y, z); }

    Pos(const Pos & pos) { copy_(pos); }

    /*! Assignment operator */
    Pos & operator = (const Pos & pos){
        if (this != & pos){ copy_(pos); } return *this; }

    /*! Assignment operator */
    Pos & operator = (const Vector < double > & v){
        if (v.size() > 2) {
            mat_[0] = v[0];
            mat_[1] = v[1];
            mat_[2] = v[2];
        } else {
            throwLengthError(WHERE_AM_I + " v.size() < 2 " + str(v.size()));
        }
        return *this;

    }

    inline double & operator [] (Index i) { return mat_[i]; }

    inline const double & operator [] (Index i) const { return mat_[i]; }

    #define DEFINE_UNARY_MOD_OPERATOR__(OP) \
    inline Pos & operator OP##= (const double & b){ mat_[0] OP##= b; mat_[1] OP##= b; mat_[2] OP##= b; return *this; }\
    inline Pos & operator OP##= (const Pos & b){ mat_[0] OP##= b[0]; mat_[1] OP##= b[1]; mat_[2] OP##= b[2]; return *this; }\

    DEFINE_UNARY_MOD_OPERATOR__(+)
    DEFINE_UNARY_MOD_OPERATOR__(-)
    DEFINE_UNARY_MOD_OPERATOR__(/)
    DEFINE_UNARY_MOD_OPERATOR__(*)

    #undef DEFINE_UNARY_MOD_OPERATOR__

    Pos operator - () const { return Pos(-mat_[0], -mat_[1], -mat_[2]); }

    inline void setValid(bool valid) { valid_ = valid; }
    inline bool valid() const { return valid_; }

    inline void assign(const double & x, const double & y, const double & z) {
        mat_[0] = x; mat_[1] = y; mat_[2] = z;
      //  x_ = x; y_ = y; z_ = z;
    }

    inline const double & x() const { return mat_[0]; }
    inline const double & y() const { return mat_[1]; }
    inline const double & z() const { return mat_[2]; }
    inline void setX(double x) { mat_[0] = x; }
    inline void setY(double y) { mat_[1] = y; }
    inline void setZ(double z) { mat_[2] = z; }

    /*! Set a value. Throws out of range exception if index check fails. */
    inline void setVal(const double & val, Index i) {
        if (i < 3) {
            mat_[i] = val;
        } else {
            throwRangeError(WHERE_AM_I, i, 0, 3);
        }
    }

    /*! Get a value. Throws out of range exception if index check fails. */
    inline const double & getVal(Index i) const {
        if (i < 3) {
            return mat_[i];
        } else {
            throwRangeError(WHERE_AM_I, i, 0, 3);
        }
        return mat_[0];
    }

    /*! Inline round to tol.*/
    inline Pos & round(double tol){
        mat_[0] = rint(mat_[0] / tol) * tol;
        mat_[1] = rint(mat_[1] / tol) * tol;
        mat_[2] = rint(mat_[2] / tol) * tol;
        return *this;
    }

    /*! Return the squared distance to p.*/
    inline double distSquared(const Pos & p) const {
        return  ((mat_[0] - p[0]) * (mat_[0] - p[0]) +
                   (mat_[1] - p[1]) * (mat_[1] - p[1]) +
                   (mat_[2] - p[2]) * (mat_[2] - p[2]));
    }

    /*! Return the distance to p.*/
    inline double dist(const Pos & p) const { return std::sqrt(distSquared(p)); }

    /*! Return the distance to p.*/
    inline double distance(const Pos & p) const { return dist(p); }

    /*! Return the squared length of this position vector.*/
    inline double distSquared() const {
        return  mat_[0] * mat_[0] + mat_[1] * mat_[1] + mat_[2] * mat_[2];
    }
    /*! Return the length of this position vector, same as length.*/
    inline double abs() const { return length(); }

    /*! Return the length of this position vector.*/
    inline double length() const { return std::sqrt(distSquared()); }

    /*! Return the angle between (this, (origin), p).*/
    double angle(const Pos & p) const;

    /*! Return the angle between (p1, this, p2).*/
    double angle(const RVector3 & p1, const RVector3 & p3) const;

    inline double dot(const Pos & p) const {
        return mat_[0] * p[0] + mat_[1] * p[1] + mat_[2] * p[2];
    }

    inline double sum() const {
        return mat_[0] + mat_[1] + mat_[2];
    }

    Pos norm(const Pos & p1, const Pos & p2) const;

    /*! Return normalised copy of this Pos. */
    Pos norm() const {
        Pos p(*this);
        p.normalize();
        return p;
    }

    /*! Normalize this Pos and return itself. */
    Pos & normalize(){
        double t = this->abs();
        if (t > TOLERANCE) *this /= t;
        return *this;
    }

    /*!DEPRECATED Normalise for backward compatibility.*/
    Pos & normalise(){ return normalize(); }

    Pos cross(const Pos & p) const;

    Pos normXY(const Pos & p) const;

    template < class Matrix > Pos & transform(const Matrix & wm){
        double x = mat_[0], y = mat_[1], z = mat_[2];

        mat_[0] = x * wm[0][0] + y * wm[0][1] + z * wm[0][2];
        mat_[1] = x * wm[1][0] + y * wm[1][1] + z * wm[1][2];
        mat_[2] = x * wm[2][0] + y * wm[2][1] + z * wm[2][2];
        return *this;
    }

    inline Pos & rotateX(double phi){
        double mat[3][3] ={{ 1.0,           0.0,            0.0},
                           { 0.0, std::cos(phi), -std::sin(phi)},
                           { 0.0, std::sin(phi),  std::cos(phi)} };
        return this->transform(mat);
    }
    inline Pos & rotateY(double phi){
        double mat[3][3] =  {{std::cos(phi),  0.0, std::sin(phi)},
                             {0.0,            1.0,           0.0},
                             {-std::sin(phi), 0.0, std::cos(phi)}};

        return this->transform(mat);
    }
    inline Pos & rotateZ(double phi){
        double mat[3][3] = {{std::cos(phi), -std::sin(phi), 0.0},
                            {std::sin(phi),  std::cos(phi), 0.0},
                            {          0.0,            0.0, 1.0}};
        return this->transform(mat);
    }

    inline Pos & rotate(const RVector3 & r){
        return this->rotateX(r[0]).rotateY(r[1]).rotateZ(r[2]);
    }
    inline Pos & rotate(double phiX, double phiY, double phiZ){
        return this->rotateX(phiX).rotateY(phiY).rotateZ(phiZ);
    }

    inline Pos & scale(const RVector3 & s){ return (*this) *= s;}

    inline Pos & translate(const RVector3 & t){ return (*this) += t;}

    inline Pos & translate(double x, double y=0.0, double z=0.0){
        mat_[0] += x; mat_[1] += y; mat_[2] += z; return *this; }


    RVector vec() const {
        RVector tmp(3);
        tmp[0] = mat_[0];
        tmp[1] = mat_[1];
        tmp[2] = mat_[2];
        return tmp;
    }

    Index hash() const {
        return GIMLI::hash(mat_[0], mat_[1], mat_[2], valid_);
    }

    /*!Swap coordinates i and j [0..3) */
    void swap(Index i, Index j){
        ASSERT_RANGE(i, 0, 3);
        ASSERT_RANGE(j, 0, 3);
        if (i != j){
            double tmp = mat_[i];
            mat_[i] =  mat_[j];
            mat_[j] =  tmp;
        }
    }
protected:

    inline void copy_(const Pos & pos) {
        valid_ = pos.valid(); assign(pos[0], pos[1], pos[2]);
    }

    bool valid_;

    double mat_[3];

};

// template < class double > const RVector3 RVector3::ZERO(0.0, 0.0, 0.0);
// template < class double > const RVector3 RVector3::UNIT_X(1.0, 0.0, 0.0);
// template < class double > const RVector3 RVector3::UNIT_Y(0.0, 1.0, 0.0);
// template < class double > const RVector3 RVector3::UNIT_Z(0.0, 0.0, 1.0);
// template < class double > const RVector3 RVector3::NEGATIVE_UNIT_X(-1.0,  0.0,  0.0);
// template < class double > const RVector3 RVector3::NEGATIVE_UNIT_Y(0.0, -1.0,  0.0);
// template < class double > const RVector3 RVector3::NEGATIVE_UNIT_Z(0.0,  0.0, -1.0);
// template < class double > const RVector3 RVector3::UNIT_SCALE(1.0, 1.0, 1.0);

inline bool operator == (const RVector3 & a , const RVector3 & b){
    if (a.valid() != b.valid()) return false;
    if (a.distSquared(b) < TOLERANCE) return true; else return false;
}
inline bool operator != (const RVector3 & a , const RVector3 & b){
    return !(a == b);
}
inline bool operator < (const RVector3 & a , const RVector3 & b){
    return a.distSquared() < b.distSquared();
}
inline bool operator <= (const RVector3 & a , const RVector3 & b){
    return a.distSquared() <= b.distSquared();
}
inline bool operator > (const RVector3 & a , const RVector3 & b){
    return a.distSquared() > b.distSquared();
}
inline bool operator >= (const RVector3 & a , const RVector3 & b){
    return a.distSquared() >= b.distSquared();
}
inline RVector3 RINT(const RVector3 & a) {
    std::cout << WHERE_AM_I << std::endl;
    return RVector3();
}

#define DEFINE_POS_BIN_OPERATOR__(OP)                      \
inline Pos operator OP (const Pos & a, const Pos & b){ \
    Pos tmp(a); return tmp OP##= b; } \
inline Pos operator OP (const Pos & a, const double & b){ \
    Pos tmp(a); return tmp OP##= b; } \
inline Pos operator OP (const double & a, const Pos & b){ \
    Pos tmp(a, a, a); return tmp OP##= b; } \
\

DEFINE_POS_BIN_OPERATOR__(+)
DEFINE_POS_BIN_OPERATOR__(*)
DEFINE_POS_BIN_OPERATOR__(/)
DEFINE_POS_BIN_OPERATOR__(-)


inline std::ostream & operator << (std::ostream & str, const Pos & pos){
  if (pos.valid()){
    str << pos[0] << "\t" << pos[1] << "\t" << pos[2];
  } else {
    str << " pos is not valid";
  }
  return str;
}

inline std::istream & operator >> (std::istream & is, Pos & pos){
    THROW_TO_IMPL
    return is;
}

/*! Sort increasing x and decreasing y*/
inline bool posLesserXrY(const RVector3 & a, const RVector3 & b){
    if (abs(a[0] - b[0]) < TOLERANCE) {
        if (abs(a[1] - b[1]) < TOLERANCE) {
            return a[2] < b[2];
        } else {
            return a[1] > b[1];
        }
    } else return a[0] < b[0];
}

/*! Sort increasing x and decreasing y*/
inline bool posLesserXYrZ(const RVector3 & a, const RVector3 & b){
    if (abs(a[0] - b[0]) < TOLERANCE) {
        if (abs(a[1] - b[1]) < TOLERANCE) {
            return a[2] > b[2];
        } else {
            return a[1] < b[1];
        }
    } else return a[0] < b[0];
}

inline bool posLesserX(const RVector3 & a, const RVector3 & b){
    if (abs(a[0] - b[0]) < TOLERANCE) {
        if (abs(a[1] - b[1]) < TOLERANCE) {
            return a[2] < b[2];
        } else {
            return a[1] < b[1];
        }
    } else return a[0] < b[0];
}

} // namespace GIMLI;

#ifndef PYGIMLI_CAST
namespace std {
    template<> struct hash<GIMLI::Pos> {
        GIMLI::Index operator()(const GIMLI::Pos & p) const noexcept {
            return p.hash();
        }
    };
    template<> struct hash< GIMLI::PosVector > {
        GIMLI::Index operator()(const GIMLI::PosVector & p) const noexcept {
            return p.hash();
        }
    };
}
#endif // PYGIMLI_CAST
#endif // _GIMLI_POS__H
