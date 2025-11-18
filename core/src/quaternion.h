/******************************************************************************
 *   Copyright (C) 2008-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_QUATERNION__H
#define _GIMLI_QUATERNION__H

#include <iostream>
#include "pos.h"

namespace GIMLI{

template < class ValueType > class Quaternion;

typedef Quaternion< double > RQuaternion;

DLLEXPORT RMatrix getRotation(const Pos & src, const Pos & dest);

template < class ValueType > class DLLEXPORT Quaternion{
public:
    Quaternion(ValueType w = 1.0, ValueType i = 0.0, ValueType j = 0.0, ValueType k = 0.0)
        : re_(w){
        im_[0] = i; im_[1] = j; im_[2] = k;
    }

    Quaternion(ValueType re, const Pos & im)
        : re_(re), im_(im){ }

    Quaternion(const Pos & xAxis, const Pos & yAxis,
                const Pos & zAxis){
        THROW_TO_IMPL
    }

    Quaternion(const Quaternion < ValueType > & q) { copy_(q); }

    Quaternion & operator = (const Quaternion < ValueType > & q) {
        if (this != & q){
            copy_(q);
        }
        return *this;
    }

    const ValueType operator [] (const size_t i) const {
        if (i == 0) return re_; else return im_[i - 1];
    }

    ValueType & operator [] (const size_t i){ if (i == 0) return re_; else return im_[i - 1]; }

    Quaternion & operator /= (const ValueType & s) { re_ /= s; im_ /= s; return * this;}
    Quaternion & operator *= (const ValueType & s) { re_ *= s; im_ /= s; return * this;}

    void createFromAxisAngle(const Pos & axis, double angle){
        //** The quaternion representing the rotation is
        //**   q = cos(a/2)+sin(a/2)*(x*i+y*j+z*k)
        double ah = 0.5 * angle;
        re_ = std::cos(ah);
        im_ = axis * std::sin(ah);
    }

    template < class Matrix > void rotMatrix(Matrix & rot) const {
        ValueType x  = 2.0 * im_[0], y  = 2.0 * im_[1], z  = 2.0 * im_[2];

        ValueType wx = x * re_, wy = y * re_, wz = z * re_;
        ValueType xx = x * im_[0], xy = y * im_[0], xz = z * im_[0];
        ValueType yy = y * im_[1], yz = z * im_[1], zz = z * im_[2];

        rot[0][0] = 1.0 - (yy + zz);
        rot[0][1] = xy - wz;
        rot[0][2] = xz + wy;
        rot[1][0] = xy + wz;
        rot[1][1] = 1.0 - (xx + zz);
        rot[1][2] = yz - wx;
        rot[2][0] = xz - wy;
        rot[2][1] = yz + wx;
        rot[2][2] = 1.0 - (xx + yy);
    }

    inline ValueType norm() const { return re_ * re_ + im_.distSquared(); }

    inline ValueType length() const { return std::sqrt(norm()); }

    inline void normalise(){ *this /= length(); }

    inline void setRe(const ValueType re){ re_ = re; }
    inline ValueType re() const { return re_; }

    inline void setIm(const Pos & im) {  im_ = im; }
    inline Pos im() const { return im_; }

protected:

    void copy_(const Quaternion < ValueType > & q){ re_ = q.re(); im_ = q.im(); }

    ValueType         re_;
    Pos im_;
};

DLLEXPORT std::ostream & operator << (std::ostream & str, const RQuaternion & q);

} // namespace GIMLI
#endif
