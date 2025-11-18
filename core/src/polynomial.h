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

#ifndef _GIMLI_POLYNOMIAL__H
#define _GIMLI_POLYNOMIAL__H

#include "gimli.h"
#include "matrix.h"
#include "modellingbase.h"
#include "regionManager.h"

#include "numericbase.h"

#include <list>

namespace GIMLI{

template< class ValueType > class PolynomialElement {
public:
    PolynomialElement(Index i, Index j, Index k, const ValueType & val)
    :   i_(i), j_(j), k_(k), val_(val){
    }

    inline ValueType operator () (const Pos & xyz) const {
        return val_ * powInt(xyz[0], i_) * powInt(xyz[1], j_) * powInt(xyz[2], k_);
    }

    Index i_, j_, k_;
    ValueType val_;
};

/*! waste to satisfy python bindings */
template < class ValueType > bool operator < (const PolynomialElement < ValueType > & a, 
                                              const PolynomialElement < ValueType > & b){
    return ((a.i_ < b.i_) && (a.j_ < b.j_) && (a.k_ < b.k_));
}

template < class ValueType > bool operator != (const PolynomialElement < ValueType > & a, 
                                               const PolynomialElement < ValueType > & b){
    return !(a == b);
}

template < class ValueType > bool operator == (const PolynomialElement < ValueType > & a, 
                                               const PolynomialElement < ValueType > & b){
    return ((a.i_ == b.i_) && (a.j_ == b.j_) && (a.k_ == b.k_) && (a.val_ == b.val_));
}


/*! Three dimensional polynomial function. For symbolic calculation and derivation.
 I.e. f(x,y,z) = 1 + x + y + z + xy + xz + zz ...
*/
template< class ValueType > class DLLEXPORT PolynomialFunction {
public:

    /*! Create empty polynomial */
    PolynomialFunction(uint size = 0){
        init_(Vector< ValueType >(size, 0.0), 
              Vector< ValueType >(size, 0.0), 
              Vector< ValueType >(size, 0.0));
    }

    /*! Create a polynomial \f$ f(x,y,z) \f$ from one dimensional polynomial coefficient array ax:
        \f$ f(x,y,z) = ax[0] + ax[1]x + ax[2]x^2 ... \f$
     */
    PolynomialFunction(const Vector < ValueType > & ax){
        init_(ax, Vector< ValueType >(0, 0.0), 
                  Vector< ValueType >(0, 0.0));
    }

    /*! Create a polynomial function \f$ f(x,y,z) \f$ from one dimensional polynomial coefficient array ax and ay:
        \f$ f(x,y,z) = ax[0] + ax[1]x + ax[2]x^2 + ... + ay[0] + ay[1]y + ay[2]y^2 + ... \f$.
        If u need a composed polynomial use f(x,y,z)*g(x,y,z).
     */
    PolynomialFunction(const Vector < ValueType > & ax, const Vector < ValueType > & ay){
        init_(ax, ay, Vector< ValueType >(0, 0.0));
    }

    /*! Create a polynomial function \f$ f(x,y,z) \f$ from one dimensional polynomial coefficient array ax, ay, az:
        \f$ f(x,y,z) = ax[0] + ax[1]x + ax[2]x^2 + ... + ay[0] + ay[1]y + ay[2]y^2 + ... + az[0] + az[1]z + az[2]z^2 + ... \f$
        If u need a composed polynomial use f(x,y,z)*g(x,y,z)*h(x,y,z).
     */
    PolynomialFunction(const Vector < ValueType > & ax, 
                       const Vector < ValueType > & ay, 
                       const Vector < ValueType > & az){
        init_(ax, ay, az);
    }

    /*! Return f(x, y, z[k]) polynomial matrix*/
    RMatrix & operator [] (Index k){ return mat_[k]; }

    /*! Return f(x, y, z[k]) polynomial matrix*/
    const RMatrix & operator [] (Index k) const { return mat_[k]; }

    /*! Return reference to f(x, y, z[k]) polynomial matrix. For python only*/
    RMatrix & matR(Index k){ return mat_[k]; }

    /*! Evaluate f(x,y,z) */
    ValueType operator () (const Pos & xyz) const {
        ValueType ret = 0.0;

        for (typename std::vector <PolynomialElement < ValueType > >::const_iterator it = elementList_.begin();
                it != elementList_.end(); it ++){
                ret += (*it)(xyz);
        }
        return ret;
    }

    /*! Evaluate f(x_i,y_i,z_i) for i = [0 .. xyz.size()). */
    Vector < ValueType > operator () (const std::vector < Pos > & xyz) const {
        Vector < ValueType > ret(xyz.size(), 0.0);

        for (Index i = 0 ; i < ret.size(); i ++) ret [i] = (*this)(xyz[i]);

        return ret;
    }

    /*! Return the size of this polynomial. */
    Index size() const { return mat_.size(); }

    /*! Create new polynomial function for \f$ \frac{\partial f(x,y,z)}{\partial \mathrm{dim}}\f$
     @param dim = 0 \f$ \frac{\partial f(x,y,z)}{\partial x}\f$
     @param dim = 1 \f$ \frac{\partial f(x,y,z)}{\partial y}\f$
     @param dim = 2 \f$ \frac{\partial f(x,y,z)}{\partial z}\f$
     */
    PolynomialFunction < ValueType > derive(uint dim) const {

        PolynomialFunction< ValueType > ret(RVector(this->size(), 0.0));

        if (dim == 0){
            for (Index k = 0; k < mat_.size(); k ++){ // z
                for (Index j = 0; j < mat_[k].rows(); j ++){ // y
                    for (Index i = 0; i < mat_[k][j].size() -1; i ++){ // x
                        ret[k][i][j] = mat_[k][i + 1][j] * (i + 1.0);
                    }
                }
            }
        } else if(dim == 1){
            for (Index k = 0; k < mat_.size(); k ++){ // z
                for (Index j = 0; j < mat_[k].rows() - 1; j ++){ // y
                    for (Index i = 0; i < mat_[k][j].size(); i ++){ // x
                        ret[k][i][j] = mat_[k][i][j + 1] * (j + 1.0);
                    }
                }
            }
        } else if (dim == 2){
            for (Index k = 0; k < mat_.size() -1 ; k ++){ // z
                for (Index j = 0; j < mat_[k].rows(); j ++){ // y
                    for (Index i = 0; i < mat_[k][j].size(); i ++){ // x
                        ret[k][i][j] = mat_[k + 1][i][j] * (k + 1.0);
                    }
                }
            }
        }
        ret.fillElementList();
        return ret;
    }


    void fillElementList(){
        elementList_.clear();

        for (Index k = 0; k < mat_.size(); k ++){ // z
            for (Index j = 0; j < mat_[k].rows(); j ++){ // y
                for (Index i = 0; i < mat_[k][j].size(); i ++){ // x
                    if (::fabs(mat_[k][i][j]) > TOLERANCE){
                        elementList_.push_back(PolynomialElement< ValueType > (i, j, k, mat_[k][i][j]));
                    }
                }
            }
        }
    }

    const std::vector <PolynomialElement < ValueType > > & elements() const { return elementList_; }

    void clear() {
        elementList_.clear();
        for (Index k = 0; k < mat_.size(); k ++){ mat_[k] *= 0.0; }
    }

    /*!
     * Fill the parameter coefficients from array. If c.size() == this.size()^3.
     * mat_[k][i][j] = c[k*(size() * size())+ j * size() + i]
     * and return the \ref PolynomialFunction itself.
     * If c.size() is size of elementList_, assume that only the values from elementList_ will be exchanged.
     * Please note, all values of c will be snapped to tolerance.
     *
     */
    PolynomialFunction <  ValueType > & fill(const Vector < ValueType > & c){
        if (c.size() == powInt(mat_.size(), 3)) {
            for (Index k = 0; k < mat_.size(); k ++){ // z
                for (Index j = 0; j < mat_[k].rows(); j ++){ // y
                    for (Index i = 0; i < mat_[k][j].size(); i ++){ // x
                        if (::fabs(c[k*(size() * size())+ j * size() + i]) > TOLERANCE){

//                         std::cout.precision(14);
//                         std::cout << i << " " << j << " " << k << " " << c[k*(size() * size())+ j * size() + i] << std::endl;

                            //mat_[k][i][j] = round(c[k*(size() * size())+ j * size() + i], 1e-12);

                            mat_[k][i][j] = c[k*(size() * size())+ j * size() + i];
                        } else {
                            mat_[k][i][j] = 0.0;
                        }
                    }
                }
            }
        } else if (c.size() == elementList_.size()){
            for (Index k = 0; k < mat_.size(); k ++){ // z
                mat_[k] *= 0.0;
            }
            for (Index i = 0; i < c.size(); i ++){
                PolynomialElement< ValueType > e = elementList_[i];
                mat_[e.k_][e.i_][e.j_] = c[i];
            }
        } else {
            throwLengthError(WHERE_AM_I + " c size out of range " +
                                str(c.size()) + " needed " + str(powInt(mat_.size(), 3)) +
                                                    " or " + str(elementList_.size()));
        }
        this->fillElementList();
        return *this;
    }

    /*! Return all coefficients. */
    RVector coeff() const {
        RVector c(powInt(mat_.size(), 3));

        for (Index k = 0; k < mat_.size(); k ++){ // z
            for (Index j = 0; j < mat_[k].rows(); j ++){ // y
                for (Index i = 0; i < mat_[k][j].size(); i ++){ // x
                    c[k*(size() * size())+ j * size() + i] = mat_[k][i][j];
                }
            }
        }
        return c;
    }

protected:

    void init_(const Vector < ValueType > & ax, 
               const Vector < ValueType > & ay, 
               const Vector < ValueType > & az){
        Index maxDim = max(max(ax.size(), ay.size()), az.size());

        for (Index k = 0; k < maxDim; k ++){
            mat_.push_back(RMatrix(maxDim, maxDim));
            mat_[k] *= 0.0;
        }

        for (Index i = 0; i < ax.size(); i ++) mat_[0][i][0] = ax[i];
        for (Index j = 0; j < ay.size(); j ++) mat_[0][0][j] = ay[j];
        for (Index k = 0; k < az.size(); k ++) mat_[k][0][0] = az[k];
        fillElementList();
    }

    std::vector < Matrix < ValueType > > mat_;

    std::vector <PolynomialElement < ValueType > > elementList_;
};

template < class ValueType > bool operator == (const PolynomialFunction < ValueType > & a, const PolynomialFunction < ValueType > & b){
    if (a.elements().size() == b.elements().size()){
        typename std::vector <PolynomialElement < ValueType > >::const_iterator ita = a.elements().begin();
        typename std::vector <PolynomialElement < ValueType > >::const_iterator itb = b.elements().begin();

        for (; ita != a.elements().end(); ita ++, itb ++){
            if (*ita != *itb) return false;
        }
    } else {
        return false;
    }

    return true;
}

template < class ValueType > PolynomialFunction < ValueType >
operator - (const PolynomialFunction < ValueType > & f){

    PolynomialFunction < ValueType > h(Vector<ValueType>(f.size(), 0.0));

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = -f[k][i][j];
            }
        }
    }
    h.fillElementList();
    return h;
}

template < class ValueType > PolynomialFunction < ValueType >
operator * (const ValueType & val, const PolynomialFunction < ValueType > & f){

    PolynomialFunction < ValueType > h(RVector(f.size()), 0.0);

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = f[k][i][j] * val;
            }
        }
    }
    h.fillElementList();
    return h;
}

template < class ValueType > PolynomialFunction < ValueType >
operator * (const PolynomialFunction < ValueType > & f, const ValueType & val){

    PolynomialFunction < ValueType > h(RVector(f.size()), 0.0);

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = f[k][i][j] * val;
            }
        }
    }
    h.fillElementList();
    return h;
}

template < class ValueType > PolynomialFunction < ValueType >
operator + (const ValueType & val, const PolynomialFunction < ValueType > & f){
    PolynomialFunction < ValueType > h(RVector(f.size(), 0.0));

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = f[k][i][j];
            }
        }
    }
    h[0][0][0] += val;
    h.fillElementList();
    return h;
}

template < class ValueType > PolynomialFunction < ValueType >
operator + (const PolynomialFunction < ValueType > & f, const ValueType & val){
    PolynomialFunction < ValueType > h(RVector(f.size(), 0.0));

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = f[k][i][j];
            }
        }
    }
    h[0][0][0] += val;
    h.fillElementList();
    return h;
}

template < class ValueType > PolynomialFunction < ValueType >
operator + (const PolynomialFunction < ValueType > & f, const PolynomialFunction < ValueType > & g){

    PolynomialFunction < ValueType > h(RVector(max(f.size(),g.size()), 0.0));

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = f[k][i][j];
            }
        }
    }
    for (Index k = 0; k < g.size(); k ++){ // z
        for (Index j = 0; j < g[k].rows(); j ++){ // y
            for (Index i = 0; i < g[k][j].size(); i ++){ // x
                h[k][i][j] += g[k][i][j];
            }
        }
    }
    h.fillElementList();
    return h;
}

template < class ValueType > PolynomialFunction < ValueType >
operator - (const PolynomialFunction < ValueType > & f, const PolynomialFunction < ValueType > & g){

    PolynomialFunction < ValueType > h(RVector(max(f.size(),g.size()), 0.0));

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                h[k][i][j] = f[k][i][j];
            }
        }
    }
    for (Index k = 0; k < g.size(); k ++){ // z
        for (Index j = 0; j < g[k].rows(); j ++){ // y
            for (Index i = 0; i < g[k][j].size(); i ++){ // x
                h[k][i][j] -= g[k][i][j];
            }
        }
    }
    h.fillElementList();
    return h;
}

/*! Create new polynomial function for f(x,y,z) * g(x,y,z).  pls refactor with expressions */
template < class ValueType > PolynomialFunction < ValueType >
operator * (const PolynomialFunction < ValueType > & f, const PolynomialFunction < ValueType > & g){


    //TODO & TEST das muss besser werden weil es nicht immer f.size() + g.size() aufspannt. bsp. x*y statt x*x
    PolynomialFunction < ValueType > h(RVector(f.size() + g.size(), 0.0));

    for (Index k = 0; k < f.size(); k ++){ // z
        for (Index j = 0; j < f[k].rows(); j ++){ // y
            for (Index i = 0; i < f[k][j].size(); i ++){ // x
                for (Index gk = 0; gk < g.size(); gk ++){ // z
                    for (Index gj = 0; gj < g[gk].rows(); gj ++){ // y
                        for (Index gi = 0; gi < g[gk][gj].size(); gi ++){ // x
                            h[k + gk][i + gi][j + gj] += f[k][i][j] * g[gk][gi][gj];
                        }
                    }
                }
            }
        }
    }
    h.fillElementList();
    return h;
}

template < class ValueType > std::ostream & operator << (std::ostream & os, const PolynomialFunction < ValueType > & p){
    os << "f(x,y,z) = ";

    if (p.size() == 0) os << "0";

    bool first = true;
    for (Index k = 0; k < p.size() ; k ++){ // z
        for (Index j = 0; j < p[k].rows(); j ++){ // y
            for (Index i = 0; i < p[k][j].size(); i ++){ // x
                if (i > 0) first = false;

                std::string si = "";

                ValueType v = p[k][i][j];
                if (v == 0.0) continue;
                if (v > 0.0) si = "+";

                std::string xterm, yterm, zterm, val;
                if (i == 0) xterm = ""; else if (i == 1) xterm = "x"; else xterm = "x^" + str(i);
                if (j == 0) yterm = ""; else if (j == 1) yterm = "y"; else yterm = "y^" + str(j);
                if (k == 0) zterm = ""; else if (k == 1) zterm = "z"; else zterm = "z^" + str(k);

                val = str(v);

                if (first){
                    if (v == -1.0) val = "-1";
                    else if (v == 1.0) val = "1";
                } else {
                    if (::fabs(v + 1.0) < TOLERANCE) val = "-";
                    else if (::fabs(v - 1.0) < TOLERANCE) val = "";
                }

                os << si << val << xterm;
                if (yterm != "") os << yterm;
                if (zterm != "") os << zterm;
            }
        }
    }

    return os;
}


} // namespace GIMLI

#endif // _GIMLI_POLYNOMIAL__H
