/******************************************************************************
 *   Copyright (C) 2009-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_CURVEFITTING__H
#define _GIMLI_CURVEFITTING__H

#include "gimli.h"

#include "polynomial.h"

namespace GIMLI{

/*! Abstract function class, i.e., functor */
template < class ArgType , class ValueType > class DLLEXPORT Function {
public:
    Function() { }

    virtual ~Function() { }

    Function(const Function & funct){ copy_(funct); }

    Function & operator = (const Function & funct){
        if (this != & funct){
            copy_(funct);
        } return * this;
    }

    inline ValueType operator() (const ArgType & arg) const { return this->getValue(arg); }

    inline Vector < ValueType > operator() (const Vector < ArgType > & arg) const { return this->getValue(arg); }

    virtual ValueType getValue(const ArgType & arg) const = 0;

    virtual Vector < ValueType > getValue(const Vector < ArgType > & arg) const = 0;

protected:
    virtual void copy_(const Function & funct) { };
};

class DLLEXPORT HarmonicFunction : public Function< double, double >{
public:
    HarmonicFunction(const RVector & coeff, double xmin, double xmax);

    virtual ~HarmonicFunction();

    virtual double getValue(const double & arg) const;

    virtual RVector getValue(const RVector & arg) const;

    void setCoefficients(const RVector & coeff);

    inline const RVector & coefficients() const { return coeff_; }

    inline void setXMin(double xmin){ xMin_ = xmin; }

    inline double xMin() const { return xMin_; }

    inline void setXMax(double xmax){ xMax_ = xmax; }

    inline double xMax() const { return xMax_; }

protected:
    virtual void copy_(const HarmonicFunction & funct);

    RVector coeff_;
    size_t nHarmonic_;
    double xMin_;
    double xMax_;
};

class DLLEXPORT HarmonicModelling : public ModellingBase {
public:
    /*! constructor, nh: number of coefficients, xvec: abscissa, */
    HarmonicModelling(size_t nh, const RVector & tvec, bool verbose = false);

    virtual ~HarmonicModelling(){ }

    /*! the main thing - the forward operator: return f(x) */
    virtual RVector response(const RVector & par);

    /*! an additional forward operator for another time basis */
    virtual RVector response(const RVector & par, const RVector tvec);

    /*! optional: generation of Jacobian matrix, uncomment for default behavior (brute force) */
    virtual void createJacobian(const RVector & model);

    /*! Define the start model */
    inline virtual RVector startModel(){ return RVector(np_, 0.0); }

protected:

    RVector t_; //! abscissa vector x
    RMatrix A_; //! function matrix
    double tMin_;
    double tMax_;
    size_t nh_, nt_, np_;
};

//! Multidimensional polynomial modelling class.
/*! Multidimensional polynomial modelling class. 1, 2 or 3 dimensional polynomial modelling operator based on \ref PolynomialFunction */
class DLLEXPORT PolynomialModelling : public ModellingBase {
public:

    PolynomialModelling(uint dim, uint nCoeffizient, const std::vector < RVector3 > & referencePoints,
                        const RVector & startModel);

    virtual RVector response(const RVector & par);

    /*! Create starting model. The dimension is recognized here \n
        * one-dimensional:         p[0][0][*] = 1 \n
        * two-dimensional:         p[0][*][*] = 1 \n
        * three-dimensional:       p[*][*][*] = 1 \n
    */
    virtual RVector startModel();

    /*! Return read only reference to the Polynomial Functor. */
    const PolynomialFunction< double > & polynomialFunction() const { return f_; }

    /*! Constrain to polynomial function based on Pascal's triangle, i.e., forbid x^iy^jz^k with (i +j +k) > size */
    void setPascalsStyle(bool is) { pascalTriangle_ = is; }

    /*! Constrain the polynomial function serendipity style Pascal's triangle, i.e., forbid x^iy^jz^k with (i +j +k) > size+1.
     *  Just work if setPascalsStyle is set.
     */
    void setSerendipityStyle(bool is) { serendipityStyle_ = is; }

    void setPowCombinationTmp(int i) { powCombination_ = i; }

protected:
    uint dim_;
    std::vector < RVector3 > referencePoints_;

    PolynomialFunction< double > f_;

    bool pascalTriangle_;
    bool serendipityStyle_;
    uint powCombination_;
};

} // namespace GIMLI{

#endif // _GIMLI_CURVEFITTING__H
