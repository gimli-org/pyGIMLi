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

#include "curvefitting.h"

#include "pos.h"

#include "regionManager.h"
#include "vectortemplates.h"

namespace GIMLI {

HarmonicFunction::HarmonicFunction(const RVector & coeff, double xmin, double xmax)
: xMin_(xmin), xMax_(xmax){
    setCoefficients(coeff);
}

HarmonicFunction::~HarmonicFunction(){
}

void HarmonicFunction::setCoefficients(const RVector & coeff){

    nHarmonic_ = coeff.size() / 2;
    if (((double)coeff.size() / 2.0 - nHarmonic_) > TOLERANCE){
        throwError(WHERE_AM_I + " coefficients size is uneven" + str(coeff.size()));
    }
    coeff_ = coeff;
}

double HarmonicFunction::getValue(const double & arg) const {
    double ret = coeff_[0];

    double tOne = (arg - xMin_) / (xMax_ - xMin_);

    ret += tOne * coeff_[1];

    for (size_t j = 1 ; j < nHarmonic_; j++){
        ret += ::cos(tOne * PI2 * j) * coeff_[j * 2];
        ret += ::sin(tOne * PI2 * j) * coeff_[j * 2 + 1];
    }

    return ret;
}

RVector HarmonicFunction::getValue(const RVector & arg) const {
    RVector ret(arg.size(), coeff_[0]);

    RVector tOne((arg - xMin_) / (xMax_ - xMin_));

    ret += tOne * coeff_[1];

    for (size_t j = 1 ; j < nHarmonic_; j++){
        ret += cos(tOne * PI2 * j) * coeff_[j * 2];
        ret += sin(tOne * PI2 * j) * coeff_[j * 2 + 1];
    }

    return ret;
}

void HarmonicFunction::copy_(const HarmonicFunction & funct){
    xMin_ = funct.xMin();
    xMax_ = funct.xMax();
    this->setCoefficients(funct.coefficients());

}


HarmonicModelling::HarmonicModelling(size_t nh, const RVector & tvec, bool verbose)
: ModellingBase(verbose),
    t_(tvec), tMin_(min(tvec)), tMax_(max(tvec)), nh_(nh), np_(2 * nh + 2) {

    regionManager().setParameterCount(np_);
    A_.clear();
    nt_ = tvec.size();

    //! constant vector of 1 -- offset
    RVector one(nt_, 1.0);
    A_.push_back(one); //! const

    //! vector linearly ascending from 0 (tmin) to 1 (tmax) -- drift
    double tMin = min(tvec), tMax = max(tvec);
    RVector tOne((t_ - tMin) / (tMax - tMin)); //** wieso nicht so:
//    RVector tOne((t_ - tMin_) / (tMax_ - tMin_));

    A_.push_back(tOne);

    //! harmonic functions cos/sin(n pi t)
    for (size_t j = 1 ; j <= nh_ ; j++){
        one = cos(tOne * PI2 * j);
        A_.push_back(one);
        one = sin(tOne * PI2 * j);
        A_.push_back(one);
    }
}

RVector HarmonicModelling::response(const RVector & par){
    return transMult (A_ , par);
}

RVector HarmonicModelling::response(const RVector & par, const RVector tvec){
    RVector ret(tvec.size(), par[0]);

    RVector tOne((tvec - tMin_) / (tMax_ - tMin_));

    ret += tOne * par[1];

    for (size_t j = 1 ; j <= nh_ ; j++){
        ret += cos(tOne * PI2 * j) * par[j * 2];
        ret += sin(tOne * PI2 * j) * par[j * 2 + 1];
    }
    return ret;
}

void HarmonicModelling::createJacobian(const RVector & model) {
    //!! jacobian = transpose(A);
    RMatrix * jacobian = dynamic_cast < RMatrix * > (jacobian_);

    if (jacobian->rows() != nt_ || jacobian->cols() != np_) {
        jacobian->resize(nt_, np_);

        for (size_t i = 0 ; i < np_ ; i++){
            for (size_t j = 0 ; j < nt_ ; j++){
                (*jacobian)[j][i] = A_[i][j];
            }
        }
    }
}

PolynomialModelling::PolynomialModelling(uint dim, uint nCoeffizient,
                                         const std::vector < RVector3 > & referencePoints,
                                         const RVector & startModel)
    : dim_(dim), referencePoints_(referencePoints) 
       //f_(PolynomialFunction < double >(nCoeffizient)) 
       {
    f_ = PolynomialFunction < double >(nCoeffizient);
    pascalTriangle_ = false;
    serendipityStyle_ = false;
    startModel_ = startModel;
    powCombination_ = 0;
    this->regionManager().setParameterCount(powInt(nCoeffizient, 3));
}

RVector PolynomialModelling::response(const RVector & par){
    return f_.fill(round(par, TOLERANCE))(referencePoints_);
}

RVector PolynomialModelling::startModel(){
    if (startModel_.size() == powInt(f_.size(), 3)) return startModel_;

    RVector p(powInt(f_.size(), 3), 0.0);
    f_.clear();

    p.setVal(1.0, 0, (SIndex)(powInt(f_.size(), dim_)));

    if (pascalTriangle_){
        for (Index k = 0; k < f_.size(); k ++){
            for (Index j = 0; j < f_.size(); j ++){
                for (Index i = 0; i < f_.size(); i ++){
                    //**  remove elements outside of the Pascal's triangle
                    if (powCombination_ > 0) {
                        if ((i + j + k) > powCombination_) p[k*(f_.size() * f_.size()) + j * f_.size() + i] = 0.0;
                    } else {
                        if ((i + j + k) >= (f_.size() + serendipityStyle_ * (dim_ -1))) p[k*(f_.size() * f_.size()) + j * f_.size() + i] = 0.0;
                    }
                }
            }
        }
    }

    return p;
}

} // namespace GIMLI{
