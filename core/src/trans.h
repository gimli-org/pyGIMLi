/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Thomas Günther thomas@resistivity.net                                    *
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

#ifndef _GIMLI_TRANS__H
#define _GIMLI_TRANS__H

#include "gimli.h"
#include "vector.h"
#include "numericbase.h"

#define TRANSTOL 1e-8

namespace GIMLI{

/*! Basis Transform vector (unity), Identity transformation,
that means it do nothing. Just for prototyping f(a). */
template< class Vec > class Trans {
public:
//
    /*! Default constructur. */
    Trans() { }

    /*! Default destructor.*/
    virtual ~Trans() { }

    /*! Default call operator. \ref trans() */
    virtual Vec operator()(const Vec & a) const { return trans(a); }

    /*! Default call operator. \ref trans() */
    virtual double operator()(double x) const { return trans(x); }

    /*! Return forward transformation. Alias for \ref trans .*/
    inline Vec fwd(const Vec & f) const { return this->trans(f); }

    /*! Return forward transformation. Alias for \ref trans .*/
    inline double fwd(double f) const { return this->trans(f); }

    /*! Return inverse transformation. Shortcut for \ref invTrans .*/
    inline Vec inv(const Vec & f) const { return this->invTrans(f); }

    /*! For compatibility if inv is called with a single value. */
    inline double inv(double f) const { return this->inv(Vec(1, f))[0]; }


    /*! For compatibility if trans is called with a single value. */
    virtual double trans(double x) const { return this->trans(Vec(1, x))[0]; }

    /*! Apply and return the transformation. Return \f$ f(x) \f$*/
    virtual Vec trans(const Vec & x) const { return x; }

    /*! Return inverse transformation \f$ x(f) \f$*/
    virtual Vec invTrans(const Vec & f) const { return f; }

    /*! For compatibility if invTrans is called with a single value. */
    virtual double invTrans(double f) const { return this->invTrans(Vec(1, f))[0]; }

    /*! Return derivative \f$ \frac{\partial f}{\partial x}(x) \f$ */
    virtual Vec deriv(const Vec & x) const { return Vec(x.size(), 1.0); }

    /*! Update parameter by df: invTrans(f(a) + df). \n
    intrinsic function that have never to be overloaded. */
    Vec update(const Vec & a, const Vec & b) const {
        return (invTrans(trans(a) + b));
    }

    /*! Error of f(a) calculated by a and relative error \f$ da/a: df = | da * df/da | \f$ \n
    intrinsic function that could be overloaded */
    Vec error(const Vec & a, const Vec & daBya) const {
        if (daBya == Vec(a.size(), 0.0)) return Vec(a.size(), 1.0);
        return abs(Vec(a * daBya * deriv(a)));
    }

    /*! Alternative Version - brute force: df = | f(a + da) - f(a) | \n
    intrinsic functions that have never to be overloaded */
    Vec error_brute(const Vec & a, const Vec & daBya) const {
        return abs(trans(a * (daBya + 1.0)) - trans(a));
    }
};

/*! Base class for non-invertible transformations, e.g. transMult and transPlus \n
    making inverse transform using Newton's method */
template< class Vec > class TransNewton : public Trans < Vec >{
public:
    TransNewton(const int maxiter=10)
    : maxiter_(maxiter) { }

    virtual ~TransNewton(){ }

    virtual Vec trans(const Vec & a) const { }

    virtual Vec deriv(const Vec & a) const { }

    virtual Vec invTrans(const Vec & a) const {
        Vec p(Vec(a.size(), 0.01)); // guess vector
        Vec dm(trans(p) - a); // function to be nulled
        //std::cout << "Newton iter norm:";
        for(size_t i = 0; i < 10; i++) {
            p = p - dm / deriv(p); // Newton method
            dm = trans(p) - a;
            //std::cout << " " << norml2(dm);
        }
        //std::cout << std::endl;
        return p;
    }
protected:
    int maxiter_;
};

//** Fundamental transformation operations: linear, inverse, log, power, exp:

/*! Linear multiplication and addition
   \f$ f(x) = b + a \cdot x \f$. With factor \f$ a \f$ and offset \f$ b \f$. */
template< class Vec > class TransLinear : public Trans < Vec > {
public:
    TransLinear(const Vec & factor, double offset=0.0)
        : factor_(factor), offset_(Vec(factor.size(), offset)) {}

    virtual ~TransLinear() { }

    /*! Return \f$ f(x) = b + a \cdot x \f$.  */
    virtual Vec trans(const Vec & x) const {
//        std::cout << a.size() << " " << factor_.size() << " "
//        << offset_.size()<< std::endl;
        return x * factor_ + offset_;
    }

    /*! Return \f$ x(f) = (f - b)/a \f$ */
    virtual Vec invTrans(const Vec & f) const { return (f - offset_) / factor_;}

    /*! Return \f$ a = \frac{\partial}{\partial x} (b + a * x) \f$ */
    virtual Vec deriv(const Vec & x) const { return factor_; }

protected:
    Vec factor_;
    Vec offset_;
};

/*! Linear multiplication and addition
   f(a) = a * factor + offset. implementation with double */
template< class Vec > class TransLin : public Trans < Vec > {
public:
    TransLin(double factor=1.0, double offset=0.0)
        : factor_(factor), offset_(offset) {}

    virtual ~TransLin() { }

    virtual Vec trans(const Vec & a) const { return a * factor_ + offset_; }

    virtual Vec invTrans(const Vec & a) const { return (a - offset_) / factor_;}

    virtual Vec deriv(const Vec & a) const { return Vec(a.size(), factor_); }

protected:
    double factor_;
    double offset_;
};

/*! Power transformation f(a) = pow(a / a_0, n). */
template< class Vec > class TransPower : public Trans < Vec > {
public:
    TransPower(double npower=-1.0, double a0=1.0)
        : npower_(npower), a0_(a0) {}
    virtual ~TransPower() { }

    /*! Return \f$ f(x) = (x/a)^b \f$.  */
    virtual Vec trans(const Vec & x) const {
        return pow(x / a0_, npower_);}

    /*! Return \f$ x(f) = f^(1/b) \cdot a \f$. */
    virtual Vec invTrans(const Vec & f) const {
        return pow(f, 1.0 / npower_) * a0_;}

    /*! Return \f$ \f$. */
    virtual Vec deriv(const Vec & a) const {
        return pow(a / a0_, npower_ - 1.0) * npower_ / a0_;}
protected:
    double npower_;
    double a0_;
};

/*! Exponential transformation f(a) = f_0 exp(- a / a_0). */
template< class Vec > class TransExp : public Trans < Vec > {
public:
    /*!*/
    TransExp(double a0=1.0, double f0=1.0)
        : a0_(a0), f0_(f0) {}

    /*!*/
    virtual ~TransExp() {}

    /*!*/
    virtual Vec trans(const Vec & a) const {
        return exp(a / a0_ * (-1.0)) * f0_; }

    /*!*/
    virtual Vec invTrans(const Vec & f) const {
        return log(f / f0_) * a0_ * (-1.0); }

    /*! */
    virtual Vec deriv(const Vec & a) const {
        return trans(a) / a0_ * (-1.0); }
protected:
    double a0_;
    double f0_;
};

/*! Inverse of the parameter (e.g. velocity/slowness) f(a) = 1 / a. */
template< class Vec > class TransInv : public Trans < Vec > {
public:
    TransInv() { }
    virtual ~TransInv() { }

    virtual Vec trans(const Vec & a) const { return 1.0 / a; }
    virtual Vec invTrans(const Vec & f) const { return 1.0 / f; }
    virtual Vec deriv(const Vec & a) const { return -1.0 / a / a; }

};

/*! Logarithm of the parameter with the natural bound 0 or a self-defined. */
template< class Vec > class TransLog : public Trans < Vec > {
public:
    TransLog(double lowerbound=0.0) : lowerbound_(lowerbound) { }
    virtual ~TransLog() { }

    virtual Vec trans(const Vec & a) const {
        double lb1 = lowerbound_ * (1.0 + TRANSTOL);
        if (min(a) < lb1) {
            std::cerr << WHERE_AM_I << " Warning! " << min(a)
                      << " <=" << lowerbound_ << " lowerbound" << std::endl;
            Vec tmp(a);
            for (uint i = 0; i < a.size(); i ++){
                tmp[i] = max(a[i], lb1 );
            }
            return log(tmp - lowerbound_);
        }
        return log(a - lowerbound_);
    }

    virtual Vec invTrans(const Vec & f) const { return exp(f) + lowerbound_; }

    virtual Vec deriv(const Vec & a) const {
        double lb1 = lowerbound_ * (1.0 + TRANSTOL);
        if (min(a) < lb1) {
            std::cerr << WHERE_AM_I << " Warning! " << min(a)
                      << " <=" << lowerbound_ << " lowerbound" << std::endl;
            Vec tmp(a);
            for (uint i = 0; i < a.size(); i ++){
                tmp[i] = max(a[i], lb1 );
            }
            return 1.0 / (tmp - lowerbound_);
        }
        return 1.0 / (a - lowerbound_);
    }

    //** suggested by Friedel(2003), but deactivated since inherited by transLogLU
//    virtual Vec error(const Vec & a, const Vec & daBya) const { return log(1.0 + daBya); }

    inline void setLowerBound(double lb) { lowerbound_ = lb; }

    inline double lowerBound() const { return lowerbound_; }

protected:
    double lowerbound_;
};


/*! Range constraint transform functions: */
/*! Logarithmic barrier with upper and lower bound f(a) = log(a - a_u) - log(a_u - a). */
template< class Vec > class TransLogLU : public TransLog < Vec > {
public:
    TransLogLU(double lowerbound=0.0, double upperbound=0.0)
        : TransLog< Vec >(lowerbound), upperbound_(upperbound) { }
    virtual ~TransLogLU() { }

    /*!WHAT IS THIS?*/
    Vec rangify(const Vec & a) const;

    virtual Vec trans(const Vec & a) const;

    virtual Vec invTrans(const Vec & a) const;

    virtual Vec deriv(const Vec & a) const;

    inline void setUpperBound(double ub) { upperbound_ = ub; }

    inline double upperBound() const { return upperbound_; }

protected:
  double upperbound_;
};

/*! Implement specialized type traits in trans.cpp. */
template <> DLLEXPORT RVector TransLogLU<RVector>::rangify(const RVector & a) const;
template <> DLLEXPORT RVector TransLogLU<RVector>::trans(const RVector & a) const;
template <> DLLEXPORT RVector TransLogLU<RVector>::invTrans(const RVector & a) const;
template <> DLLEXPORT RVector TransLogLU<RVector>::deriv(const RVector & a) const;


/*! Cotangens barrier method, e.g. for water content (NMR) */
template< class Vec > class TransCotLU : public Trans < Vec > {
public:

    TransCotLU(double lowerbound=0.0, double upperbound=0.0)
        : lowerbound_(lowerbound), upperbound_(upperbound) { }

    virtual ~TransCotLU() { }

    virtual Vec trans(const Vec & a) const {
        Vec tmp(a);
        double fak = 1.00001;
        if (min(a) <= lowerbound_ ){
            std::cerr << WHERE_AM_I << " Warning! " << min(a) << " < " << lowerbound_ << " = lowerbound" << std::endl;
            for (uint i = 0; i < a.size(); i ++){
                tmp[i] = max(a[i], lowerbound_ * fak);
            }
        }
        if (max(a) >= upperbound_ ){
            std::cerr << WHERE_AM_I << " Warning! " << max(a) << " > " << upperbound_ << " = upperbound" << std::endl;
            for (uint i = 0; i < a.size(); i ++){
                tmp[i] = min(a[i], upperbound_ / fak);
            }
        }
        return cot((tmp - lowerbound_) / (upperbound_ - lowerbound_) * PI) * (-1.0);
    }

    virtual Vec invTrans(const Vec & a) const {
//    return acot(a * (-1.0)) * (upperbound_ - lowerbound_) / PI + lowerbound_;
        return atan(a) * (upperbound_ - lowerbound_) / PI + (lowerbound_ + upperbound_) / 2;
    }

    virtual Vec deriv(const Vec & a) const {
        return ((trans(a) * trans(a) + 1.0) * PI / (upperbound_ - lowerbound_));
    }

protected:
    double lowerbound_;
    double upperbound_;
};

/*! Different combined transformations: */
/*! Two nested TF's, e.g. log(G * resistance) or LogLU(CRIM(slowness)) */
template< class Vec > class TransNest : public Trans < Vec >{
public:
    TransNest(Trans< Vec > & T1, Trans< Vec > & T2) : T1_(&T1), T2_(&T2){
    }
    virtual ~TransNest(){ }

    virtual Vec trans(const Vec & a) const {
        return T1_->trans(T2_->trans(a));
    }

    virtual Vec invTrans(const Vec & a) const {
        return T2_->invTrans(T1_->invTrans(a));
    }

    virtual Vec deriv(const Vec & a) const {
        return T1_->deriv(T2_->trans(a)) * T2_->deriv(a);
    }

protected:
    Trans< Vec > * T1_;
    Trans< Vec > * T2_;
};

/*! Two added TF's f(a) = f1(a) + f2(a) */
template< class Vec > class TransAdd : public TransNewton < Vec >{
public:
    TransAdd(Trans< Vec > & T1, Trans< Vec > & T2) : T1_(&T1), T2_(&T2){ }
    virtual ~TransAdd(){ }

    virtual Vec trans(const Vec & a) const {
        return T1_->trans(a) + T2_->trans(a);
    }

    virtual Vec deriv(const Vec & a) const {
        return T1_->deriv(a) + T2_->deriv(a);
    }

protected:
    Trans< Vec > * T1_;
    Trans< Vec > * T2_;
};

/*! Two multiplied TF's f(a) = f1(a) * f2(a) (not necessarily stable). */
template< class Vec > class TransMult : public TransNewton < Vec >{
public:
    TransMult(Trans< Vec > & T1, Trans< Vec > & T2) : T1_(&T1), T2_(&T2){ }
    virtual ~TransMult(){ }

    virtual Vec trans(const Vec & a) const {
        return T1_->trans(a) * T2_->trans(a);
    }

    virtual Vec deriv(const Vec & a) const {
        return T1_->deriv(a) * T2_->trans(a) + T2_->deriv(a) * T1_->trans(a);
    }

protected:
    Trans< Vec > * T1_;
    Trans< Vec > * T2_;
};

/*! Petrophysical transform functions (specializations&combinations of other) */
/*! Complex refractive index model (CRIM) based on transLinear */
template< class Vec > class TransCRIM : public TransLinear < Vec > {
public:
    /*! petrophysical constructor using porosity and permittivity of matrix and water as input */
    TransCRIM(const Vec & phi, double ematrix = 5.0, double ewater = 81.0)
        : TransLinear< Vec >(1.0 / phi / 2.99e8 / (std::sqrt(ewater) - 1.0),
          (phi + (phi - 1.0) * std::sqrt(ematrix)) / (phi * (std::sqrt(ewater) - 1.0))){
    }
    virtual ~TransCRIM() { }
};

/*! Archie's law transforming resistivity into porosity. */
template< class Vec > class TransArchie : public TransPower < Vec > {
public:
    /*! petrophysical constructor using fluid resistivity as input */
    TransArchie(double resfluid, double mexp = 2.0)
        :TransPower< Vec >(resfluid, 1 / mexp){
    }

    virtual ~TransArchie() { }
};

/*! obsolete transform functions, e.g. for test reasons */
/*! transQuadrad, test class for checking transNewton */
template< class Vec > class TransQuadrat : public TransNewton < Vec >{
public:
    TransQuadrat(const int maxiter = 10) : TransNewton< Vec >(maxiter){ }
    virtual ~TransQuadrat(){ }
    virtual Vec trans(const Vec & a) const { return a * a; }
    virtual Vec deriv(const Vec & a) const { return a * 2.0; }
};

/*! Logarithm of parameter times factor (e.g. log(G * resistance)), better by transNest(transLog,transLinea) */
template< class Vec > class TransLogMult : public TransLog < Vec > {
public:
    TransLogMult(const Vec & factor, double lowerbound = 0.0)
        : TransLog< Vec >(lowerbound), factor_(factor) { }
    virtual ~TransLogMult() { }

    virtual Vec trans(const Vec & a) const { return TransLog< Vec >::trans(a * factor_); }
    virtual Vec invTrans(const Vec & a) const { return TransLog< Vec >::invTrans(a) / factor_; }
    virtual Vec deriv(const Vec & a) const { return TransLog< Vec >::deriv(a * factor_) * factor_; }

protected:
    Vec factor_;
};

/*! Tangens barrier method, replaced by TransCotLU, same results */
template< class Vec > class TransTanLU : public Trans < Vec > {
public:
    TransTanLU(double lowerbound=0.0, double upperbound=0.0)
        : lowerbound_(lowerbound), upperbound_(upperbound) {
    }
    virtual ~TransTanLU() { }

    virtual Vec trans(const Vec & a) const {
        Vec tmp(a);
        double tol = 1e-5;
        if (a.min() < lowerbound_ ){
            std::cerr << WHERE_AM_I << " Warning! " << a.min() << " < " << lowerbound_ << " lowerbound" << std::endl;
            for (int i = 0; i < a.size(); i ++){
                tmp[i] = max(a[i], lowerbound_) + tol;
            }
        }
        if (a.max() > upperbound_ ){
            std::cerr << WHERE_AM_I << " Warning! " << a.max() << " > " << upperbound_ << " lowerbound" << std::endl;
            Vec tmp(a);
            for (int i = 0; i < a.size(); i ++){
                tmp[i] = min(a[i], upperbound_) - tol;
            }
        }
        return tan(((a - lowerbound_) / (upperbound_ - lowerbound_) - 0.5) * PI);
    }

    virtual Vec invTrans(const Vec & a) const {
        return atan(a) / PI * (upperbound_ - lowerbound_) + (lowerbound_ + upperbound_) / 2;
    }

    virtual Vec deriv(const Vec & a) const {
        return (pow(trans(a), 2.0) + 1.0);
    }

protected:
    double lowerbound_;
    double upperbound_;
};

/*! Logarithm barrier of parameter times factor, better by transNest(transLogLU,transLinear) */
template< class Vec > class TransLogLUMult : public TransLogLU < Vec > {
public:
    TransLogLUMult(const Vec & factor,
                   double lowerbound=0.0,
                   double upperbound=0.0)
        : TransLogLU< Vec >(lowerbound, upperbound), factor_(factor){ }

    virtual ~TransLogLUMult() { }

    virtual Vec trans(const Vec & a) const {
        return TransLogLU< Vec >::trans(a * factor_);
    }

    virtual Vec invTrans(const Vec & a) const {
        return TransLogLU< Vec >::invTrans(a) / factor_ ;
    }

    virtual Vec deriv(const Vec & a) const {
        return TransLogLU< Vec >::deriv(a * factor_) * factor_;
    }

protected:
    Vec factor_;
};

/*! Auxillary cumulative transformation function using a vector of transformations*/
//! Very Slow. Refactor it!!
template< class Vec > class TransCumulative : public Trans < Vec >{
public:
    TransCumulative() { }

    virtual ~TransCumulative() { }

    virtual Vec trans(const Vec & a) const;

    virtual Vec invTrans(const Vec & a) const;

    virtual Vec deriv(const Vec & a) const;

    Index size() const { return transVec_.size(); }

    void clear() {
        transVec_.clear();
        slice_.clear();
        indices_.clear();
    }

    void add(Trans< Vec > & trans, Index size);

    void add(Trans< Vec > & trans, Index start, Index end);

    void add(Trans< Vec > & trans, const IVector & indices);

    /*!Return the containing transformationb object.*/
    Trans < Vec > & at(Index i) { return *transVec_.at(i);}

    const std::pair< Index, Index> & slice(Index i) const { return slice_.at(i);}

    const IVector & indices(Index i) const { return indices_.at(i);}

protected:
    std::vector < Trans< Vec > * > transVec_;
    std::vector < std::pair< Index, Index> > slice_;
    std::vector < IVector > indices_;
};

template <> DLLEXPORT RVector TransCumulative < RVector >::trans(const RVector & a) const;
template <> DLLEXPORT RVector TransCumulative < RVector >::invTrans(const RVector & a) const;
template <> DLLEXPORT RVector TransCumulative < RVector >::deriv(const RVector & a) const;

template <> DLLEXPORT void TransCumulative < RVector >::add(Trans< RVector > & trans, Index size);
template <> DLLEXPORT void TransCumulative < RVector >::add(Trans< RVector > & trans, Index start, Index end);
template <> DLLEXPORT void TransCumulative < RVector >::add(Trans< RVector > & trans, const IVector & indices);

typedef Trans < RVector > RTrans;
typedef TransLinear < RVector > RTransLinear;
typedef TransLin < RVector > RTransLin;
typedef TransPower < RVector > RTransPower;
typedef TransLog < RVector > RTransLog;
typedef TransLogLU < RVector > RTransLogLU;
typedef TransCotLU < RVector > RTransCotLU;
typedef TransCumulative< RVector > RTransCumulative;


} // namespace GIMLI


#endif // _GIMLI_TRANS__H

