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

#ifndef GIMLI_EXPRESSIONS__H
#define GIMLI_EXPRESSIONS__H

//** Idea taken from
// @Article{Veldhuizen95b,
//   author =       "Todd L. Veldhuizen",
//   title =        "Expression templates",
//   journal =      "C++ Report",
//   volume =       "7",
//   number =       "5",
//   pages =        "26--31",
//   month =        jun,
//   year =         "1995",
//   note =         "Reprinted in C++ Gems, ed. Stanley Lippman"
// }

#include <iostream>
#include <cmath>

namespace std{
    // make std::less happy .. this wont work without this namespace
    inline bool operator < (const std::complex< double > & a,
                            const std::complex< double > & b){
        return a.real() < b.real() || (!(b.real() < a.real()) && a.imag() < b.imag());
    }

    inline bool operator > (const std::complex< double > & a,
                            const std::complex< double > & b){
        return !(a < b);
    }
    inline bool operator >= (const std::complex< double > & a,
                             const std::complex< double > & b){
        return (a > b) || (a==b);
    }
    inline bool operator <= (const std::complex< double > & a,
                             const std::complex< double > & b){
        return (a < b) || (a==b);
    }
}

namespace GIMLI{

class ExprIdentity;
template< class A > class Expr;
template< class A > class ExprLiteral;

struct XAxis__ {};
struct YAxis__ {};
struct ZAxis__ {};

template <class T> struct Variable {
};

template<> struct Variable< XAxis__ > {
    double operator()(const double x) const { return x; }
    double operator()(const double x, const double y) const { return x; }
    double operator()(const double x, const double y, const double z) const { return x; }
};

template<> struct Variable< YAxis__ > {
    double operator()(const double x, const double y) const { return y; }
    double operator()(const double x, const double y, const double z) const { return y; }
};

template<> struct Variable< ZAxis__ > {
    double operator()(const double x, const double y, const double z) const { return z; }
};

typedef Expr< ExprIdentity > DPlaceholder;
typedef Expr< ExprIdentity > IntPlaceholder;
typedef Expr< ExprIdentity > Placeholder;
typedef Expr< Variable<XAxis__> > VariableX;
typedef Expr< Variable<YAxis__> > VariableY;
typedef Expr< Variable<ZAxis__> > VariableZ;

//**! Simple show the results of the expression
template< class Ex > void show(Expr< Ex > expr, double start, double end, double step = 1.0){
    for (double i = start; i < end; i += step) std::cout << expr(i) << std::endl;
}

template< class Ex > double sum(Expr< Ex > expr, double start, double end, double step){
    double result = 0.0;
    for (double i = start; i < end; i += step) result += expr(i);
    return result;
}

template< class Ex > double sum(Expr< Ex > expr, double start, double end){
    double result = 0.0;
    for (double i = start; i < end; i += 1.0) result += expr(i);
    return result;
}
#ifndef PI
//** I redefine PI here, because vector.h should independ from alternative non std::headers, Ca
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define PI2 6.2831853071795864769252867665590057683942387987502116419498891846
#define PI_2 1.5707963267948967384626433682795028841971193993751058209749445923
#endif

#define DEFINE_BINARY_OPERATOR__(OP, FUNCT) \
struct OP { template < class T1, class T2 > inline T1 operator()(const T1 & a, const T2 & b) const { return a FUNCT b; } };\

DEFINE_BINARY_OPERATOR__(PLUS, +)
DEFINE_BINARY_OPERATOR__(MINUS, -)
DEFINE_BINARY_OPERATOR__(MULT, *)
DEFINE_BINARY_OPERATOR__(DIVID, /)

#undef DEFINE_BINARY_OPERATOR__

// inline bool isEqual(const Complex & a, const Complex & b) { return (a == b); }
//inline bool isNonEqual(const double & a, const double & b) { return !isEqual(a, b); }
//inline bool isNonEqual(const Complex & a, const Complex & b) { return (a != b); }
template < class T > bool isEqual(const T & a, const T & b) { return a==b; }
template < > inline bool isEqual(const double & a, const double & b) {
    return (::fabs(a - b) < TOLERANCE); }
template < > inline bool isEqual(const Complex & a, const Complex & b) {
    return (::fabs(a.real() - b.real()) < TOLERANCE && 
            ::fabs(a.imag() - b.imag()) < TOLERANCE); }
template < class T > bool isNonEqual(const T & a, const T & b) { return !isEqual(a, b); }

template < class T > bool isLesser(const T & a, const T & b) { return a < b; }
template < class T > bool isGreater(const T & a, const T & b) { return a > b; }
template < class T > bool isLesserEqual(const T & a, const T & b) { return a < b || isEqual(a, b); }
template < class T > bool isGreaterEqual(const T & a, const T & b) { return a > b || isEqual(a, b); }


#ifdef _MSC_VER
template < class T > bool isInfNaN(const T & a){ return (isinf(a) || isnan(a)); }
template < class T > bool isInf(const T & a){ return (isinf(a));}
template < class T > bool isNaN(const T & a){ return (isnan(a));}
inline bool isNaN(const Complex & a){ return (isnan(a.real()) || isnan(a.imag())); }
inline bool isInf(const Complex & a){ return (isinf(a.real()) || isinf(a.imag())); }
#else
template < class T > bool isInfNaN(const T & a){ return (std::isinf(a) || std::isnan(a)); }
template < class T > bool isInf(const T & a){ return (std::isinf(a));}
template < class T > bool isNaN(const T & a){ return (std::isnan(a));}
inline bool isNaN(const Complex & a){ return (std::isnan(a.real()) || std::isnan(a.imag())); }
inline bool isInf(const Complex & a){ return (std::isinf(a.real()) || std::isinf(a.imag())); }
#endif

inline bool isInfNaN(const Complex & a){ return (isInfNaN(a.real()) || isInfNaN(a.imag())); }

inline double abs(const double a) { return std::fabs(a); }
inline double abs(const Complex & a) { return std::abs(a); }

inline double conj(const double & a) { return a; }
inline Complex conj(const Complex & a) { return std::conj(a); }

inline Complex RINT(const Complex & a) { THROW_TO_IMPL; return Complex(0); }
inline double RINT(const double & a) { return rint(a); }

template < class T > T roundTo(const T & a, const T & tol){ return RINT(a / tol) * tol; }

template < class T > T square(const T & a){ return a * a;}

inline double cot(const double & a) { return 1.0 / std::tan(a); }
inline double acot(const double & a) { return PI / 2.0 * std::atan(a); }
inline double sign(const double & a) { return a > 0.0 ? 1.0 : (a < 0.0 ? -1.0 : 0.0); }
inline double exp10(const double & a) { return std::pow(10.0, a); }


// template < class T > inline T square(const T & a) { return a * a; }
// template < class T > inline T cot(const T & a) { return 1.0 / std::tan(a); }
// template < class T > inline T acot(const T & a) { return PI / 2.0 * std::atan(a); }
#define DEFINE_UNARY_OPERATOR__(OP, FUNCT) \
struct OP { template < class T > T operator()(const T & a) const { return FUNCT(a); } };\

DEFINE_UNARY_OPERATOR__(ACOT, acot)
DEFINE_UNARY_OPERATOR__(ABS_, std::fabs)
DEFINE_UNARY_OPERATOR__(SIN , std::sin)
DEFINE_UNARY_OPERATOR__(COS , std::cos)
DEFINE_UNARY_OPERATOR__(COT , cot)
DEFINE_UNARY_OPERATOR__(TAN , std::tan)
DEFINE_UNARY_OPERATOR__(ATAN, std::atan)
DEFINE_UNARY_OPERATOR__(TANH, std::tanh)
DEFINE_UNARY_OPERATOR__(LOG , std::log)
DEFINE_UNARY_OPERATOR__(LOG10, std::log10)
DEFINE_UNARY_OPERATOR__(EXP , std::exp)
DEFINE_UNARY_OPERATOR__(EXP10, exp10)
DEFINE_UNARY_OPERATOR__(SQRT, std::sqrt)
DEFINE_UNARY_OPERATOR__(SIGN, sign)
DEFINE_UNARY_OPERATOR__(SQR,  square)

#undef DEFINE_UNARY_OPERATOR__

#define DEFINE_UNARY_IF_FUNCTION__(OP, FUNCT) \
struct OP { template < class T > bool operator()(const T & a) const { return FUNCT(a); } }; \

DEFINE_UNARY_IF_FUNCTION__(ISINF, isInf)
DEFINE_UNARY_IF_FUNCTION__(ISNAN, isNaN)
DEFINE_UNARY_IF_FUNCTION__(ISINFNAN, isInfNaN)

#undef DEFINE_UNARY_IF_FUNCTION__

/*! Expr is a wrapper class which contains a more interesting expression type,
such as ExprIdentity, ExprLiteral, unary (UnaryExprOp) or  binaray expression operator (BinaryExprOp). */
template< class A > class Expr {
public:
    Expr() : a_(A()){ }

    Expr(const A & a) : a_(a){ }

    template < class ValueType > inline ValueType operator()(const ValueType & x) const { return a_(x); }
    template < class ValueType > inline ValueType operator()(const ValueType & x, const ValueType & y) const { return a_(x, y); }
    template < class ValueType > inline ValueType operator()(const ValueType & x, const ValueType & y, const ValueType & z) const { return a_(x, y, z); }
private:
    A a_;
};

class ExprIdentity {
public:
    ExprIdentity(){ }

    template < class ValueType > inline ValueType operator()(const ValueType & x) const { return x; }
};

template < class ValueType > class ExprLiteral {
public:
    ExprLiteral(const ValueType & value): value_(value){ }

    inline ValueType operator()(const ValueType & x) const { return value_; }
    inline ValueType operator()(const ValueType & x, const ValueType & y) const { return value_; }
    inline ValueType operator()(const ValueType & x, const ValueType & y, const ValueType & z) const { return value_; }
private:
    ValueType value_;
};

template< class A, class Op > class UnaryExprOp {
public:
    UnaryExprOp(const A & a) : a_(a){ }
    template < class ValueType > ValueType operator() (const ValueType & x) const {
        return Op()(a_(x));
    }
    template < class ValueType > ValueType operator() (const ValueType & x, const ValueType & y) const {
        return Op()(a_(x, y));
    }
    template < class ValueType > ValueType operator() (const ValueType & x, const ValueType & y, const ValueType & z) const {
        return Op()(a_(x, y, z));
    }
private:
    A a_;
};

template< class A, class B, class Op > class BinaryExprOp {
public:
    BinaryExprOp(const A & a, const B & b) : a_(a), b_(b){ }

    template < class ValueType > ValueType operator() (const ValueType & x) const {
        return Op()(a_(x), b_(x));
    }
    template < class ValueType > ValueType operator() (const ValueType & x, const ValueType & y) const {
        return Op()(a_(x, y), b_(x, y));
    }
    template < class ValueType > ValueType operator() (const ValueType & x, const ValueType & y, const ValueType & z) const {
        return Op()(a_(x, y, z), b_(x, y, z));
    }
private:
    A a_;
    B b_;
};

#define DEFINE_UNARY_EXPR_OPERATOR__(OP, FUNCT) \
template< class A > Expr < UnaryExprOp < Expr< A >, FUNCT > >\
    OP(const Expr< A > & a) { \
    typedef UnaryExprOp< Expr< A >, FUNCT > ExprT; \
    return Expr< ExprT >(ExprT(a));\
}\

DEFINE_UNARY_EXPR_OPERATOR__(abs,   ABS_)
DEFINE_UNARY_EXPR_OPERATOR__(acot,  ACOT)
DEFINE_UNARY_EXPR_OPERATOR__(atan,  ATAN)
DEFINE_UNARY_EXPR_OPERATOR__(cos,   COS)
DEFINE_UNARY_EXPR_OPERATOR__(cot,   COT)
DEFINE_UNARY_EXPR_OPERATOR__(exp,   EXP)
DEFINE_UNARY_EXPR_OPERATOR__(exp10, EXP10)
DEFINE_UNARY_EXPR_OPERATOR__(fabs,  ABS_)
DEFINE_UNARY_EXPR_OPERATOR__(log,   LOG)
DEFINE_UNARY_EXPR_OPERATOR__(log10, LOG10)
DEFINE_UNARY_EXPR_OPERATOR__(sign, SIGN)
DEFINE_UNARY_EXPR_OPERATOR__(sin,   SIN)
DEFINE_UNARY_EXPR_OPERATOR__(sqrt,  SQRT)
DEFINE_UNARY_EXPR_OPERATOR__(square, SQR)
DEFINE_UNARY_EXPR_OPERATOR__(tan,   TAN)
DEFINE_UNARY_EXPR_OPERATOR__(tanh,  TANH)

#undef DEFINE_UNARY_EXPR_OPERATOR__


#define DEFINE_EXPR_OPERATOR__(OP, FUNCT) \
template< class ValueType, class A > Expr < BinaryExprOp < ExprLiteral< ValueType >, Expr< A >, FUNCT > >\
    operator OP (const ValueType & x, const Expr< A > & a) {\
    typedef BinaryExprOp< ExprLiteral< ValueType >, Expr< A >, FUNCT > ExprT;\
    return Expr< ExprT >(ExprT(ExprLiteral< ValueType >(x), a));\
}\
template< class ValueType, class A > Expr < BinaryExprOp < Expr< A >, ExprLiteral< ValueType >, FUNCT > >\
    operator OP (const Expr< A > & a, const ValueType & x) {\
    typedef BinaryExprOp< Expr< A >, ExprLiteral< ValueType >, FUNCT > ExprT;\
    return Expr< ExprT >(ExprT(a, ExprLiteral< ValueType >(x)));\
}\
template < class A, class B > Expr< BinaryExprOp< Expr< A >, Expr< B >, FUNCT > >\
    operator OP (const Expr< A > & a, const Expr< B > & b){\
    typedef BinaryExprOp< Expr< A >, Expr< B >, FUNCT > ExprT;\
    return Expr< ExprT >(ExprT(a, b));\
}\

DEFINE_EXPR_OPERATOR__(+, PLUS)
DEFINE_EXPR_OPERATOR__(-, MINUS)
DEFINE_EXPR_OPERATOR__(*, MULT)
DEFINE_EXPR_OPERATOR__(/, DIVID)

#undef DEFINE_EXPR_OPERATOR__

} // namespace GIMLI

#endif
