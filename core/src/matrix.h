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

#pragma once

// #if USE_EIGEN3
//     #include <Eigen/Dense>
// #endif


#include "cassert"
#include "gimli.h"
#include "pos.h"
#include "vector.h"

#include <cstddef>
#include <memory>

#include <cstring>
#include <fstream>
#include <iostream>
#include <cerrno>

#ifdef USE_THREADS
    #if USE_BOOST_THREAD
        #include <boost/thread.hpp>
    #endif // USE_BOOST_THREAD
#endif // USE_THREADS


namespace GIMLI{

/*
// #if USE_EIGEN3

// // Matrix<int, Dynamic, Dynamic, RowMajor> RowMatrixXi;
//     typedef Eigen::Matrix < double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > SmallMatrix;
//     //typedef Matrix < double > SmallMatrix;
//     //typedef Matrix < double > SmallMatrix;

//     void toEigenMatrix(const Matrix < double > & m, SmallMatrix & r);
//     void toRMatrix(const SmallMatrix & m, Matrix < double > & r);
//     void toRVector(const Eigen::VectorXd & m, RVector & r, double b=0.0);

// #define SET_MAT_ROW_SLICE(A, row, RVec, b, start, end) \
//     A(row, Eigen::seq(start, end-1)) = \
//                 Eigen::Map <const Eigen::VectorXd>(&RVec[0], RVec.size()) * b;
// #define ADD_MAT_ROW_SLICE(A, row, RVec, b, start, end) \
//     A(row, Eigen::seq(start, end-1)) += \
//                 Eigen::Map <const Eigen::VectorXd>(&RVec[0], RVec.size()) * b;
// #define MAT_TRANS_ADD(A, B) A += B.transpose();
// #define MAT_ROW_IMUL(A, r, b) A(r, Eigen::all) *= b;
// // #define RVEC_ASSIGN_MAT_ROW_MUL(c, A, k, b) toRVector(A.row(k) * b, c);
// // #define RVEC_IADD_MAT_ROW_MUL(c, A, k, b) toRVector(A.row(k) * b, c, 1.0);

// #else
    // typedef Matrix < double > SmallMatrix;

// #define SET_MAT_ROW_SLICE(A, row, RVec, b, start, end) \
//                            A(row).setVal(RVec * b, start, end);
// #define ADD_MAT_ROW_SLICE(A, row, RVec, b, start, end) \
//                            A(row).addVal(RVec * b, start, end);
// #define MAT_TRANS_ADD(A, B) A.transAdd(B);
// #define MAT_ROW_IMUL(A, r, b) A(r) *= b;
// #define RVEC_ASSIGN_MAT_ROW_MUL(c, A, k, b) c = A(k) * b;
// #define RVEC_IADD_MAT_ROW_MUL(c, A, k, b) c += A(k) * b;
*/
// #endif

template < class ValueType >
bool load(Matrix < ValueType > & A, const std::string & filename);

template < class ValueType > class DLLEXPORT Matrix3 {
public:
    ValueType mat_[9];

    Matrix3()
        : valid_(false){ clear(); }

    Matrix3(const Matrix3 < ValueType > & m ){
        mat_[0] = m.mat_[0]; mat_[1] = m.mat_[1]; mat_[2] = m.mat_[2];
        mat_[3] = m.mat_[3]; mat_[4] = m.mat_[4]; mat_[5] = m.mat_[5];
        mat_[6] = m.mat_[6]; mat_[7] = m.mat_[7]; mat_[8] = m.mat_[8];
    }

    inline ValueType & operator [](Index i){ return mat_[i];}
    inline const ValueType & operator [](Index i) const { return mat_[i];}

    inline const ValueType & operator ()(Index i, Index j) const {
        THROW_TO_IMPL
        return mat_[0];
        }

    #define DEFINE_UNARY_MOD_OPERATOR__(OP, NAME) \
    inline Matrix3 < ValueType > & operator OP##= (const ValueType & val) { \
        for (Index i = 0; i < 9; i += 3) {\
            mat_[i] OP##= val; mat_[i+1] OP##= val; mat_[i+2] OP##= val; } return *this; \
    }

    DEFINE_UNARY_MOD_OPERATOR__(+, PLUS)
    DEFINE_UNARY_MOD_OPERATOR__(-, MINUS)
    DEFINE_UNARY_MOD_OPERATOR__(/, DIVID)
    DEFINE_UNARY_MOD_OPERATOR__(*, MULT)

    #undef DEFINE_UNARY_MOD_OPERATOR__

    void clear(){
        mat_[0] = 0.0; mat_[1] = 0.0; mat_[2] = 0.0;
        mat_[3] = 0.0; mat_[4] = 0.0; mat_[5] = 0.0;
        mat_[6] = 0.0; mat_[7] = 0.0; mat_[8] = 0.0;
    }
    inline Index rows() const { return 3; }
    inline Index cols() const { return 3; }

    inline const Vector < ValueType > col(Index i) const {
        Vector < ValueType >ret(3);
        ret[0] = mat_[i]; ret[1] = mat_[3 + i]; ret[2] = mat_[6 + i];
        return ret;
    }

    inline const Vector < ValueType > row(Index i) const {
        Vector < ValueType >ret(3);
        ret[0] = mat_[i * 3]; ret[1] = mat_[i * 3 + 1]; ret[2] = mat_[i * 3 + 2];
        return ret;
    }

    // inline void setVal(const RVector & v, Index i){
    //     __M
    //     log(Warning, "deprecated");
    //     mat_[i * 3] = v[0]; mat_[i * 3 + 1] = v[1]; mat_[i * 3 + 2] = v[2];
    // }
    inline void setVal(Index i, const RVector & v){
        mat_[i * 3] = v[0]; mat_[i * 3 + 1] = v[1]; mat_[i * 3 + 2] = v[2];
    }

    inline void setValid(bool v){valid_ = v;}

    inline bool valid() const {return valid_;}

    inline ValueType det() const {
        return mat_[0] * (mat_[4] * mat_[8] - mat_[5] * mat_[7]) -
               mat_[1] * (mat_[3] * mat_[8] - mat_[5] * mat_[6]) +
               mat_[2] * (mat_[3] * mat_[7] - mat_[4] * mat_[6]);
    }

protected:

    bool valid_;
};


template < class ValueType >
std::ostream & operator << (std::ostream & str, const Matrix3 < ValueType > & vec){
    for (Index i = 0; i < 3; i ++)
        str << vec[i] << " ";
    str << std::endl;
    for (Index i = 0; i < 3; i ++)
        str << vec[3 + i] << " ";
    str << std::endl;
    for (Index i = 0; i < 3; i ++)
        str << vec[6 + i] << " ";
    str << std::endl;
    return str;
}

template < class ValueType > \
Pos operator * (const Matrix3 < ValueType > & A, const Pos & b) {
    return Pos (A[0] * b[0] + A[1] * b[1] + A[2] * b[2],
                             A[3] * b[0] + A[4] * b[1] + A[5] * b[2],
                             A[6] * b[0] + A[7] * b[1] + A[8] * b[2]);
}

#define DEFINE_MATRIX_MV_INTERFACE(MAT, VALTYPE) \
/*! c = alpha * A*b + beta * c */ \
void DLLEXPORT mult(const MAT & A,  \
                    const Vector<VALTYPE> & b, Vector<VALTYPE> & c,\
                    const VALTYPE & alpha=1.0, const VALTYPE & beta=0.0,\
                    Index bOff = 0, Index cOff = 0);\
/*! c = alpha * A.T*b + beta * c */ \
void DLLEXPORT transMult(const MAT & A,  \
                    const Vector<VALTYPE> & b, Vector<VALTYPE> & c,\
                    const VALTYPE & alpha=1.0, const VALTYPE & beta=0.0,\
                    Index bOff = 0, Index cOff = 0);\

DEFINE_MATRIX_MV_INTERFACE(RDenseMatrix, double)
DEFINE_MATRIX_MV_INTERFACE(CDenseMatrix, Complex)
DEFINE_MATRIX_MV_INTERFACE(RMatrix, double)
DEFINE_MATRIX_MV_INTERFACE(CMatrix, Complex)
DEFINE_MATRIX_MV_INTERFACE(RSparseMapMatrix, double)
DEFINE_MATRIX_MV_INTERFACE(CSparseMapMatrix, Complex)
DEFINE_MATRIX_MV_INTERFACE(RSparseMatrix, double)
DEFINE_MATRIX_MV_INTERFACE(CSparseMatrix, Complex)
DEFINE_MATRIX_MV_INTERFACE(RBlockMatrix, double)
DEFINE_MATRIX_MV_INTERFACE(CBlockMatrix, Complex)
#undef DEFINE_MATRIX_MV_INTERFACE

#define DEFINE_MATRIX_MM_INTERFACE(MAT, VALTYPE) \
/*! C = alpha * A*B + beta * C */ \
void DLLEXPORT mult(const MAT & A,  \
                    const MAT & b, MAT & c,\
                    const VALTYPE & alpha=1.0, const VALTYPE & beta=0.0);\
void DLLEXPORT transMult(const MAT & A,  \
                    const MAT & b, MAT & c,\
                    const VALTYPE & alpha=1.0, const VALTYPE & beta=0.0);\

DEFINE_MATRIX_MM_INTERFACE(RDenseMatrix, double)
DEFINE_MATRIX_MM_INTERFACE(CDenseMatrix, Complex)
DEFINE_MATRIX_MM_INTERFACE(RMatrix, double)
DEFINE_MATRIX_MM_INTERFACE(CMatrix, Complex)

#undef DEFINE_MATRIX_MM_INTERFACE

//! Interface class for matrices.
/*! Pure virtual interface class for matrices.
 * If you want your own Jacobian matrix to be used in \ref Inversion or \ref ModellingBase
 you have to derive your matrix from this class and implement all necessary members. */
class DLLEXPORT MatrixBase{
public:

    /*! Default constructor. */
    MatrixBase(bool verbose=false)
        : verbose_(verbose), _rows(0), _cols(0) {}

    /*! Default constructor. */
    MatrixBase(Index rows, Index cols, bool verbose=false)
        : verbose_(verbose), _rows(rows), _cols(cols) {}

    /*! Default destructor. */
    virtual ~MatrixBase(){}

    /*! Return entity rtti value. */
    virtual uint rtti() const { return GIMLI_MATRIXBASE_RTTI; }

    void setVerbose(bool verbose){ verbose_ = verbose; }

    bool verbose() const { return verbose_; }

    /*! Return number of rows */
    inline Index size() const { return this->rows(); }

    /*! Return number of cols */
    inline Index rows() const {
       return this->_rows;
    }

    /*! Return number of colums */
    inline Index cols() const {
        return this->_cols;
    }

    /*! Resize this matrix to rows, cols */
    virtual void resize(Index rows, Index cols){
       this->_rows = rows;
       this->_cols = cols;
    }

    /*! Fill Vector with 0.0. Don't change size.*/
    virtual void clean() {
        this->_rows = 0;
        this->_cols = 0;
        log(Warning, "no clean() implemented for: ", typeid(*this).name());
    }
    /*! Fill Vector with 0.0. Don't change size. For compatibility with eigen interface.*/
    inline void setZero() { this->clean(); }
    /*! Clear the data, set size to zero and frees memory. */
    virtual void clear() {
        this->_rows = 0;
        this->_cols = 0;
        log(Warning, "no clear() implemented for: ", typeid(*this).name());
    }

    /*! Return this * a. For being numpy api-compatible  */
    virtual RVector dot(const RVector & a) const {
        return this->mult(a);
    }

    /*! Multiplication c = alpha * (A*b) + beta * c. */
    virtual void mult(const RVector & b, RVector & c,
                     const double & alpha=1.0,
                     const double & beta=0.0,
                     Index bOff=0, Index cOff=0) const {
        log(Warning, "no mult(const & b, const & b, alpha, beta, bOff, cOff) implemented for: ", typeid(*this).name());
    }
    /*! Multiplication c = alpha * (A*b) + beta * c. */
    virtual void transMult(const RVector & b, RVector & c,
                     const double & alpha=1.0,
                     const double & beta=0.0,
                     Index bOff=0, Index cOff=0) const {
        log(Warning, "no transMult(const & b, const & b, alpha, beta, bOff, cOff) implemented for: ", typeid(*this).name());
    }
    /*! Multiplication c = alpha * (A*b) + beta * c. */
    virtual void mult(const CVector & b, CVector & c,
                     const Complex & alpha=1.0,
                     const Complex & beta=0.0,
                     Index bOff=0, Index cOff=0) const {
        log(Warning, "no mult(const & b, const & b, alpha, beta, bOff, cOff) implemented for: ", typeid(*this).name());
    }
    /*! Multiplication c = alpha * (A*b) + beta * c. */
    virtual void transMult(const CVector & b, CVector & c,
                     const Complex & alpha=1.0,
                     const Complex & beta=0.0,
                     Index bOff=0, Index cOff=0) const {
        log(Warning, "no transMult(const & b, const & b, alpha, beta, bOff, cOff) implemented for: ", typeid(*this).name());
    }

    /*! Return this * a  */
    virtual RVector mult(const RVector & b) const {
        RVector ret(this->rows(), 0.0);
        this->mult(b, ret);
        return ret;
    }
    /*! Return this * a  */
    virtual CVector mult(const CVector & b) const {
        CVector ret(this->rows(), 0.0);
        this->mult(b, ret);
        return ret;
    }
    /*! Return this.T * a */
    virtual RVector transMult(const RVector & b) const {
        RVector ret(this->cols(), 0.0);
        this->transMult(b, ret);
        return ret;
    }
    /*! Return this.T * a */
    virtual CVector transMult(const CVector & b) const {
        CVector ret(this->cols(), 0.0);
        this->transMult(b, ret);
        return ret;
    }

    virtual RVector mult(const RVector & b, Index startI, Index endI) const {
        log(Warning, "no RVector mult(const RVector & b, Index startI, Index endI) implemented for: ", typeid(*this).name());
        return RVector(rows());
    }
    virtual CVector mult(const CVector & b, Index startI, Index endI) const {
        log(Warning, "no CVector mult(const CVector & b, Index startI, Index endI) implemented for: ", typeid(*this).name());
        return CVector(rows());
    }

/*
    virtual void setCol(Index col, const RVector & v) const {
        THROW_TO_IMPL
    }

    virtual void setCol(Index col, const CVector & v) const {
        THROW_TO_IMPL
    }*/

    /*! Save this matrix into the file filename given. */
    virtual void save(const std::string & filename) const {
        log(Warning, "no save(const std::string & filename) implemented for: ", typeid(*this).name());
    }

protected:
    bool verbose_;
    Index _rows;
    Index _cols;
};

class DLLEXPORT SparseMatrixBase : public MatrixBase {
public:
    SparseMatrixBase(bool verbose=false)
        : MatrixBase(verbose) {}

    /*! Default destructor. */
    virtual ~SparseMatrixBase(){}

    /*! Return entity rtti value. */
    virtual uint rtti() const { return GIMLI_SPARSEMATRIXBASE_RTTI; }

    // Virtual calls are expensive so use them with care
    #define DEFINE_ADDS(A_TYPE) \
    virtual void add(const ElementMatrix < double > & A, const A_TYPE & f, \
                     const double & scale=1.0){ \
        THROW_TO_IMPL \
    } \
    void addS(const ElementMatrix < double > & A, const A_TYPE & f, \
              const double & scale=1.0){ \
        THROW_TO_IMPL \
    } \

    DEFINE_ADDS(double)
    DEFINE_ADDS(RVector3)
    DEFINE_ADDS(RSmallMatrix)
    #undef DEFINE_ADDS

    // virtual void add(const ElementMatrix < double > & A, double scale=1.0, bool neg=false){
    //     THROW_TO_IMPL
    // }
    virtual void setVal(Index row, Index Col, const double & val){
        THROW_TO_IMPL
    }
    virtual void addVal(Index row, Index Col, const double & val){
        THROW_TO_IMPL
    }

protected:
};

//! Identity matrix: derived from matrixBase
class DLLEXPORT IdentityMatrix : public MatrixBase {
public:
    /*! Default constructor (empty matrix). */
    IdentityMatrix()
        : MatrixBase(), val_(0.0){}

    /*! Constructor with number of rows/colums. */
    IdentityMatrix(Index rows, double val=1.0)
        : MatrixBase(rows, rows), val_(val){
    }

    /*! Default destructor. */
    virtual ~IdentityMatrix(){}

    /*! Return this * a  */
    virtual RVector mult(const RVector & a) const {
        if (a.size() != _rows) {
            throwLengthError(WHERE_AM_I + " vector/matrix lengths do not match " +
                                  str(_rows) + " " + str(a.size()));
        }
        return a * val_;
    }

    /*! Return this.T * a */
    virtual RVector transMult(const RVector & a) const {
        if (a.size() != _rows) {
            throwLengthError(WHERE_AM_I + " matrix/vector lengths do not match " +
                                 str(a.size()) + " " + str(_rows));
        }
        return a * val_;
    }

protected:
    double val_;
};

Index DLLEXPORT cblasCount(bool reset=false);
double DLLEXPORT cblasSumTime(bool reset=false);
double DLLEXPORT cblasMinTime(bool reset=false);

//##############################################################################
// DenseMatrix
//##############################################################################

//! Simple row-ordered dense matrix based on continuous memory block.
/*! Simple row-ordered dense matrix based on continuous memory block. */
template < class ValueType > class DLLEXPORT DenseMatrix : public MatrixBase {
public:
    /*! Return entity rtti value. */
    virtual uint rtti() const { return GIMLI_DENSE_MATRIX_RTTI; }

    DenseMatrix()
        : MatrixBase(), _rowView(0){
        // __MS("mat(0)", this);
        this->resize(0, 0);
    }
    /*!Create Densematrix with specified dimensions.*/
    // DenseMatrix(Index rows)
    //     : MatrixBase(), _rowView(0){
    //     this->resize(rows, 0);
    //     log(Error, "Densematrix need rows and cols in constructor.");
    // }

    /*!Create Densematrix with specified dimensions.*/
    DenseMatrix(Index rows, Index cols)
        : MatrixBase(), _rowView(0){
        // __MS("mat", this);
        // __MS("mat(i,i): ", _data.use_count())
        this->resize(rows, cols);
    }
    /*!Create Densematrix with specified dimensions and copy content
    from data.*/
    DenseMatrix(Index rows, Index cols, ValueType * data)
        : MatrixBase(), _rowView(0){
        this->resize(rows, cols);
        std::memcpy(_data.get(), data, sizeof(ValueType)*length());
    }

    DenseMatrix(const DenseMatrix < ValueType > & mat)
        : MatrixBase(), _rowView(0) {
        // __MS("mat(C)", this, &mat);
        copy_(mat);
    }
    DenseMatrix(const Matrix < ValueType > & S)
        : MatrixBase(), _rowView(0){
        this->resize(S.rows(), S.cols());
        for (Index i = 0; i < S.rows(); i ++ ){
            for (Index j = 0; j < S.cols(); j ++ ){
                this->setVal(i, j, S[i][j]);
            }
        }
    }

    virtual ~DenseMatrix(){
        // __MS("~mat", this);
        free_();
    }
    DenseMatrix < ValueType > & operator=(const DenseMatrix < ValueType > & mat) {
        // __MS("mat=", this, &mat);
        if (this != & mat){
            copy_(mat);
        } return *this;
    }

    // /*!Return write access \ref Vector of borrowed memory for the ith row.*/
    // inline Vector< ValueType > & operator[](Index i) {
    //     return row(i);
    // }
    // /*!Return \ref Vector of borrowed memory for the ith row.*/
    // inline const Vector< ValueType > & operator[](Index i) const {
    //     return row(i);
    // }
    /*!Return write access \ref Vector of borrowed memory for the ith row.*/
    inline Vector< ValueType > operator[](Index i) {
        return row(i);
    }
    /*!Return \ref Vector of borrowed memory for the ith row.*/
    inline Vector< ValueType > operator[](Index i) const {
        return row(i);
    }

    /*! Read only access to matrix element i,j. */
    inline const ValueType & operator ()(Index i, Index j) const {
        // assert(i < this->_rows && j < this->_cols);

        if (i >= this->_rows || j >= this->_cols) {
            throwLengthError(WHERE_AM_I + " ASSERT_LOWER2: ");
            // weirdly this costs 80% runtime performance even if not called.
            // throwLengthError(WHERE_AM_I + " ASSERT_LOWER2: " + str(i) + " < "  + str(this->_rows) + " or " + str(j) + " < "  + str(this->_cols));
        }
        // ASSERT_LOWER2(i, this->_rows, j, this->_cols)
        return _data[this->_cols * i + j];
    }

    /*! Write access to matrix element i,j. */
    inline ValueType & operator ()(Index i, Index j) {
        // assert(i < this->_rows && j < this->_cols);

        if (i >= this->_rows || j >= this->_cols) {
            throwLengthError(WHERE_AM_I + " ASSERT_LOWER2: ");
            // weirdly this costs 80% runtime performance even if not called.
            // throwLengthError(WHERE_AM_I + " ASSERT_LOWER2: " + str(i) + " < "  + str(this->_rows) + " or " + str(j) + " < "  + str(this->_cols));
        }
        return _data[this->_cols * i + j];
    }

    inline void setVal(Index i, Index j, const ValueType & v){
        this->operator()(i,j) = v;
        // _data[this->_cols * i + j]  = v;
    }
    inline void addVal(Index i, Index j, const ValueType & v){
        this->operator()(i,j) += v;
        // _data[this->_cols * i + j] += v;
    }

    inline ValueType * pData(){ return _data.get(); }
    inline const ValueType * pData() const { return _data.get(); }
    // #ifndef PYGIMLI_CAST
    // #endif
    std::shared_ptr< ValueType [] > & data(){
        return _data;
    }

    /*! For template compatablity. Does nothing but return data buffer.*/
    ValueType * toData(ValueType * target, Index size=0) const {
        return _data.get();
    }
    /*! For template compatablity. Does nothing if src is equal data or perform memcopy. */
    void fromData(ValueType * src, Index m, Index n){
        if (src != _data.get()){
            __MS("should I be here? Check!")
            ASSERT_THIS_SIZE(m)
            ASSERT_EQUAL(n, this->_cols)
            std::memcpy(src, _data.get(), sizeof(ValueType) * length());
        }
    }
    /*! Return all values as \ref Vector of size \ref length*/
    Vector< ValueType > values() const {
        return Vector< ValueType >(_data.get(), this->length());
    }

    // /*! Return read only view to row i*/
    // const Vector< ValueType > & row(Index i) const;

    // /*! Return view to row i*/
    // Vector< ValueType > & row(Index i);

    /*! Return read only view to row i*/
    Vector< ValueType > row(Index i) const {
#ifndef PYGIMLI_CAST
        ASSERT_THIS_SIZE(i)
        // return Vector< ValueType >(&_data[this->_cols * i], this->_cols);
        return Vector< ValueType >(this->_cols, _data, this->_cols * i);
#endif
        return Vector < ValueType >(0);
    }

    /*! Return view to row i*/
    Vector< ValueType > row(Index i){
#ifndef PYGIMLI_CAST
        ASSERT_THIS_SIZE(i)
        // return Vector< ValueType >(&_data[this->_cols * i], this->_cols);
        return Vector< ValueType >(this->_cols, _data, this->_cols * i);
#endif
        return Vector < ValueType >(0);
    }

    /*! Return read only view to row i*/
    Vector< ValueType > rowView(Index i) const{
#ifndef PYGIMLI_CAST
        ASSERT_THIS_SIZE(i)
        return Vector< ValueType >(this->_cols, _data, this->_cols * i);
#endif
    return Vector < ValueType >(0);
    }

    /*! Return view to row i*/
    Vector< ValueType > rowView(Index i){
#ifndef PYGIMLI_CAST
        ASSERT_THIS_SIZE(i)
        return Vector< ValueType >(this->_cols, _data, this->_cols * i);
#endif
        return Vector < ValueType >(0);
    }

    inline void setRow(Index i, const Vector< ValueType > & r) {
        row(i).assign(r);
    }
    inline void setCol(Index j, const Vector< ValueType > & r) {
        ASSERT_LOWER(j, this->_cols)
        for (Index i = 0; i < this->_rows; i ++){
            _data[i*this->_cols + j] = r[i];
        }
    }
    inline void cleanCols(const IndexArray & c) {
        for (Index i = 0; i < this->_rows; i ++){
            for (auto j: c){
                ASSERT_LOWER(j, this->_cols)
                _data[i*this->_cols + j] = 0;
            }
        }
    }

    inline Vector< ValueType > col(Index c) const {
        ASSERT_LOWER(c, this->_cols)
        Vector < ValueType > ret(this->_rows);
        for (Index i = 0; i < this->_rows; i ++){
            ret[i] = _data[i*this->_cols + c];
        }
        return ret;
    }
    inline Vector< ValueType > back() {
        ASSERT_THIS_SIZE(1)
        return row(_rows - 1);
    }
    void push_back(const Vector< ValueType > & vec) {
        ASSERT_VEC_SIZE(vec, this->_cols)
        log(Warning, "Efficency .. push_back for dense matrix not recommanded");
        __MS("Check Refcounter!")

        if (_data.use_count() > 2){
            log(Error, "Cannot push_back on data that has been borrowed.");
        }
        Index oldLength = length();

        ValueType *tmp = new ValueType[oldLength];
        std::memcpy(tmp, &_data[0], oldLength);
        this->resize(_rows+1, _cols);
        std::memcpy(&_data[0], tmp, oldLength);

        // std::shared_ptr< ValueType [] > d(new ValueType[(_rows+1)*_cols]);
        // std::memcpy(&d[0], &_data[0], length());
        // this->_rows ++;
        // __M
        row(this->_rows-1).assign(vec);
        delete [] tmp;
        // __M
        // _data = d;
        // __M
    }

    #define DEFINE_UNARY_MOD_OPERATOR__(OP, NAME) \
    inline DenseMatrix < ValueType > & operator OP##=(const DenseMatrix < ValueType>&A){\
        if (A.rows() == this->rows() && A.cols() == this->cols()) { \
            for (Index i = 0; i < length(); i ++) {_data[i] OP##= A.pData()[i];} \
            return *this;\
        }\
        if (A.rows() == 1 && A.cols() == this->cols()) { \
            for (Index i = 0; i < this->size(); i ++) {row(i) OP##= A[0]; }\
            return *this;\
        } \
        if (A.cols() == 1 && A.rows() == this->rows()) { \
            for (Index i = 0; i < this->size(); i ++) {row(i) OP##= A[i][0];} \
            return *this;\
        } \
        if (this->rows() == 1 && this->cols() == A.cols()){ \
            Vector < ValueType > tmp(&_data[0], this->cols()); /*!need deepcopy */\
            this->resize(A.rows(), this->cols()); \
            for (Index i = 0; i < this->size(); i ++){ \
                this->setRow(i, A[i] OP tmp);} \
            return *this;\
        } \
        if (this->cols() == 1 && this->rows() == A.rows()){ \
            Vector < ValueType > tmp(this->col(0)); \
            this->resize(this->rows(), A.cols()); \
            for (Index i = 0; i < this->size(); i ++){ \
                this->setRow(i, A[i] OP tmp[i]);} \
            return *this;\
        } \
        throwLengthError(WHERE_AM_I + " Cannot operate on mishaped matrices. "+\
            "(" + str(this->rows()) + "," + str(this->cols()) + ") " + "OP" + \
            " (" + str(A.rows()) + "," + str(A.cols()) + ")" );\
        return *this;\
    }\
    inline DenseMatrix < ValueType > & operator OP##= (const ValueType & val) { \
      for (Index i = 0; i < length(); i ++){_data.get()[i] OP##= val;}return*this;}\
    inline DenseMatrix < ValueType > & operator OP##= (const Vector < ValueType > & val) { \
      for (Index i = 0; i < this->size(); i ++){row(i) OP##= val;}return*this;}\

    DEFINE_UNARY_MOD_OPERATOR__(+, PLUS)
    DEFINE_UNARY_MOD_OPERATOR__(-, MINUS)
    DEFINE_UNARY_MOD_OPERATOR__(/, DIVID)
    DEFINE_UNARY_MOD_OPERATOR__(*, MULT)

    #undef DEFINE_UNARY_MOD_OPERATOR__

    /*! A += a.T*/
    DenseMatrix < ValueType > & transAdd(const DenseMatrix < ValueType > & a);

    /*! Multiplication c = alpha * (A*b) + beta * c. */
    inline void mult(const Vector < ValueType > & b,
                      Vector < ValueType >& c,
                      const ValueType & alpha=1.0,
                      const ValueType & beta=0.0,
                      Index bOff=0, Index cOff=0) const {
        return GIMLI::mult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Multiplication c = alpha * (A.T*b) + beta * c. */
    inline void transMult(const Vector < ValueType > & b,
                          Vector < ValueType > & c,
                          const ValueType & alpha=1.0,
                          const ValueType & beta=0.0,
                          Index bOff=0, Index cOff=0) const {
        return GIMLI::transMult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Multiplication c = alpha * (A*b) + beta * c.
    Same as mult but slightly faster to call from python because of unique name. */
    inline void multMV(const Vector < ValueType > & b,
                       Vector < ValueType >& c,
                       const ValueType & alpha=1.0,
                       const ValueType & beta=0.0,
                       Index bOff=0, Index cOff=0) const {
        return GIMLI::mult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Multiplication c = alpha * (A.T*b) + beta * c.
    Same like transMult but slightly faster to call from python because of unique name. */
    inline void transMultMV(const Vector < ValueType > & b,
                            Vector < ValueType > & c,
                            const ValueType & alpha=1.0,
                            const ValueType & beta=0.0,
                            Index bOff=0, Index cOff=0) const {
        return GIMLI::transMult(*this, b, c, alpha, beta, bOff, cOff);
    }

    /*! Return this * a  */
    inline Vector < ValueType > mult(const Vector < ValueType > & b) const {
        Vector < ValueType > ret(this->rows(), 0.0);
        this->mult(b, ret);
        return ret;
    }
    /*! Return this.T * a */
    inline Vector < ValueType > transMult(const Vector < ValueType > & b) const {
        Vector < ValueType > ret(this->cols(), 0.0);
        this->transMult(b, ret);
        return ret;
    }
    /*! Multiplication C = alpha * (A*B) + beta * C. */
    inline void mult(const DenseMatrix < ValueType > & B,
                     DenseMatrix < ValueType > & C,
                     const ValueType & alpha = 1.0,
                     const ValueType & beta = 0.0) const {
        GIMLI::mult((*this), B, C, alpha, beta);
    }
    /*! Multiplication C = alpha * (A.T*B) + beta * C. */
    inline void transMult(const DenseMatrix < ValueType > & B,
                          DenseMatrix < ValueType > & C,
                          const ValueType & alpha = 1.0,
                          const ValueType & beta = 0.0) const {
        GIMLI::transMult((*this), B, C, alpha, beta);
    }

    /*! Multiplication C = alpha * (A*B) + beta * C.
    Same like mult but slightly faster to call from python because of unique name. */
    inline void multMM(const DenseMatrix < ValueType > & B,
                     DenseMatrix < ValueType > & C,
                     const ValueType & alpha = 1.0,
                     const ValueType & beta = 0.0) const {
        GIMLI::mult((*this), B, C, alpha, beta);
    }


    /*! Multiplication C = alpha * (A.T*B) + beta * C. Same like transMult but slightly faster to call from python because of unique name. */
    inline void transMultMM(const DenseMatrix < ValueType > & B,
                          DenseMatrix < ValueType > & C,
                          const ValueType & alpha = 1.0,
                          const ValueType & beta = 0.0) const {
        GIMLI::transMult((*this), B, C, alpha, beta);
    }

    /*! Return sum of all values.*/
    ValueType sum() const {
        return std::accumulate(&_data[0], &_data[_rows*_cols], ValueType(0));
    }

    /*! Round all values of this matrix to tol.*/
    void round(const ValueType & tol);

    /*! Resize the matrix to rows x cols. */
    virtual void resize(Index rows, Index cols){
        // content not more valid anymore
        allocate_(rows, cols);
    }

    /*! Fill Vector with 0.0. Don't change size.*/
    inline void clean() {
        std::memset(_data.get(), '\0', sizeof(ValueType) * _rows*_cols);
    }

    /*! Clear the matrix and free memory. */
    inline void clear() { this->free_(); }

    inline Index length() const {return this->_rows * this->_cols;}

    std::shared_ptr< ValueType [] > _data;

    mutable Vector < ValueType > * _rowView;
protected:
    void copy_(const DenseMatrix< ValueType > & mat){
        resize(mat.rows(), mat.cols());
        std::memcpy(_data.get(), mat.pData(), sizeof(ValueType) * length());
    }

    void allocate_(Index rows, Index cols){
        if (rows * cols > _rows * _cols){

            if (_data.use_count() > 2){
                __MS(_data.use_count())
               throwError("Matrix data are in use and can't be new allocated.");
            }

            // use_count == 2 because default rowView use one additional
            if (_data.use_count() == 2 || _data != nullptr){
                free_();
            }

            _data = std::shared_ptr< ValueType [] >(new ValueType[rows*cols]);

            std::memset(_data.get(), '\0', sizeof(ValueType) * rows*cols);

        }

        if (cols != _cols || _rowView == nullptr){
            if (_rowView != nullptr){
                delete _rowView;
            }
            if (_data != nullptr){
#ifndef PYGIMLI_CAST
                _rowView = new Vector< ValueType >(cols, _data, 0);
#endif
            }
        }

        _rows = rows;
        _cols = cols;
    }

    void free_(){
        _rows = 0;
        _cols = 0;

        if (_rowView != nullptr){
            //** --> _data.use_count() == 2
            delete _rowView;
            //** --> _data.use_count() == 1
            _rowView = 0;
        }

        if (_data.use_count() > 1){
            __MS(_data.use_count())
            log(Error, "Matrix data are in use and can't be deallocated.", this);
        }
        if (_data.use_count() == 1){
            _data.reset();
        }
    }

};

#define DEFINE_BINARY_OPERATOR__(OP, NAME) \
template < class ValueType > \
DenseMatrix < ValueType > operator OP (const DenseMatrix < ValueType > & A, \
                                       const DenseMatrix < ValueType > & B) { \
DenseMatrix < ValueType > tmp(A); \
return tmp OP##= B; } \
template < class ValueType > \
DenseMatrix < ValueType > operator OP (const DenseMatrix < ValueType > & A, \
                                       const ValueType & v) { \
DenseMatrix < ValueType > tmp(A); \
return tmp OP##= v; }

DEFINE_BINARY_OPERATOR__(+, PLUS)
DEFINE_BINARY_OPERATOR__(-, MINUS)
DEFINE_BINARY_OPERATOR__(/, DIVID)
DEFINE_BINARY_OPERATOR__(*, MULT)

#undef DEFINE_BINARY_OPERATOR__

// template <> DLLEXPORT const Vector< double > &
// DenseMatrix< double >::row(Index i) const;
// template <> DLLEXPORT Vector< double > &
// DenseMatrix< double >::row(Index i);

// template <> DLLEXPORT const Vector< Complex > &
// DenseMatrix< Complex >::row(Index i) const;
// template <> DLLEXPORT Vector< Complex > &
// DenseMatrix< Complex >::row(Index i);

template <> DLLEXPORT void
DenseMatrix< double >::round(const double & v);
template <> DLLEXPORT void
DenseMatrix< Complex >::round(const Complex & v);

template <> DLLEXPORT DenseMatrix<double> &
DenseMatrix<double>::transAdd(const DenseMatrix < double > & a);
template <> DLLEXPORT DenseMatrix<Complex> &
DenseMatrix<Complex>::transAdd(const DenseMatrix < Complex > & a);


//##############################################################################
// Matrix
//##############################################################################

//! Simple row-based dense matrix based on \ref Vector
/*! Simple row-based dense matrix based on \ref Vector */
template < class ValueType > class DLLEXPORT Matrix : public MatrixBase {
public:
    /*! Constructs an empty matrix with the dimension rows x cols. Content of the matrix is zero. */
    Matrix()
        : MatrixBase() {
        resize(0, 0);
    }
    // Matrix(Index rows)
    //     : MatrixBase() {
    //     resize(rows, 0);
    // }
    // no default arg here .. pygimli@win64 linker bug
    Matrix(Index rows, Index cols)
        : MatrixBase() {
        resize(rows, cols);
    }
    Matrix(Index rows, Index cols, ValueType *src)
        : MatrixBase() {
        fromData(src, rows, cols);
    }
    // Matrix(const Vector< ValueType > & r)
    //     : MatrixBase() {
    //     resize(1, r.size());
    //     mat_[0] = r;
    // }
    /*! Copy constructor */

    Matrix(const std::vector < Vector< ValueType > > & mat)
        : MatrixBase(){ copy_(mat); }

    /*! Constructor, read matrix from file see \ref load(Matrix < ValueType > & A, const std::string & filename). */
    Matrix(const std::string & fileName)
        : MatrixBase() { this->load(fileName); }

    /*! Copy constructor */
    Matrix(const Matrix < ValueType > & mat)
        : MatrixBase() { copy_(mat); }

    /*! Assignment operator */
    Matrix < ValueType > & operator = (const Matrix< ValueType > & mat){
        if (this != & mat){
            copy_(mat);
        } return *this;
    }

    /*! Destruct matrix and free memory. */
    virtual ~Matrix(){}

    /*! Force the copy of the matrix entries. */
    inline void copy(const Matrix < ValueType > & mat){ copy_(mat); }

    /*! Load content of file. */
    inline bool load(const std::string & fileName){
        return GIMLI::load(*this, fileName);
    }

    /*! Return entity rtti value. */
    virtual uint rtti() const { return GIMLI_MATRIX_RTTI; }

    #define DEFINE_UNARY_MOD_OPERATOR__(OP, NAME) \
    inline Matrix < ValueType > & operator OP##=(const Matrix < ValueType>&A){\
        if (A.rows() == this->rows() && A.cols() == this->cols()) { \
            for (Index i = 0; i < mat_.size(); i ++) {mat_[i] OP##= A[i];} \
            return *this;\
        }\
        if (A.rows() == 1 && A.cols() == this->cols()) { \
            for (Index i = 0; i < mat_.size(); i ++) {mat_[i] OP##= A[0]; }\
            return *this;\
        } \
        if (A.cols() == 1 && A.rows() == this->rows()) { \
            for (Index i = 0; i < mat_.size(); i ++) {mat_[i] OP##= A[i][0];} \
            return *this;\
        } \
        if (this->rows() == 1 && this->cols() == A.cols()){ \
            Vector < ValueType > tmp(this->row(0)); \
            this->resize(A.rows(), this->cols()); \
            for (Index i = 0; i < mat_.size(); i ++){ \
                this->setRow(i, A[i] OP tmp);} \
            return *this;\
        } \
        if (this->cols() == 1 && this->rows() == A.rows()){ \
            Vector < ValueType > tmp(this->col(0)); \
            this->resize(this->rows(), A.cols()); \
            for (Index i = 0; i < mat_.size(); i ++){ \
                this->setRow(i, A[i] OP tmp[i]);} \
            return *this;\
        } \
        throwLengthError(WHERE_AM_I + " Cannot operate on mishaped matrices. "+\
            "(" + str(this->rows()) + "," + str(this->cols()) + ") " + "OP" + \
            " (" + str(A.rows()) + "," + str(A.cols()) + ")" );\
        return *this;\
    }\
    inline Matrix < ValueType > & operator OP##= (const ValueType & val) { \
      for (Index i = 0; i < mat_.size(); i ++){mat_[i] OP##= val;}return*this;}\
    inline Matrix < ValueType > & operator OP##= (const Vector < ValueType > & val) { \
      for (Index i = 0; i < mat_.size(); i ++){mat_[i] OP##= val;}return*this;}\

    DEFINE_UNARY_MOD_OPERATOR__(+, PLUS)
    DEFINE_UNARY_MOD_OPERATOR__(-, MINUS)
    DEFINE_UNARY_MOD_OPERATOR__(/, DIVID)
    DEFINE_UNARY_MOD_OPERATOR__(*, MULT)

    #undef DEFINE_UNARY_MOD_OPERATOR__

//     Index col = cols();
//         for (Index i = 0; i < mat_.size(); i ++) {
//             ValueType * Aj = &mat_[i][0];
//             ValueType * Aje = &mat_[i][col];
//             for (; Aj != Aje;) *Aj++ OP##= val;
//         }   return *this; }

    /*! Index operator for write operations without boundary check. */
    Vector< ValueType > & operator [] (Index i) {
        return rowRef(i);
    }

    /*! Read only C style index operator, without boundary check. */
    const Vector< ValueType > & operator [] (Index i) const {
        return row(i);
    }

    /*! Read only access to matrix element i,j. */
    inline const ValueType & operator ()(Index i, Index j) const {
        // ASSERT_THIS_SIZE(i)
        // ASSERT_LOWER(j, this->_cols)
        return mat_[i][j];
    }

    /*! Write access to matrix element i,j. */
    inline ValueType & operator ()(Index i, Index j) {
        // ASSERT_THIS_SIZE(i) // doubles runtime
        // ASSERT_LOWER(j, this->_cols)
        return mat_[i][j];
    }

    /*! Read only access to matrix row i. */
    inline const Vector< ValueType > & operator ()(Index i) const {
        return mat_[i];}

    /*! Write access to matrix row i */
    inline Vector< ValueType > & operator ()(Index i) {
        return mat_[i];}

    /*! Implicite type converter. */
    template < class T > operator Matrix< T >(){
        Matrix< T > f(this->rows());
        for (uint i = 0; i < this->rows(); i ++){ f[i] = Vector < T >(mat_[i]); }
        return f;
    }

    /*! Resize the matrix to rows x cols. */
    virtual void resize(Index rows, Index cols){ allocate_(rows, cols); }

    /*! Fill Vector with 0.0. Don't change size.*/
    inline void clean() {
        for (Index i = 0; i < mat_.size(); i ++) mat_[i].clean();
    }

    /*! Clear the matrix and free memory. */
    inline void clear() {
        mat_.clear();
        _cols = 0; _rows = 0;
    }

    inline Index length() const {return rows() * cols();}

    /*! Set a value. Throws out of range exception if index check fails. */
    // inline void setRow(const Vector < ValueType > & val, Index i) {
    //     log(Warning, "deprecated use setRow(i, val)");
    //     ASSERT_THIS_SIZE(i)
    //     mat_[i] = val;
    // }
    inline void setRow(Index i, const Vector < ValueType > & val) {
        ASSERT_THIS_SIZE(i)
        mat_[i] = val;
        this->_cols = max(this->_cols, val.size());
    }

    /*! Set a value. Throws out of range exception if index check fails. */
    // inline void setVal(const Vector < ValueType > & val, Index i) {
    //     log(Warning, "deprecated, use setVal(i, val)");
    //     return setRow(i, val);
    // }
    inline void setVal(Index i, const Vector < ValueType > & val) {
        return this->setRow(i, val);
    }
    inline void setVal(Index i, Index j, const ValueType & val) {
        this->rowRef(i).setVal(val, j);
    }
    inline void addVal(Index i, Index j, const ValueType & val) {
        this->rowRef(i).addVal(val, j);
    }

    /*! Return reference to row. Used for pygimli. */
    inline Vector < ValueType > & row(Index i) {
        ASSERT_THIS_SIZE(i)
        return mat_[i];
    }
    /*! Readonly getter. */
    inline const Vector < ValueType > & row(Index i) const {
        ASSERT_THIS_SIZE(i)
        return mat_[i];
    }
    /*! Return reference to row. Used for pygimli. */
    inline Vector < ValueType > & rowRef(Index i) {
        ASSERT_THIS_SIZE(i)
        return mat_[i];
    }

    /*! Readonly column entry of matrix, with boundary check. Probably slow.*/
    virtual const Vector< ValueType > col(Index i) const {
        //__M
        if (i < 0 || i > this->cols()-1) {
            throwLengthError(WHERE_AM_I + " col bounds out of range " +
                                str(i) + " " + str(this->cols())) ;
        }
        Vector < ValueType > col(this->rows());
        for (Index j = 0, jmax = rows(); j < jmax; j ++) col[j] = mat_[j][i];
        return col;
    }


    /*! Add another row vector at the end. */
    inline void push_back(const Vector < ValueType > & vec) {
        //** push_back reallocates mem and copies content
        this->_cols = max(this->_cols, vec.size());
        mat_.push_back(vec);
        rowFlag_.resize(rowFlag_.size() + 1);
        this->_rows = mat_.size();

        // don't use vec here .. if vec was content of this matrix is isnt
        // valid anymore after any reallocation

        // if (vec.size() > 10){
        //     print(&vec);
        //     for (auto &v: mat_){
        //         print(&v, v);
        //     }
        // }

    }

    /*! Return last row vector. */
    inline Vector< ValueType > & back() { return mat_.back(); }

    /*! Set one specific column */
    inline void setCol(Index col, const Vector < ValueType > & v){
        if (col < 0 || col > this->cols()-1) {
            throwLengthError(WHERE_AM_I + " col bounds out of range " +
                                str(col) + " " + str(this->cols())) ;
        }
        if (v.size() > this->rows()) {
            throwLengthError(WHERE_AM_I + " rows bounds out of range " +
                                str(v.size()) + " " + str(this->rows())) ;
        }
        for (Index i = 0; i < v.size(); i ++) mat_[i][col] = v[i];
    }

    /*! Add one specific column */
    inline void addCol(Index col, const Vector < ValueType > & v){
        if (col < 0 || col > this->cols()-1) {
            throwLengthError(WHERE_AM_I + " col bounds out of range " +
                                str(col) + " " + str(this->cols())) ;
        }
        if (v.size() > this->rows()) {
            throwLengthError(WHERE_AM_I + " rows bounds out of range " +
                                str(v.size()) + " " + str(this->rows())) ;
        }
        for (Index i = 0; i < v.size(); i ++) mat_[i][col] += v[i];
    }

    /*! Return reference to row flag vector. Maybee you can check if the rows are valid. Size is set automatic to the amount of rows. */
    BVector & rowFlag(){ return rowFlag_; }

    /*! A += a*/
    inline Matrix < ValueType > & add(const Matrix < ValueType > & a){
        return (*this)+=a;
    }

    /*! A += a.T*/
    Matrix < ValueType > & transAdd(const Matrix < ValueType > & a);

    /*! Multiplication c = alpha * (A*b) + beta * c. */
    inline void mult(const Vector < ValueType > & b,
                      Vector < ValueType >& c,
                      const ValueType & alpha=1.0,
                      const ValueType & beta=0.0,
                      Index bOff=0, Index cOff=0) const {
        return GIMLI::mult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Multiplication c = alpha * (A.T*b) + beta * c. */
    inline void transMult(const Vector < ValueType > & b,
                          Vector < ValueType > & c,
                          const ValueType & alpha=1.0,
                          const ValueType & beta=0.0,
                          Index bOff=0, Index cOff=0) const {
        return GIMLI::transMult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Multiplication c = alpha * (A*b) + beta * c. Same like mult but slightly faster to call from python because of unique name. */
    inline void multMV(const Vector < ValueType > & b,
                      Vector < ValueType >& c,
                      const ValueType & alpha=1.0,
                      const ValueType & beta=0.0,
                      Index bOff=0, Index cOff=0) const {
        return GIMLI::mult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Multiplication c = alpha * (A.T*b) + beta * c.  Same like transMult but slightly faster to call from python because of unique name. */
    inline void transMultMV(const Vector < ValueType > & b,
                          Vector < ValueType > & c,
                          const ValueType & alpha=1.0,
                          const ValueType & beta=0.0,
                          Index bOff=0, Index cOff=0) const {
        return GIMLI::transMult(*this, b, c, alpha, beta, bOff, cOff);
    }
    /*! Return this * a  */
    inline Vector < ValueType > mult(const Vector < ValueType > & b) const {
        Vector < ValueType > ret(this->rows(), 0.0);
        this->mult(b, ret);
        return ret;
    }
    /*! Return this.T * a */
    inline Vector < ValueType > transMult(const Vector < ValueType > & b) const {
        Vector < ValueType > ret(this->cols(), 0.0);
        this->transMult(b, ret);
        return ret;
    }

    /*! Multiplication C = alpha * (A*B) + beta * C. */
    inline void mult(const Matrix < ValueType > & B,
                     Matrix < ValueType > & C,
                     const ValueType & alpha = 1.0,
                     const ValueType & beta = 0.0) const {
        GIMLI::mult((*this), B, C, alpha, beta);
    }
    /*! Multiplication C = alpha * (A.T*B) + beta * C. */
    inline void transMult(const Matrix < ValueType > & B,
                          Matrix < ValueType > & C,
                          const ValueType & alpha = 1.0,
                          const ValueType & beta = 0.0) const {
        GIMLI::transMult((*this), B, C, alpha, beta);
    }
    /*! Multiplication C = alpha * (A*B) + beta * C. Same like mult but slightly faster to call from python because of unique name. */
    inline void multMM(const Matrix < ValueType > & B,
                     Matrix < ValueType > & C,
                     const ValueType & alpha = 1.0,
                     const ValueType & beta = 0.0) const {
        GIMLI::mult((*this), B, C, alpha, beta);
    }
    /*! Multiplication C = alpha * (A.T*B) + beta * C. Same like transMult but slightly faster to call from python because of unique name. */
    inline void transMultMM(const Matrix < ValueType > & B,
                          Matrix < ValueType > & C,
                          const ValueType & alpha = 1.0,
                          const ValueType & beta = 0.0) const {
        GIMLI::transMult((*this), B, C, alpha, beta);
    }

    /*! Save matrix to file. */
    virtual void save(const std::string & filename) const {
        saveMatrix(*this, filename);
    }

    /*! Return sum of all values.*/
    ValueType sum() const {
        ValueType ret(0);
        for (Index i = 0; i < mat_.size(); i ++) ret += GIMLI::sum(mat_[i]);
        return ret;
    }

    /*! Round each matrix element to a given tolerance. */
    void round(const ValueType & tolerance){
        for (Index i = 0; i < mat_.size(); i ++) mat_[i].round(tolerance);
        // ??? std::for_each(mat_.begin, mat_.end, boost::bind(&Vector< ValueType >::round, tolerance));
    }

    /*!Copy matrix content into buf of length size. If size is not sufficient,  a new buffer is created which need to be deleted manually. */
    ValueType * toData(ValueType * buf, Index size=0) const{
        ValueType * b_;

        if (size < length()){
            b_ = new ValueType[length()];
        } else {
            b_ = buf;
        }
        Index N = sizeof(ValueType) * this->cols();

        for (Index i = 0; i < mat_.size(); i ++) {
            std::memcpy(&b_[i*this->cols()], &mat_[i][0], N);
        }
        return b_;
    }
    void fromData(ValueType * src, Index m, Index n){
        this->resize(m, n);
        Index N = sizeof(ValueType) * n;
        for (Index i = 0; i < m; i ++) {
            std::memcpy(&mat_[i][0], &src[i*n], N);
        }
    }

    Vector< ValueType > values() const {
        Vector< ValueType > b_(this->rows()*this->cols());
        Index N = sizeof(ValueType) * this->cols();

        for (Index i = 0; i < mat_.size(); i ++) {
            std::memcpy(&b_[i*this->cols()], &mat_[i][0], N);
        }
        return b_;
    }

    ValueType trace() const {
        if (this->_cols != this->_rows){
            log(Error, "trace not defined for non quadratic matrix",
            this->_cols, this->_rows);
        };
        ValueType ret(0) ;
        for (Index i = 0; i < mat_.size(); i ++) {
            ret += mat_[i][i];
        }
        return ret;
    }
	std::vector < Vector< ValueType > > mat_;

protected:

    void allocate_(Index rows, Index cols){

        this->_cols = cols;
        this->_rows = rows;

        if (mat_.size() != rows) mat_.resize(rows);
        for (Index i = 0; i < mat_.size(); i ++) {
//             __MS(this << " " << &mat_[i] << " "<< cols)
            mat_[i].resize(cols);
//             __MS(&mat_[i] << " " <<mat_[i].size())
        }
        rowFlag_.resize(rows);
    }

    void copy_(const Matrix < ValueType > & mat){
        allocate_(mat.rows(), mat.cols());
        for (Index i = 0; i < mat_.size(); i ++) mat_[i] = mat[i];
    }

    /*! BVector flag(rows) for free use, e.g., check if rows are set valid. */
    BVector rowFlag_;
};

template <> DLLEXPORT Matrix<double> &
Matrix<double>::transAdd(const Matrix < double > & a);
template <> DLLEXPORT Matrix<Complex> &
Matrix<Complex>::transAdd(const Matrix < Complex > & a);


#define DEFINE_BINARY_OPERATOR__(OP, NAME) \
template < class ValueType > \
Matrix < ValueType > operator OP (const Matrix < ValueType > & A, \
                                  const Matrix < ValueType > & B) { \
Matrix < ValueType > tmp(A); \
return tmp OP##= B; } \
template < class ValueType > \
Matrix < ValueType > operator OP (const Matrix < ValueType > & A, \
                                  const ValueType & v) { \
Matrix < ValueType > tmp(A); \
return tmp OP##= v; }

DEFINE_BINARY_OPERATOR__(+, PLUS)
DEFINE_BINARY_OPERATOR__(-, MINUS)
DEFINE_BINARY_OPERATOR__(/, DIVID)
DEFINE_BINARY_OPERATOR__(*, MULT)

#undef DEFINE_BINARY_OPERATOR__

template < class Mat >
bool operator == (const Mat & A, const Mat & B){

    if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
    for (Index i = 0; i < A.rows(); i ++){
        // __MS(i, A.rows(), A[i] != B[i])
        if (A[i] != B[i]) return false;
    }
    return true;
}

template < class ValueType, template < typename > class Mat >
void scaleMatrix(Mat < ValueType >& A,
                 const Vector < ValueType > & l,
                 const Vector < ValueType > & r){
    Index rows = A.rows();
    Index cols = A.cols();
    if (rows != l.size()){
        throwLengthError(WHERE_AM_I + " " + str(rows) + " != " + str(l.size()));
    };
    if (cols != r.size()){
        throwLengthError(WHERE_AM_I + " " + str(cols) + " != " + str(r.size()));
    }

    for (Index i = 0 ; i < rows ; i++) {
        //for (Index j = 0 ; j < cols ; j++) A[i][j] *= (l[i] * r[j]);
        A[i] *= r * l[i];
    }
}

template < class ValueType, template < typename > class Mat >
void rank1Update(Mat < ValueType > & A,
                 const Vector < ValueType > & u,
                 const Vector < ValueType > & v) {
    Index rows = A.rows();
    Index cols = A.cols();
    if (rows != u.size()){
        throwLengthError(WHERE_AM_I + " " + str(rows) + " != " + str(u.size()));
    };

    if (cols != v.size()){
        throwLengthError(WHERE_AM_I + " " + str(cols) + " != " + str(v.size()));
    }

    for (Index i = 0 ; i < rows ; i++) {
        //for (Index j = 0 ; j < ncols ; j++) A[i][j] += (u[i] * v[j]);
        A[i] += v * u[i];
    }
    return;
}

template < class ValueType, template < typename > class Mat >
Mat < ValueType > fliplr(const Mat< ValueType > & m){
    Mat < ValueType > ret(m.rows(), m.cols());
    for (Index i = 0; i < m.rows(); i ++) ret[i] = fliplr(m[i]);
    return ret;
}

template < class ValueType >
Matrix < ValueType > real(const Matrix < std::complex< ValueType > > & cv){
    Matrix < ValueType > v(cv.rows(), cv.cols());
    for (Index i = 0; i < cv.rows(); i ++) v[i] = real(cv[i]);
    return v;
}

template < class ValueType >
Matrix < ValueType > imag(const Matrix < std::complex< ValueType > > & cv){
    Matrix < ValueType > v(cv.rows(), cv.cols());
    for (Index i = 0; i < cv.rows(); i ++) v[i] = imag(cv[i]);
    return v;
}

//********************* MATRIX I/O *******************************

/*! Save matrix into a file (Binary).
    File suffix ($MATRIXBINSUFFIX) will be append if none given.
    Format: rows(uint32) cols(uint32) vals(rows*cols(ValueType))
    If IOFormat == Ascii matrix will be saved in Ascii format, See: \ref saveMatrixRow
*/
template < class ValueType >
bool saveMatrix(const Matrix < ValueType > & A, const std::string & filename, IOFormat format = Binary){
    if (format == Ascii) return saveMatrixRow(A, filename);
    std::string fname(filename);
    if (fname.rfind('.') == std::string::npos) fname += MATRIXBINSUFFIX;

    FILE *file; file = fopen(fname.c_str(), "w+b");

    if (!file){
		std::cerr << fname << ": " << strerror(errno) << " " << errno << std::endl;
        return false;
    }

    uint32 rows = A.rows();
    uint ret = fwrite(& rows, sizeof(uint32), 1, file);
    if (ret == 0) {
        fclose(file);
        return false;
    }
    uint32 cols = A.cols();
    ret = fwrite(& cols, sizeof(uint32), 1, file);

    for (uint i = 0; i < rows; i ++){
        for (uint j = 0; j < cols; j ++){
            ret = fwrite(&A[i][j], sizeof(ValueType), 1, file);
        }
    }
    fclose(file);
    return true;
}

/*! Load matrix from a single or multiple files (Binary).
    File suffix (\ref MATRIXBINSUFFIX, ".matrix", ".mat") given or not -- loads single datafile, else try to load matrix from multiple binary vector files.
    Single format: see \ref save(const Matrix < ValueType > & A, const std::string & filename)
*/
template < class ValueType >
bool load(Matrix < ValueType > & A, const std::string & filename){

    //!* First check if filename suffix is ".matrix", ".mat", \ref MATRIXBINSUFFIX;
    if (filename.rfind(".matrix") != std::string::npos ||
         filename.rfind(".mat") != std::string::npos ||
         filename.rfind(MATRIXBINSUFFIX) != std::string::npos) {
        //!** yes, load \ref loadMatrixSingleBin(filename)
        return loadMatrixSingleBin(A, filename);
    }

    //!* no: check if filename is expandable with suffix ".matrix" or ".mat";
    if (fileExist(filename + ".matrix")) return loadMatrixSingleBin(A, filename + ".matrix");
    if (fileExist(filename + ".mat"))    return loadMatrixSingleBin(A, filename + ".mat");
    if (fileExist(filename + MATRIXBINSUFFIX))
        //!** yes , load \ref loadMatrixSingleBin(filename + \ref MATRIXBINSUFFIX)
        return loadMatrixSingleBin(A, filename + MATRIXBINSUFFIX);

    //!* no: try to load matrix from multiple binary vectors;
    return loadMatrixVectorsBin(A, filename);
}

DLLEXPORT bool loadMatrixSingleBin(RMatrix & A, const std::string & filename);
DLLEXPORT bool loadMatrixSingleBin(CMatrix & A, const std::string & filename);

/*! Force to load multiple binary vector files into one matrix (row-based). File name will be determined from filenamebody + successive increased number (read while files exist). \n
e.g. read "filename.0.* ... filename.n.* -> Matrix[0--n)[0..vector.size())\n
kCount can be given to use as subcounter. \n
e.g. read "filename.0_0.* ... filename.n_0.* ... filename.0_kCount-1.* ... filename.n_kCount-1.* ->
Matrix[0--n*kCount)[0..vector.size())
*/
DLLEXPORT bool loadMatrixVectorsBin(RMatrix & A, const std::string & filenameBody, uint kCount=1);
DLLEXPORT bool loadMatrixVectorsBin(CMatrix & A, const std::string & filenameBody, uint kCount=1);

/*! Save Matrix into Ascii File (column based). */
template < class ValueType >
bool saveMatrixCol(const Matrix < ValueType > & A, const std::string & filename){
    return saveMatrixCol(A, filename, "");
}

/*! Save Matrix into Ascii File (column based)  with optional comments header line. */
template < class ValueType >
bool saveMatrixCol(const Matrix < ValueType > & A, const std::string & filename,
                    const std::string & comments){
    std::fstream file; openOutFile(filename, & file, true);
    if (comments.length() > 0){
        file << "#" << comments << std::endl;
    }

    for (uint i = 0; i < A.cols(); i ++){
        for (uint j = 0; j < A.rows(); j ++){
            file << A[j][i] << "\t";
        }
        file << std::endl;
    }
    file.close();
    return true;
}

/*! Load Matrix from Ascii File (column based). */
template < class ValueType >
bool loadMatrixCol(Matrix < ValueType > & A, const std::string & filename){
    std::vector < std::string > comments;
    return loadMatrixCol(A, filename, comments);
}

/*! Load Matrix from Ascii File (column based), with optional comments header line. */
template < class ValueType >
bool loadMatrixCol(Matrix < ValueType > & A, const std::string & filename,
                    std::vector < std::string > & comments){

    uint commentCount = 0;
    uint cols = countColumnsInFile(filename, commentCount);
//     Index rows = countRowsInFile(filename);
//     // get length of file:
//     std::fstream file; openInFile(filename, & file, true);
//     Index length = fileLength(file);
//
//     // allocate memory:
//     char * buffer = new char[length];
//     file.read(buffer, length);
//     file.close();
//
//     delete buffer;
//     return true;

    Vector < ValueType > tmp;
    Vector < ValueType > row(cols);
    std::fstream file; openInFile(filename, & file, true);
    for (uint i = 0; i < commentCount; i ++) {
        std::string str;
        getline(file, str);
        comments = getSubstrings(str.substr(str.find('#'), -1));
    }

    double val;
    while(file >> val) tmp.push_back(val);

    file.close();
    Index rows = tmp.size() / cols ;
    A.resize(cols, rows);

    for (uint i = 0; i < rows; i ++){
        for (uint j = 0; j < cols; j ++){
            A[j][i] = tmp[i * cols + j];
        }
    }
    return true;
}

/*! Save Matrix into Ascii File (row based). */
template < class ValueType >
bool saveMatrixRow(const Matrix < ValueType > & A, const std::string & filename){
    return saveMatrixRow(A, filename, "");
}

/*! Save Matrix into Ascii File (row based)  with optional comments header line. */
template < class ValueType >
bool saveMatrixRow(const Matrix < ValueType > & A, const std::string & filename,
                    const std::string & comments){
    std::fstream file; openOutFile(filename, & file, true);
    if (comments.length() > 0){
        file << "#" << comments << std::endl;
    }

    for (uint i = 0; i < A.rows(); i ++){
        for (uint j = 0; j < A.cols(); j ++){
            file << A[i][j] << "\t";
        }
        file << std::endl;
    }
    file.close();
    return true;
}

/*! Load Matrix from Ascii File (row based). */
template < class ValueType >
bool loadMatrixRow(Matrix < ValueType > & A, const std::string & filename){

    std::vector < std::string > comments;
    return loadMatrixRow(A, filename, comments);
}

/*! Load Matrix from Ascii File (row based), with optional comments header line. */
template < class ValueType >
bool loadMatrixRow(Matrix < ValueType > & A,
                    const std::string & filename,
                    std::vector < std::string > & comments){

    uint commentCount = 0;
    uint cols = countColumnsInFile(filename, commentCount);

    Vector < ValueType > row(cols);
    std::fstream file; openInFile(filename, & file, true);
    for (uint i = 0; i < commentCount; i ++) {
        std::string str;
        getline(file, str);
        comments = getSubstrings(str.substr(str.find('#'), -1));
    }

    double val;
    Vector < ValueType > tmp;
    while(file >> val) tmp.push_back(val);

    file.close();
    Index rows = tmp.size() / cols ;
    A.resize(rows, cols);

    for (uint i = 0; i < rows; i ++){
        for (uint j = 0; j < cols; j ++){
            A[i][j] = tmp[i * cols + j];
        }
    }
    return true;
}

/*!Inplace matrix calculation: $C = a * A.T * B * A$ + b*C.
Size of A is (n,m) and B need to be square (n,n), C will resized to (m,m).
AtB might be for temporary memory allocation.  */
DLLEXPORT void matMultABA(const RDenseMatrix & A,
                          const RDenseMatrix & B,
                          RDenseMatrix & C,
                          RDenseMatrix & AtB,
                          const double & a=1.0, const double & b=0.0);
DLLEXPORT void matMultABA(const RMatrix & A,
                          const RMatrix & B,
                          RMatrix & C,
                          RMatrix & AtB,
                          const double & a=1.0, const double & b=0.0);

/*!Inplace matrix calculation: $C = a*A*B + b*C$. B are transposed if needed to fit appropriate dimensions. */
// DLLEXPORT void matMult(const RDenseMatrix & A,
//                        const RDenseMatrix & B,
//                        RDenseMatrix & C,
//                        const double & a=1.0, const double & b=0.0);
// DLLEXPORT void matMult(const RMatrix & A,
//                        const RMatrix & B,
//                        RMatrix & C,
//                        const double & a=1.0, const double & b=0.0);

/*!Inplace matrix calculation: $C = a * A.T * B + b*C$. B are transposed if needed to fit appropriate dimensions. */
// DLLEXPORT void matTransMult(const RDenseMatrix & A,
//                             const RDenseMatrix & B,
//                             RDenseMatrix & C,
//                             const double & a=1.0, const double & b=0.0);
// DLLEXPORT void matTransMult(const RMatrix & A,
//                             const RMatrix & B,
//                             RMatrix & C,
//                             const double & a=1.0, const double & b=0.0);

/*! Return determinant for Matrix(2 x 2). */
template < class T > inline T det(const T & a, const T & b,
                                  const T & c, const T & d){
    return a * d - b * c;
}

/*! Return determinant for Matrix A. This function is a stub. Only Matrix dimensions of 2 and 3 are considered. */
template < class ValueType > double det(const Matrix3< ValueType > & A){
    return A.det();
}

/*! Return determinant for Matrix A. This function is a stub. Only Matrix dimensions of 2 and 3 are considered. */
template < class Matrix > double det(const Matrix & A){
    //** das geht viel schoener, aber nicht mehr heute.;
    double det = 0.0;
    switch (A.rows()){
        case 2: det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
            break;
        case 3:
            det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                  A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                  A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

            // .T
            // det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])-
            //       A[1][0] * (A[0][1] * A[2][2] - A[2][1] * A[0][2])+
            //       A[2][0] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]);

            break;
        default:
            std::cerr << WHERE_AM_I << " matrix determinant of dim not yet implemented -- dim: " << A.rows() << std::endl;
            break;
    }
    return det;
}

/*! Return the inverse of Matrix A3. This function is a stub. Only Matrix dimensions of 2 and 3 are considered. */
template < class ValueType > inline Matrix3<ValueType> inv(const Matrix3< ValueType > & A){
    Matrix3< ValueType > I;
    inv(A, I);
    return I;
}

/*! Return the inverse of Matrix3 A. */
template < class ValueType > inline void inv(const Matrix3< ValueType > & A,
                                             Matrix3< ValueType > & I){
//     __M
//     std::cout << A << std::endl;
    I[0] =  (A[4] * A[8] - A[5] * A[7]);
    I[3] = -(A[3] * A[8] - A[5] * A[6]);
    I[6] =  (A[3] * A[7] - A[4] * A[6]);
    I[1] = -(A[1] * A[8] - A[2] * A[7]);
    I[4] =  (A[0] * A[8] - A[2] * A[6]);
    I[7] = -(A[0] * A[7] - A[1] * A[6]);
    I[2] =  (A[1] * A[5] - A[2] * A[4]);
    I[5] = -(A[0] * A[5] - A[2] * A[3]);
    I[8] =  (A[0] * A[4] - A[1] * A[3]);
//     std::cout << I << std::endl;
//     std::cout << (A[0] * I[0] + A[1] * I[3] + A[2] * I[6]) << std::endl;
//     std::cout << det(A) << std::endl;
    I /= (A[0] * I[0] + A[1] * I[3] + A[2] * I[6]);
}

/*! Return the inverse of Matrix A. This function is a stub. Only Matrix dimensions of 2 and 3 are considered. */
template < class Matrix > Matrix inv(const Matrix & A){
    //** das geht viel schoener, aber nicht mehr heute.; Wie?
    Matrix I(A.rows(), A.cols());
    inv(A, I);
    return I;
}

/*! Return the inverse of Matrix A. This function is a stub. Only Matrix dimensions of 2 and 3 are considered. */
template < class Matrix > void inv(const Matrix & A, Matrix & I){
    //** das geht viel schoener, aber nicht mehr heute.; Wie?

    switch (I.rows()){
        case 2:
            I[0][0] =  A[1][1];
            I[1][0] = -A[1][0];
            I[0][1] = -A[0][1];
            I[1][1] =  A[0][0];
            break;
        case 3:
            I[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]);
            I[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]);
            I[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
            I[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]);
            I[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]);
            I[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]);
            I[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]);
            I[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]);
            I[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
            break;
        default:
            std::cerr << WHERE_AM_I << " matrix determinant of dim not yet implemented -- dim: " << A.rows() << std::endl;
            break;
    }
    I /= det(A);
}

inline void save(const MatrixBase & A, const std::string & filename){
    A.save(filename);
}

inline void save(MatrixBase & A, const std::string & filename){
    A.save(filename);
}

inline RVector operator * (const MatrixBase & A, const RVector & b){
    return A.mult(b);
}
template < class ValueType, template < typename > class Mat >
Vector < ValueType > operator * (const Mat< ValueType > & A,
                                const Vector < ValueType > & b){
    return A.mult(b);
}

inline RVector operator * (const RMatrix & A, const RVector & b){
    return A.mult(b);
}
inline CVector operator * (const CMatrix & A, const CVector & b){
    return A.mult(b);
}

inline RVector transMult(const MatrixBase & A, const RVector & b){
    return A.transMult(b);
}
template < class ValueType, template < typename > class Mat >
Vector < ValueType > mult(const Mat< ValueType > & A,
                          const Vector < ValueType > & b){
    return A.mult(b);
}
template < class ValueType, template < typename > class Mat >
Vector < ValueType > transMult(const Mat< ValueType > & A,
                               const Vector < ValueType > & b){
    return A.transMult(b);
}

inline RMatrix real(const CMatrix & A){
    RMatrix R(A.rows(), A.cols());
    for (Index i = 0; i < A.rows(); i ++) R[i] = real(A[i]);
    return R;
}
inline RMatrix imag(const CMatrix & A){
    RMatrix R(A.rows(), A.cols());
    for (Index i = 0; i < A.rows(); i ++) R[i] = imag(A[i]);
    return R;
}

template < class T >
std::ostream & operator << (std::ostream & str, const Matrix < T > & M){
    for (Index i = 0; i < M.rows(); i ++) {
        str << M[i] << std::endl;
    }
    return str;
}
template < class T >
std::ostream & operator << (std::ostream & str, const DenseMatrix < T > & M){
    for (Index i = 0; i < M.rows(); i ++) {
        str << M[i] << std::endl;
    }
    return str;
}

} //namespace GIMLI

