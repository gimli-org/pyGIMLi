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

#include "gimli.h"
// #include "vector.h"
#include "elementmatrix.h"

#include "matrix.h"
#include "vector.h"
#include "stopwatch.h"

#if OPENBLAS_CBLAS_FOUND
    #include <cblas.h>
#endif


namespace GIMLI{

template class DenseMatrix< double >;

#define _MATRIX_WS_SIZE 8192

static double _wsA[_MATRIX_WS_SIZE];
static double _wsB[_MATRIX_WS_SIZE];
static double _wsC[_MATRIX_WS_SIZE];


static double __cblasTime__ = 0.0;
static double __cblasMinTime__ = std::numeric_limits<double>::max();
static Index __cblasCount__ = 0;

Index cblasCount(bool reset){
    if (reset) {
        Index r = __cblasCount__;
        __cblasCount__ = 0;
        return r;
    }
    return __cblasCount__;
}
void _updateCblasTime_(double t){
    __cblasTime__ += t;
    __cblasCount__ ++;
    __cblasMinTime__ = min(__cblasMinTime__, t);
}
double cblasSumTime(bool reset){
    if (reset) {
        double r = __cblasTime__;
        __cblasTime__ = 0.0;
        return r;
    }
    return __cblasTime__;
}

double cblasMinTime(bool reset){
    if (reset) {
        double r = __cblasMinTime__;
        __cblasMinTime__ = std::numeric_limits<double>::max();
        return r;
    }
    return __cblasMinTime__;
}


// void toEigenMatrix(const RMatrix & m, RMatrix & r){
//     r.resize(m.rows(), m.cols());

//     for (Index i=0; i < r.rows(); i ++){
// #if USE_EIGEN3
//         r(i, Eigen::all) = Eigen::Map <const Eigen::VectorXd>(&m[i][0],
//                                                               r.cols());
// #else
//     for (Index j=0; j < r.cols(); j ++){
//         r(i,j) = m(i,j);
//     }
// #endif
//     }
// }
// void toRMatrix(const RMatrix & m, RMatrix & r){
//     r.resize(m.rows(), m.cols());

//     // optimize if in use
//     for (Index i=0; i < r.rows(); i ++){
//         for (Index j=0; j < r.cols(); j ++){
//             r(i,j) = m(i,j);
//         }
//     }
// }

// #if USE_EIGEN3
// void toRVector(const Eigen::VectorXd & m, RVector & r, double b){
//     r.resize(m.size());

//     // optimize if in use
//     if (b == 0.0){
//         for (Index i=0; i < m.size(); i ++){
//             r[i] = m(i);
//         }
//     } else if (b == 1.0){
//         for (Index i=0; i < m.size(); i ++){
//             r[i] += m(i);
//         }
//     } else if (b == -1.0){
//         for (Index i=0; i < m.size(); i ++){
//             r[i] -= m(i);
//         }
//     } else {
//         for (Index i=0; i < m.size(); i ++){
//             r[i] = m(i) + r[i] * b;
//         }
//     }
// }
// #endif

//##############################################################################
// DenseMatrix related implementations
//##############################################################################
// template <> const Vector< double > DenseMatrix< double >::rowView(Index i) const {
//     ASSERT_THIS_SIZE(i)
//     return Vector< double >(this->_cols, _data, this->_cols * i);
// }
// template <> Vector< double > DenseMatrix< double >::rowView(Index i) {
//     ASSERT_THIS_SIZE(i)
//     return Vector< double >(this->_cols, _data, this->_cols * i);
// }
// template <> Vector< Complex > DenseMatrix< Complex >::rowView(Index i) const {
//     ASSERT_THIS_SIZE(i)
//     return Vector< Complex >(this->_cols, _data, this->_cols * i);
// }
// template <> Vector< Complex > DenseMatrix< Complex >::rowView(Index i) {
//     ASSERT_THIS_SIZE(i)
//     return Vector< Complex >(this->_cols, _data, this->_cols * i);
// }

// template <> const Vector< double > & DenseMatrix< double >::row(Index i) const {
//     ASSERT_THIS_SIZE(i)
//     _rowView->setBorrowedDataOffset(this->_cols * i);
//     // __MS(i, *_rowView)
//     return *_rowView;
// }
// template <> Vector< double > & DenseMatrix< double >::row(Index i) {
//     ASSERT_THIS_SIZE(i)
//     _rowView->setBorrowedDataOffset(this->_cols * i);
//     // __MS(i, *_rowView)
//     return *_rowView;
// }
// template <> const Vector< Complex > & DenseMatrix< Complex >::row(Index i) const {
//     ASSERT_THIS_SIZE(i)
//     _rowView->setBorrowedDataOffset(this->_cols * i);
//     return *_rowView;
// }
// template <> Vector< Complex > & DenseMatrix< Complex >::row(Index i) {
//     ASSERT_THIS_SIZE(i)
//     _rowView->setBorrowedDataOffset(this->_cols * i);
//     return *_rowView;
// }

template <> void
DenseMatrix< double >::round(const double & tol){
    Vector < double > view(length(), this->_data, 0);
    view.round(tol);
}
template <> void
DenseMatrix< Complex >::round(const Complex & tol){
    THROW_TO_IMPL
}

/*! Generic fall back implementation for c = alpha*(A*b) + beta*c. */
template < class ValueType, class Mat >
void mult_MV_fallback_T(const Mat & A,
                 const Vector < ValueType > & b, Vector < ValueType > & c,
                 const ValueType & alpha, const ValueType & beta,
                 Index bOff, Index cOff, bool trans=false){

    // __MS(A.cols(), A.rows(), typeid(A).name())
    if (trans){
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.rows())
        c.resize(A.cols() + cOff);
    } else {
        ASSERT_GREATER_EQUAL(b.size() + bOff, A.cols())
        c.resize(A.rows() + cOff);
    }

    Stopwatch sw;
    ValueType _c = 0.0;

    for (Index i = 0, iMax = c.size()-cOff; i < iMax; i ++){
        if (trans){
            _c = A.col(i).mult(b, bOff);
        } else {
            _c = A.row(i).mult(b, bOff);
        }

        if (alpha != 1.0){
            _c *= alpha;
        }

        if (beta == 0.0){
            c[i + cOff] = _c;
        } else if (beta == 1.0){
            c[i + cOff] += _c;
        } else if (beta == -1.0){
            c[i + cOff] -= _c;
        } else {
            c[i + cOff] = _c + beta * c[i + cOff];
        }
    }
    // print("mult impl:", sw.duration());
}

/*! Generic fall back implementation for C = a*(A*B) + b*C*/
template <class ValueType, class Mat >
void mult_MM_fallback_T(const Mat & A, const Mat & B, Mat & C,
                        const ValueType & alpha, const ValueType & beta,
                        bool bIsTrans, Index n){

    if (bIsTrans){
        ASSERT_EQUAL(A.cols(), B.cols())
        C.resize(A.rows(), B.rows());
    } else {
        ASSERT_EQUAL(A.cols(), B.rows())
        C.resize(A.rows(), B.cols());
    }

    Stopwatch sw;
    double c = 0;
    for (Index i = 0; i < A.rows(); i ++){
        const RVector & Ai = A[i];

        for (Index j = 0; j < n; j ++){
            c = 0;

            if (bIsTrans){
                const RVector & Bj = B[j];
                for (Index k = 0; k < A.cols(); k ++){
                    c += Ai[k] * Bj[k];
                }
            } else{
                for (Index k = 0; k < A.cols(); k ++){
                    c += Ai[k] * B[k][j];
                }
            }

            if (alpha != 1.0){
                c *= alpha;
            }
            if (beta == 0.0){
                C[i][j] = c;
            } else if (beta == 1.0){
                C[i][j] += c;
            } else if (beta == -1.0){
                C[i][j] -= c;
            } else {
                C[i][j] = beta * C[i][j] + c;
            }
        }
    }
    print("mult_MM_impl:", sw.duration());
}

template <class ValueType, class Mat >
void transMult_MM_fallback_T(const Mat & A, const Mat & B, Mat & C,
                             const ValueType & alpha, const ValueType & beta,
                             bool bIsTrans, Index n){

    if (bIsTrans){
        ASSERT_EQUAL(A.rows(), B.cols())
        C.resize(A.cols(), B.rows());
    } else {
        ASSERT_EQUAL(A.rows(), B.rows())
        C.resize(A.cols(), B.cols());
    }

    double c = 0;
    for (Index i = 0; i < A.cols(); i ++){
        for (Index j = 0; j < n; j ++){
            c = 0.0;
            for (Index k = 0; k < A.rows(); k ++){
                if (bIsTrans){
                    c += A[k][i] * B[j][k];
                } else {
                    c += A[k][i] * B[k][j];
                }
            }

            if (alpha != 1.0){
                c *= alpha;
            }

            if (beta == 0.0){
                C[i][j] = c;
            } else if (beta == 1.0){
                C[i][j] += c;
            } else if (beta == -1.0){
                C[i][j] -= c;
            } else {
                C[i][j] = beta * C[i][j] + c;
            }
        }
    }
}

template < class Mat >
void mult_MM_T(const Mat & A, const Mat & B,
               Mat & C, const double & a, const double & b){
    // C = a * A*B + b*C || C = a * A*B.T + b*C

    // __MS("matMult: ", A.rows(), A.cols(), B.rows(), B.cols())

    Index m = A.rows(); // C.rows()
    Index k = A.cols();

    // B.rows()
    Index n = B.cols(); // C.cols()
    bool bIsTrans = false;
    Index bRows = n;

    if (k == B.rows()){ // A * B (k == k)

    } else if (k == B.cols()){ // A * B.T (k == k)
        bIsTrans = true;
        n = B.rows();
        bRows = k;
    } else {
        log(Error, "matMult sizes mismatch. ",
            A.cols(), "!=", B.rows());
    }
    C.resize(m, n);

#if OPENBLAS_CBLAS_FOUND

    if (noCBlas()){
        return mult_MM_fallback_T(A, B, C, a, b, bIsTrans, n);
    }

    CBLAS_TRANSPOSE aTrans = CblasNoTrans;
    CBLAS_TRANSPOSE bTrans = CblasNoTrans;

    if (bIsTrans) bTrans = CblasTrans;

    // not threadsafe at all for using Buffer for RMatrix

    double * bA = A.toData(&_wsA[0], _MATRIX_WS_SIZE);
    double * bB = B.toData(&_wsB[0], _MATRIX_WS_SIZE);
    double * bC = C.toData(&_wsC[0], _MATRIX_WS_SIZE);

    // lda ## leading dimension for a, means column for CblasRowMajor
    Stopwatch sw;
    cblas_dgemm(CblasRowMajor, aTrans, bTrans,
                m, n, k,
                a, bA, k, bB, bRows,
                b, bC, n);
    _updateCblasTime_(sw.duration());
    C.fromData(bC, m, n);

    // vector is new so buffer was to small for RMatrix, which allocate her own
    if (bA != _wsA && A.rtti() == GIMLI_MATRIX_RTTI){
        delete [] bA;
        delete [] bB;
        delete [] bC;
    }
#else
    mult_MM_fallback_T(A, B, C, a, b, bIsTrans, n);
#endif
}

template < class Mat >
void mult_MM_T(const Mat & A, const Mat & B,
               Mat & C, const Complex & a, const Complex & b){
    THROW_TO_IMPL
}

template < class Mat >
void transMult_MM_T(const Mat & A, const Mat & B,
                    Mat & C, const double & a, const double & b){

    Index k = A.rows(); // B.rows()
    Index m = A.cols(); // C.rows()
    Index n = B.cols(); // C.cols()

    bool bIsTrans = false;
    Index bRows = n;

    if (k == B.rows()){ // A.T * B

        if (b == 0.0) C.resize(m, n);

        if (C.rows() != A.cols() || C.cols() != B.cols()){
            __MS("isinuse?")
            //__MS(C.rows(), C.cols(), A.cols(), B.cols())
            //** Target array have wrong dimensions
            if (C.rows() == B.cols() && C.cols() == A.cols()){
                // C = a * B.T*A + b*C
                //** Target array seems needed to be transposed
                //** C = a*(A.T*B).T + b * C
                // __MS("ret transmult")
                return transMult_MM_T(B, A, C, a, b);

                // retTrans = true;
            } else {
                //** resize target array
                C.resize(m, n);
            }
        }
    } else if(k == B.cols()){
        bIsTrans = true;
        n = B.rows();
        bRows = k;

    } else {
        // if not A.T * B or A.T * B.T
        __MS(A.rows(), A.cols(), '\n', A)
        __MS(B.rows(), B.cols(), '\n', B)
        __MS(C.rows(), C.cols())
        throwLengthError("matTransMult sizes mismatch.");
    }

#if OPENBLAS_CBLAS_FOUND

    if (noCBlas()){
        transMult_MM_fallback_T(A, B, C, a, b, bIsTrans, n);
    }

    CBLAS_TRANSPOSE aTrans = CblasTrans;
    CBLAS_TRANSPOSE bTrans = CblasNoTrans;

    if (bIsTrans) bTrans = CblasTrans;

    // not threadsafe at all for using Buffer for RMatrix

    double * bA = A.toData(&_wsA[0], _MATRIX_WS_SIZE);
    double * bB = B.toData(&_wsB[0], _MATRIX_WS_SIZE);
    double * bC = C.toData(&_wsC[0], _MATRIX_WS_SIZE);

    Stopwatch sw;
    cblas_dgemm(CblasRowMajor, aTrans, bTrans, m, n, k,
                a, bA, m, bB, bRows,
                b, bC, n);
    _updateCblasTime_(sw.duration());

    C.fromData(bC, m, n);

    if (bA != _wsA && A.rtti() == GIMLI_MATRIX_RTTI){
        delete [] bA;
        delete [] bB;
        delete [] bC;
    }

#else
    mult_MV_fallback_T(A, B, C, a, b, bIsTrans, n);
#endif

}

template < class Mat >
void transMult_MM_T(const Mat & A, const Mat & B,
                    Mat & C, const Complex & a, const Complex & b){
    THROW_TO_IMPL
}

/*! Special blas implementation of c = alpha*(A.T*b) + beta*c. */
void mult(const RDenseMatrix & A,
          const RVector & b, RVector & c,
          const double & alpha, const double & beta,
          Index bOff, Index cOff){

#if OPENBLAS_CBLAS_FOUND
    Stopwatch sw;
    if (noCBlas()){
        mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff);
    } else {

        Index m = A.rows();
        Index n = A.cols();
        ASSERT_VEC_SIZE(b, n)
        c.resize(m);

        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha,
                    A.pData(), n, &b[0],
                    1, beta, &c[0], 1);

        _updateCblasTime_(sw.duration());
    }
#else //#if OPENBLAS_CBLAS_FOUND
    mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff);
#endif
}
void mult(const CDenseMatrix & A,
          const CVector & b, CVector & c,
          const Complex & alpha, const Complex & beta,
          Index bOff, Index cOff){

#if OPENBLAS_CBLAS_FOUND
    Stopwatch sw;
    if (noCBlas()){
        mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff);
    } else {
        Index m = A.rows();
        Index n = A.cols();
        ASSERT_VEC_SIZE(b, n)
        c.resize(m);
        THROW_TO_IMPL
        // cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha,
        //             A.pData(), n, &b[0],
        //             1, beta, &c[0], 1);

        // print("dgemv:", sw.duration());
    }
#else //#if OPENBLAS_CBLAS_FOUND
    mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff);
#endif
}

/*! Special blas implementation of c = alpha*(A.T*b) + beta*c. */
void transMult(const RDenseMatrix & A,
               const RVector & b, RVector & c,
               const double & alpha, const double & beta,
               Index bOff, Index cOff){

#if OPENBLAS_CBLAS_FOUND
    Stopwatch sw;
    if (noCBlas()){
        mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff, true);
    } else {
        Index m = A.rows();
        Index n = A.cols();
        ASSERT_VEC_SIZE(b, m)
        c.resize(n);
        cblas_dgemv(CblasRowMajor, CblasTrans, m, n, alpha,
                    A.pData(), n, &b[0],
                    1, beta, &c[0], 1);

        _updateCblasTime_(sw.duration());
    }
#else //#if OPENBLAS_CBLAS_FOUND
    mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff, true);
#endif
}

/*! Special blas implementation of c = alpha*(A.T*b) + beta*c. */
void transMult(const CDenseMatrix & A,
               const CVector & b, CVector & c,
               const Complex & alpha, const Complex & beta,
               Index bOff, Index cOff){

#if OPENBLAS_CBLAS_FOUND
    Stopwatch sw;
    if (noCBlas()){
        mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff, true);
    } else {
        Index m = A.rows();
        Index n = A.cols();
        ASSERT_VEC_SIZE(b, m)
        c.resize(n);
        THROW_TO_IMPL
        // cblas_cgemv(CblasRowMajor, CblasNoTrans, m, n, alpha,
        //             A.pData(), n, &b[0],
        //             1, beta, &c[0], 1);

        _updateCblasTime_(sw.duration());
    }
#else //#if OPENBLAS_CBLAS_FOUND
    mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff, true);
#endif
}

void mult(const RDenseMatrix & A, const RDenseMatrix & B,
          RDenseMatrix & C, const double & a, const double & b){
    return mult_MM_T(A, B, C, a, b);
}
void mult(const CDenseMatrix & A, const CDenseMatrix & B,
          CDenseMatrix & C, const Complex & a, const Complex & b){
    return mult_MM_T(A, B, C, a, b);
}
void transMult(const RDenseMatrix & A, const RDenseMatrix & B,
               RDenseMatrix & C, const double & a, const double & b){
    transMult_MM_T(A, B, C, a, b);
}
void transMult(const CDenseMatrix & A, const CDenseMatrix & B,
               CDenseMatrix & C, const Complex & a, const Complex & b){
    transMult_MM_T(A, B, C, a, b);
}


//##############################################################################
// Matrix related implementations
//##############################################################################
void mult(const RMatrix & A,
          const RVector & b, RVector & c,
          const double & alpha, const double & beta,
          Index bOff, Index cOff){
    return mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff);
    // Bufferalloc for OPENBLAS to xpensive
}
void mult(const CMatrix & A,
          const CVector & b, CVector & c,
          const Complex & alpha, const Complex & beta,
          Index bOff, Index cOff){
    return mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff);
    // Bufferalloc for OPENBLAS to xpensive
}
void transMult(const RMatrix & A,
          const RVector & b, RVector & c,
          const double & alpha, const double & beta,
          Index bOff, Index cOff){
    return mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff, true);
    // Bufferalloc for OPENBLAS to xpensive
}
void transMult(const CMatrix & A,
          const CVector & b, CVector & c,
          const Complex & alpha, const Complex & beta,
          Index bOff, Index cOff){
    return mult_MV_fallback_T(A, b, c, alpha, beta, bOff, cOff, true);
    // Bufferalloc for OPENBLAS to xpensive
}

template < class ValueType > Matrix < ValueType > &
_transAdd(Matrix < ValueType > * a, const Matrix < ValueType > & b){
    // a+=b.T
    if (a->rows() != b.cols() || a->cols() != b.rows()){
        __MS(a->rows(), a->cols(), ":",  b.rows(), b.cols())
        // print(*a);
        // print(b);
        log(Error, "Matrix _transAdd with wrong dimensions");
        // #exit(-1);
        return *a;
    }

    for (Index i = 0; i < a->rows(); i ++ ){
        for (Index j = 0; j < a->cols(); j ++ ){
            a->mat_[i][j] += b.mat_[j][i];
        }
    }
    return *a;
}
template <> Matrix < double > &
Matrix<double>::transAdd(const Matrix < double > & a){
    return _transAdd(this, a);
}
template <> Matrix < Complex > &
Matrix<Complex>::transAdd(const Matrix < Complex > & a){
    return _transAdd(this, a);
}
template < class ValueType > DenseMatrix < ValueType > &
_transAdd(DenseMatrix < ValueType > * a, const DenseMatrix < ValueType > & b){
    // a+=b.T
    if (a->rows() != b.cols() || a->cols() != b.rows()){
        __MS(a->rows(), a->cols(), ":",  b.rows(), b.cols())
        // print(*a);
        // print(b);
        log(Error, "Matrix _transAdd with wrong dimensions");
        // #exit(-1);
        return *a;
    }

    for (Index i = 0; i < a->rows(); i ++ ){
        for (Index j = 0; j < a->cols(); j ++ ){
            a->_data[a->cols() * i + j] += b._data[b.cols() * j + i];
        }
    }
    return *a;
}
template <> DenseMatrix<double> &
DenseMatrix<double>::transAdd(const DenseMatrix < double > & a) {
    return _transAdd(this, a);
}
template <> DenseMatrix<Complex> &
DenseMatrix<Complex>::transAdd(const DenseMatrix < Complex > & a){
    return _transAdd(this, a);
}

template < class Mat >
void matMultABA_T(const Mat & A, const Mat & B,
                   Mat & C, Mat & AtB, const double & a, const double & b){
    // C = a A.T * B * A + b * C
    // __MS("matMultABA: ", A.rows(), A.cols(), B.rows(), B.cols())

    if (A.rows() != B.rows()){
        log(Error, "matMultABA B sizes mismatch.", A.rows(), "!=", B.rows());
        return;
    }
    AtB.resize(A.cols(), B.rows());
    transMult_MM_T(A, B, AtB, 1.0, 0.0);
    mult_MM_T(AtB, A, C, a, b);
}

void matMultABA(const RDenseMatrix & A, const RDenseMatrix & B,
                RDenseMatrix & C, RDenseMatrix & AtB,
                const double & a, const double & b){
    return matMultABA_T(A, B, C, AtB, a, b);
}
void matMultABA(const RMatrix & A, const RMatrix & B,
                RMatrix & C, RMatrix & AtB,
                const double & a, const double & b){

    return matMultABA_T(A, B, C, AtB, a, b);
}

void mult(const CMatrix & A, const CMatrix & B, CMatrix & C,
          const Complex & alpha, const Complex & beta){
    THROW_TO_IMPL
}

void transMult(const CMatrix & A, const CMatrix & B, CMatrix & C,
               const Complex & alpha, const Complex & beta){
    THROW_TO_IMPL
}

void mult(const RMatrix & A, const RMatrix & B,
          RMatrix & C, const double & a, const double & b){
    return mult_MM_T(A, B, C, a, b);
}
void transMult(const RMatrix & A, const RMatrix & B,
               RMatrix & C, const double & a, const double & b){
    return transMult_MM_T(A, B, C, a, b);
}


// void matTransMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a, double b){
//     //** C = a * A.T*B + b*C|| C = a * A.T*B.T  + b*C
//     //** C = (a * A.T*B).T + b*C if C has the right dimension
//     //** implement with openblas dgemm too and check performance
//     // __MS("matTransMult: "<< A.rows() << " " << A.cols() << " : "
//     //     << B.rows() << " " << B.cols() << " " << a << " " << b)
//     // bool retTrans = false;

//     // A(k, m).T * B(k, n) = C(m, n)

//     Index k = A.rows(); // B.rows()
//     Index m = A.cols(); // C.rows()
//     Index n = B.cols(); // C.cols()

//     if (k == B.rows()){ // A.T * B

//         if (C.rows() != A.cols() || C.cols() != B.cols()){

//             // __MS(C.rows() << " " << C.cols() << " " << A.cols() << " " << B.cols())
//             //** Target array have wrong dimensions
//             if (C.rows() == B.cols() && C.cols() == A.cols()){
//                 // C = a * B.T*A + b*C
//                 //** Target array seems needed to be transposed
//                 //** C = a*(A.T*B).T + b * C
//                 // __MS("ret transmult")
//                 return matTransMult(B, A, C, a, b);

//                 //retTrans = true;
//             } else {
//                 //** resize target array
//                 C.resize(m, n);
//             }
//         }

// #if OPENBLAS_CBLAS_FOUND

//         double *A2 = new double[k * m];
//         double *B2 = new double[k * n];
//         double *C2 = new double[m * n];

//         A.dumpData(A2);
//         B.dumpData(B2);
//         C.dumpData(C2);

//     // std::cout << "A" << std::endl;
//     // for (Index i = 0; i < m*k; i ++ ){std::cout << A2[i] << " ";} std::cout << std::endl;
//     // std::cout << "b" << std::endl;
//     // for (Index i = 0; i < n*k; i ++ ){std::cout << B2[i] << " ";} std::cout << std::endl;
//     // std::cout << "C1" << std::endl;
//     // for (Index i = 0; i < m*n; i ++ ){std::cout << C2[i] << " ";} std::cout << std::endl;

//         cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k,
//                     a, A2, m, B2, n, b, C2, n);

//     // std::cout << "C2" << std::endl;
//     // for (Index i = 0; i < m*n; i ++ ){std::cout << C2[i] << " ";} std::cout << std::endl;

//         C.fromData(C2, m, n);

//     // std::cout << "C3" << std::endl;
//     // std::cout << C << std::endl;

//         delete [] A2;
//         delete [] B2;
//         delete [] C2;

// #else

//         for (Index i = 0; i < A.cols(); i ++){
//             for (Index j = 0; j < B.cols(); j ++){
//                 double c = 0;
//                 for (Index k = 0; k < A.rows(); k ++){
//                     c += A[k][i] * B[k][j];
//                 }
//                 // if (retTrans){
//                 if (0){
//                     THROW_TO_IMPL
//                     C[j][i] = a * c + C[j][i]*b;
//                 } else {
//                     C[i][j] = a * c + C[i][j]*b;
//                 }
//             }
//         }
// #endif
//     } else { // A.T * B.T
//         __MS(A)
//         __MS(B)
//         log(Error, "matTransMult sizes mismatch.", A.rows(), "!=", B.rows());
//     }
// }

/*! Force to load single matrix binary file.
    Format: see \ref save(const Matrix < ValueType > & A, const std::string & filename). */
template < class ValueType >
bool loadMatrixSingleBin_T(Matrix < ValueType > & A,
                          const std::string & filename){

    std::ifstream testFile(filename, std::ios::binary);
    const auto begin = testFile.tellg();
    testFile.seekg (0, std::ios::end);
    const auto end = testFile.tellg();
    SIndex fsize = (end-begin);
    testFile.close();
    FILE *file; file = fopen(filename.c_str(), "r+b");
    if (!file) {
        throwError(WHERE_AM_I + " " +
                   filename + ": " + strerror(errno));
    }

    Index ret;
    uint32 rows = 0;
    ret = fread(&rows, sizeof(uint32), 1, file);
    if (ret == 0) throwError("fail reading file " + filename);
    uint32 cols = 0;
    ret = fread(&cols, sizeof(uint32), 1, file);
    if (ret == 0) throwError("fail reading file " + filename);

    if ((SIndex)(rows*cols*sizeof(ValueType) + 2*sizeof(uint32)) != fsize){
        __MS("rows: ", str(rows), " cols: ", cols, " fsize: ",fsize)
        __MS(" filesize needed: ", rows*cols*sizeof(ValueType)+2*sizeof(uint32))
        fclose(file);
        throwError(WHERE_AM_I + " " + filename + ": size invalid");
    }

    A.resize(rows, cols);
    for (uint32 i = 0; i < rows; i ++){
        for (uint32 j = 0; j < cols; j ++){
            ret = fread((char*)&A[i][j], sizeof(ValueType), 1, file);
            if (ret == 0) throwError("fail reading file " + filename);
        }
    }
    fclose(file);
    A.rowFlag().fill(1);
    return true;
}

bool loadMatrixSingleBin(RMatrix & A, const std::string & filename){
    return loadMatrixSingleBin_T(A, filename);
}
bool loadMatrixSingleBin(CMatrix & A, const std::string & filename){
    return loadMatrixSingleBin_T(A, filename);
}

template < class ValueType >
bool loadMatrixVectorsBin_T(Matrix < ValueType > & A,
                            const std::string & filenameBody, uint kCount = 1){

    A.clear();
    Vector < ValueType > tmp;
    std::string filename;

    for (Index i = 0; i < kCount; i++){
        Index count = 0;
        while (1){ // load as long as posible
            if (kCount > 1){
                filename = filenameBody + "." + str(count) + "_" + str(i) + ".pot";
            } else {
                filename = filenameBody + "." + str(count) + ".pot";
            }

            if (!fileExist(filename)){
                filename = filenameBody + "." + str(count);
                if (!fileExist(filename)){
                    if (count == 0) {
	                    std::cerr << "Can't found: " << filename << std::endl;
                    }
                    break;
                }
            }
            if (GIMLI::load(tmp, filename, Binary)) A.push_back(tmp);
            count ++;
        } // while files exist
    } // for each k count
    return true;
}

bool loadMatrixVectorsBin(RMatrix & A,
                           const std::string & filenameBody, uint kCount){
    return loadMatrixVectorsBin_T(A, filenameBody, kCount);
}
bool loadMatrixVectorsBin(CMatrix & A,
                           const std::string & filenameBody, uint kCount){
    return loadMatrixVectorsBin_T(A, filenameBody, kCount);
}


} // namespace GIMLI{
