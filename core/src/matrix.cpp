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

template < class ValueType > Vector < ValueType >
_mult(const Matrix< ValueType > & M, const Vector < ValueType > & b) {
    Index cols = M.cols();
    Index rows = M.rows();

    Vector < ValueType > ret(rows, 0.0);

    //ValueType tmpval = 0;
    if (b.size() == cols){
        for (Index i = 0; i < rows; ++i){
            ret[i] = sum(M.mat_[i] * b);
            // for (Index j = 0; j < cols; j++){
            //     ret[i] += M.mat_[i][j] * b[j];
            // }
        }
    } else {
        throwLengthError(WHERE_AM_I + " " + str(cols) + " != " + str(b.size()));
    }
    return ret;
}

template<> Vector < double >
Matrix< double >::mult(const Vector < double > & b) const { return _mult((*this), b); }
template<> Vector < Complex >
Matrix< Complex >::mult(const Vector < Complex > & b) const { return _mult((*this), b); }

template < class ValueType > Vector < ValueType >
_mult(const Matrix< ValueType > & M, const Vector < ValueType > & b, Index startI, Index endI) {
    Index cols = M.cols();
    Index rows = M.rows();
    Index bsize = Index(endI - startI);

    if (bsize != cols) {
        throwLengthError(WHERE_AM_I + " " + str(cols) + " < " + str(endI) + "-" + str(startI));
    }
    Vector < ValueType > ret(rows, 0.0);
    for (Index i = 0; i < rows; ++i){
        for (Index j = startI; j < endI; j++) {
            ret[i] += M.mat_[i][j] * b[j];
        }
    }
    return ret;
}

template<> Vector < double >
Matrix< double >::mult(const Vector < double > & b, Index startI, Index endI) const {
    return _mult((*this), b, startI, endI);
}
template<> Vector < Complex >
Matrix< Complex >::mult(const Vector < Complex > & b, Index startI, Index endI) const {
    return _mult((*this), b, startI, endI);
}

template < class ValueType > Vector < ValueType >
_transMult(const Matrix< ValueType > & M, const Vector < ValueType > & b) {
    Index cols = M.cols();
    Index rows = M.rows();
    Vector < ValueType > ret(cols, 0.0);

    if (b.size() == rows){
        for (Index i = 0; i < rows; i++){
            // ret += M.mat_[i] * b[i];
            for (Index j = 0; j < cols; j++){
                ret[j] += M.mat_[i][j] * b[i];
            }
        }
    } else {
        throwLengthError(WHERE_AM_I + " " + str(rows) + " != " + str(b.size()));
    }
    return ret;
}

template<> Vector< double >
Matrix< double >::transMult(const Vector < double > & b) const {
    return _transMult((*this), b);
}
template<> Vector< Complex >
Matrix< Complex >::transMult(const Vector < Complex > & b) const {
    return _transMult((*this), b);
}

template < class ValueType > Matrix < ValueType > &
_transAdd(Matrix < ValueType > * a, const Matrix < ValueType > & b){
    if (a->rows() != b.cols() || a->cols() != b.rows()){
        __MS(a->rows() << " " << b.cols() << " " << a->cols() << " " << b.rows())
        log(Error, "Matrix trans add with wrong dimensions");
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


void matMultABA(const RMatrix & A, const RMatrix & B, RMatrix & C,
                RMatrix & AtB, double a, double b){
    // C = a A.T * B * A + b * C
    // __MS("matMultABA: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())

    if (A.rows() != B.rows()){
        log(Error, "matMultABA B sizes mismatch.", A.rows(), "!=", B.rows());
        return;
    }
    AtB.resize(A.cols(), B.rows());
    matTransMult(A, B, AtB, 1.0, 0.0);
    matMult(AtB, A, C, a, b);
}

void matMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a, double b){
    // C = a * A*B + b *C || C += a * A*B.T + b*C
    // __MS("matMult: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())
    Index m = A.rows(); // C.rows()
    Index n = B.cols(); // C.cols()
    Index k = A.cols(); // B.rows()

    if (k == B.rows()){ // A * B (k == k)
        C.resize(m, n);

#if OPENBLAS_CBLAS_FOUND

        double *A2 = new double[m * k];
        double *B2 = new double[k * n];
        double *C2 = new double[m * n];

        A.dumpData(A2);
        B.dumpData(B2);
        C.dumpData(C2);
        // lda ## leading dimension for a, means column for CblasRowMajor
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                    a, A2, k, B2, n,
                    b, C2, n);

        C.fromData(C2, m, n);

        delete [] A2;
        delete [] B2;
        delete [] C2;
#else
    // __MS("matMult: "<< A.rows() << " " << A.cols() << " : " << B.rows() << " " << B.cols())

        for (Index i = 0; i < A.rows(); i ++){
            for (Index j = 0; j < B.cols(); j ++){
                double c = 0;
                for (Index k = 0; k < A.cols(); k ++){
                    c += A[i][k] * B[k][j];
                }
                C[i][j] = b * C[i][j] + a * c;
            }
        }
#endif
    } else { // A * B.T
        log(Error, "matMult sizes mismatch. implement fallback A*.B.T", A.cols(), "!=", B.rows());
    }
}

void matTransMult(const RMatrix & A, const RMatrix & B, RMatrix & C, double a, double b){
    //** C = a * A.T*B + b*C|| C = a * A.T*B.T  + b*C
    //** C = (a * A.T*B).T + b*C if C has the right dimension
    //** implement with openblas dgemm too and check performance
    // __MS("matTransMult: "<< A.rows() << " " << A.cols() << " : "
    //     << B.rows() << " " << B.cols() << " " << a << " " << b)
    // bool retTrans = false;

    // A(k, m).T * B(k, n) = C(m, n)

    Index k = A.rows(); // B.rows()
    Index m = A.cols(); // C.rows()
    Index n = B.cols(); // C.cols()

    if (k == B.rows()){ // A.T * B

        if (C.rows() != A.cols() || C.cols() != B.cols()){

            // __MS(C.rows() << " " << C.cols() << " " << A.cols() << " " << B.cols())
            //** Target array have wrong dimensions
            if (C.rows() == B.cols() && C.cols() == A.cols()){
                // C = a * B.T*A + b*C
                //** Target array seems needed to be transposed
                //** C = a*(A.T*B).T + b * C
                // __MS("ret transmult")
                return matTransMult(B, A, C, a, b);

                //retTrans = true;
            } else {
                //** resize target array
                C.resize(m, n);
            }
        }

#if OPENBLAS_CBLAS_FOUND

        double *A2 = new double[k * m];
        double *B2 = new double[k * n];
        double *C2 = new double[m * n];

        A.dumpData(A2);
        B.dumpData(B2);
        C.dumpData(C2);

    // std::cout << "A" << std::endl;
    // for (Index i = 0; i < m*k; i ++ ){std::cout << A2[i] << " ";} std::cout << std::endl;
    // std::cout << "b" << std::endl;
    // for (Index i = 0; i < n*k; i ++ ){std::cout << B2[i] << " ";} std::cout << std::endl;
    // std::cout << "C1" << std::endl;
    // for (Index i = 0; i < m*n; i ++ ){std::cout << C2[i] << " ";} std::cout << std::endl;

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k,
                    a, A2, m, B2, n, b, C2, n);

    // std::cout << "C2" << std::endl;
    // for (Index i = 0; i < m*n; i ++ ){std::cout << C2[i] << " ";} std::cout << std::endl;

        C.fromData(C2, m, n);

    // std::cout << "C3" << std::endl;
    // std::cout << C << std::endl;

        delete [] A2;
        delete [] B2;
        delete [] C2;

#else

        for (Index i = 0; i < A.cols(); i ++){
            for (Index j = 0; j < B.cols(); j ++){
                double c = 0;
                for (Index k = 0; k < A.rows(); k ++){
                    c += A[k][i] * B[k][j];
                }
                // if (retTrans){
                if (0){
                    THROW_TO_IMPL
                    C[j][i] = a * c + C[j][i]*b;
                } else {
                    C[i][j] = a * c + C[i][j]*b;
                }
            }
        }
#endif
    } else { // A.T * B.T
        __MS(A)
        __MS(B)
        log(Error, "matTransMult sizes mismatch.", A.rows(), "!=", B.rows());
    }
}

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
        __MS("rows: " << str(rows) << " cols: " <<  str(cols) << " fsize: " << str(fsize))
        __MS(" filesize needed: " << str(rows*cols*sizeof(ValueType)+2*sizeof(uint32)))
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
            if (loadVec(tmp, filename, Binary )) A.push_back(tmp);
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
