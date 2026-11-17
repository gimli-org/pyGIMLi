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

#ifndef GIMLI_SPARSEMATRIX__H
#define GIMLI_SPARSEMATRIX__H

#include "gimli.h"
#include "elementmatrix.h"
#include "vector.h"
#include "vectortemplates.h"
#include "matrix.h"
#include "mesh.h"
#include "meshentities.h"
#include "node.h"
#include "stopwatch.h"
#include "sparsemapmatrix.h"

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>

namespace GIMLI{


#define SPARSE_NOT_VALID throwError(WHERE_AM_I + " no data/or sparsity pattern defined.");

//! Sparse matrix in compressed row storage (CRS) form
/*! Sparse matrix in compressed row storage (CRS) form.
* IF you need native CCS format you need to transpose CRS
* Symmetry type: 0 = nonsymmetric, -1 symmetric lower part, 1 symmetric upper part.*/

template < class ValueType > class SparseMatrix : public SparseMatrixBase{
public:

  /*! Default constructor. Builds invalid sparse matrix */
    SparseMatrix()
        : SparseMatrixBase(), valid_(false), stype_(0){ }

    /*! Copy constructor. */
    SparseMatrix(const SparseMatrix < ValueType > & S)
        : SparseMatrixBase(),
          colPtr_(S.vecColPtr()),
          rowIdx_(S.vecRowIdx()),
          vals_(S.vecVals()), valid_(true), stype_(0){
          _rows = S.rows();
          _cols = S.cols();
    }

    /*! Copy constructor. */
    #ifndef PYGIMLI_CAST // disallow automatic type conversion in python
    SparseMatrix(const SparseMapMatrix< ValueType, Index > & S)
        : SparseMatrixBase(), valid_(true){
        this->copy_(S);
    }
    #endif
    SparseMatrix(const IndexArray & colPtr,
                 const IndexArray & rowIdx,
                 const Vector < ValueType > vals, int stype=0)
        : SparseMatrixBase(){
        WITH_TICTOC("SparseMatrix(c, r, v).1")
        colPtr_ = std::vector < int >(colPtr.size());
        // __MS(colPtr)
        // __MS(rowIdx)
        // __MS(vals)
        rowIdx_ = std::vector < int >(rowIdx.size());
        for (Index i = 0; i < colPtr_.size(); i ++ ) colPtr_[i] = colPtr[i];
        for (Index i = 0; i < rowIdx_.size(); i ++ ) rowIdx_[i] = rowIdx[i];
        vals_   = vals;
        stype_  = stype;
        valid_  = true;
        _cols = max(rowIdx_) + 1;
        _rows = colPtr_.size() - 1;
    }

    SparseMatrix(const std::vector < int > & colPtr,
                 const std::vector < int > & rowIdx,
                 const Vector < ValueType > vals, int stype=0)
        : SparseMatrixBase(){
        WITH_TICTOC("SparseMatrix(c, r, v).2")
          colPtr_ = colPtr;
          rowIdx_ = rowIdx;
          vals_   = vals;
          stype_  = stype;
          valid_  = true;
          _cols = max(rowIdx_) + 1;
          _rows = colPtr_.size() - 1;
    }

    /*! Destructor */
    virtual ~SparseMatrix(){ }

    /*! Copy assignment operator. */
    SparseMatrix < ValueType > & operator = (const SparseMatrix < ValueType > & S){
        WITH_TICTOC("SparseMatrix=().1")
        if (this != &S){
            this->copy(S);
        } return *this;
    }

    void copy_( const SparseMatrix< ValueType > & S){
        colPtr_ = S.vecColPtr();
        rowIdx_ = S.vecRowIdx();
        vals_   = S.vecVals();
        stype_  = S.stype();
        valid_  = true;
        _cols = S.cols();
        _rows = S.rows();
    }
    /*! Inline copy.  */
    void copy(const SparseMatrix< double > & S){ this->copy_(S); }
    /*! Inline copy. */
    void copy(const SparseMatrix< Complex > & S){ this->copy_(S); }

    void copy_(const SparseMapMatrix< double, Index > & S);

    void copy_(const SparseMapMatrix< Complex, Index > & S);

    /*! Inline copy.  */
    void copy(const SparseMapMatrix< double, Index > & S){ this->copy_(S); }
    /*! Inline copy. */
    void copy(const SparseMapMatrix< Complex, Index > & S){ this->copy_(S); }

    SparseMatrix < ValueType > & operator = (const SparseMapMatrix< ValueType, Index > & S){
        WITH_TICTOC("SparseMatrix=().2")
        this->copy_(S);
        return *this;
    }

    #define DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(OP, FUNCT) \
        void FUNCT(int i, int j, ValueType val){ \
            ASSERT_NON_EMPTY(colPtr_) \
            if ((stype_ < 0 && i > j) || (stype_ > 0 && i < j)) return; \
            if (abs(val) > TOLERANCE || 1){ \
                for (int k = colPtr_[i]; k < colPtr_[i + 1]; k ++){ \
                    if (rowIdx_[k] == j) { \
                        vals_[k] OP##= val; return; \
                    } \
                } \
                std::cerr << WHERE_AM_I << " pos " << i << " " << j << " is not part of the sparsity pattern " << std::endl; \
            } \
        }\
        SparseMatrix< ValueType > & operator OP##= (const ValueType & v){ \
            vals_ OP##= v; \
            return *this;\
        }\

        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(+, addVal)
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(-, subVal)
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(*, mulVal)
        DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__(/, divVal)

    #undef DEFINE_SPARSEMATRIX_UNARY_MOD_OPERATOR__

    SparseMatrix< ValueType > & operator += (const SparseMatrix< ValueType > & A){
        // pattern check to expensive
        ASSERT_EQUAL_SIZE(this->values(), A.values())
        vals_ += A.values();
        return *this;
    }
    SparseMatrix< ValueType > & operator -= (const SparseMatrix< ValueType > & A){
        ASSERT_EQUAL_SIZE(this->values(), A.values())
        vals_ -= A.values();
        return *this;
    }

    SparseMatrix< ValueType > & operator += (const ElementMatrix< double > & A){
        if (!valid_) SPARSE_NOT_VALID;
        Index As = A.size();
        Index Ai = 0;
        for (Index i = 0; i < As; i++){
            Ai = A.idx(i);
            for (Index j = 0; j < As; j++){
                addVal(Ai, A.idx(j), A.getVal(i, j));
            }
        }
        return *this;
    }

    virtual uint rtti() const { return GIMLI_SPARSE_CRS_MATRIX_RTTI; }


    // void add(const ElementMatrix< double > & A, const double & scale){
    //     return add(A, ValueType(scale), 1.0scale);
    // }
    virtual void add(const ElementMatrix< double > & A,
                     const ValueType & f, const double & scale=1.0);
    virtual void add(const ElementMatrix< double > & A,
                     const Pos & f, const double & scale=1.0);
    virtual void add(const ElementMatrix< double > & A,
                     const SmallMatrix < ValueType > & f, const double & scale=1.0);

    /*! Perftest .. maybe optimizer problem. Single calls outside calls of addVall suffer from polymorphism. */
    void addS(const ElementMatrix< double > & A, const ValueType & f,
              const double & scale){

        ValueType b = f*scale;

        A.integrate();

        for (Index i = 0, imax = A.rows(); i < imax; i++){
            for (Index j = 0, jmax = A.cols(); j < jmax; j++){
                this->addVal(A.rowIDs()[i], A.colIDs()[j], b * A.getVal(i, j));
            }
        }
    }

    /*! Set all values to zero but keep sparsity pattern.*/
    void clean(){ this->vals_*=ValueType(0); }

    /*! Clear the whole matrix and remove values and sparsity pattern.*/
    void clear(){
        colPtr_.clear();
        rowIdx_.clear();
        vals_.clear();
        valid_ = false;
        _cols = 0;
        _rows = 0;
    }

    virtual void setVal(Index i, Index j, const ValueType & val);

    /*!Get matrix value at i,j. If i and j is not part of the matrix
     * sparsity pattern return 0 and print a warning.
     * This warning can be disabled by setting warn to false.*/
    virtual ValueType getVal(Index i, Index j, bool warn=true) const;

    void cleanRow(int row){
        ASSERT_RANGE(row, 0, (int)this->rows())
        for (int col = colPtr_[row]; col < colPtr_[row + 1]; col ++){
            vals_[col] = ValueType(0);
        }
    }

    void cleanCol(int col){
        ASSERT_RANGE(col, 0, (int)this->cols())
        for (int i = 0; i < (int)this->rowIdx_.size(); i++){
            if (rowIdx_[i] == col) {
                vals_[i] = ValueType(0);
            }
        }
    }

    /**
     * @brief Sets the values in the sparse matrix based on the provided mask and values.
     *
     * This function updates the values in the sparse matrix using the indices specified
     * in the mask and the corresponding values provided in the vals vector. The size of
     * the mask and vals vectors must be the same.
     *
     * @param mask A vector of indices where the values need to be updated.
     * @param vals A vector of values to be set at the specified indices in the mask.
     *
     * @throws std::runtime_error if the sizes of mask and vals do not match.
     */
    void setMaskValues(const IndexArray & mask,
                       const Vector < ValueType > & vals){
        ASSERT_EQUAL(mask.size(), vals.size())
        Index i = 0;
        for (auto id: mask){
            vals_[id] = vals[i];
            i++;
        }
    }
    /**
     * @brief Sets the values in the sparse matrix to a specified value based on a mask.
     *
     * This function iterates over the indices provided in the mask and sets the corresponding
     * values in the sparse matrix to the specified value.
     *
     * @param mask An array of indices indicating which elements in the sparse matrix should be set.
     * @param val The value to set for all elements specified by mask.
     */
    void setMaskValues(const IndexArray & mask, const ValueType & val){
        for (auto id: mask){
            vals_[id] = val;
        }
    }

    /**
     * @brief Reduces the sparse matrix by removing rows and columns specified by the given indices.
     *
     * @param ids A vector of indices specifying which rows and columns to remove.
     * @param keepDiag A boolean flag indicating whether to keep the diagonal elements.
     *                 If true, diagonal elements are preserved even if their indices are in the ids vector.
     */
    void reduce(const IVector & ids, bool keepDiag=true);

    /*! Create value mask to reduce rows and columns.*/
    /**
     * @brief Creates a reduction mask based on the provided indices.
     *
     * This function generates an IndexArray that serves as a mask for reducing
     * the elements of a given vector or matrix. The mask is created based on
     * the indices provided in the input vector `ids`.
     *
     * @param ids A vector of indices used to create the reduction mask.
     * @return An IndexArray representing the reduction mask.
     */
    IndexArray createReduceMask(const IVector & ids);

    /**
     * @brief Creates a mask for the diagonal elements of the sparse matrix.
     *
     * This function generates an index array that represents the positions of the
     * diagonal elements in the sparse matrix. It can be used to efficiently access
     * or modify the diagonal elements.
     *
     * @return IndexArray An array containing the indices of the diagonal elements.
     */
    IndexArray createDiagonalMask();

    void buildSparsityPattern(const Mesh & mesh);
    void buildSparsityPattern(const std::vector < std::set< Index > > & idxMap);
    void addSparsityPattern(const std::vector < std::set< Index > > & idxMap);

    void fillStiffnessMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern=true);
    void fillStiffnessMatrix(const Mesh & mesh){
        RVector a(mesh.cellCount(), 1.0);
        fillStiffnessMatrix(mesh, a);
    }

    void fillMassMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern=true);
    void fillMassMatrix(const Mesh & mesh){
        RVector a(mesh.cellCount(), 1.0);
        fillMassMatrix(mesh, a);
    }

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

    /*! symmetric type. 0 = nonsymmetric, -1 symmetric lower part, 1 symmetric upper part.*/
    inline int stype() const {return stype_;}

    inline int * colPtr() { if (valid_) return &colPtr_[0]; else
                                        SPARSE_NOT_VALID;  return 0; }
    inline const int & colPtr() const { if (valid_) return colPtr_[0]; else
                                       SPARSE_NOT_VALID; return colPtr_[0]; }
    inline std::vector < int > & vecColPtr() { return colPtr_; }
    inline const std::vector < int > & vecColPtr() const { return colPtr_; }

    inline int * rowIdx() { if (valid_) return &rowIdx_[0]; else
                            SPARSE_NOT_VALID; return 0; }
    inline const int & rowIdx() const { if (valid_) return rowIdx_[0]; else
                                        SPARSE_NOT_VALID; return rowIdx_[0]; }
    inline std::vector < int > & vecRowIdx() { return rowIdx_; }
    inline const std::vector < int > & vecRowIdx() const { return rowIdx_; }

    inline void fillIndices(IndexArray & ids) const {
        ids.resize(rowIdx_.size());
        for (Index i = 0; i < ids.size(); i ++) {ids[i] = Index(rowIdx_[i]);}
    }
    inline void fillIndptr(IndexArray & ptr) const {
        ptr.resize(colPtr_.size());
        for (Index i = 0; i < ptr.size(); i ++) {ptr[i] = Index(colPtr_[i]);}
    }

    //!! check if needed
    inline ValueType * pVals() { if (valid_) return &vals_[0]; else SPARSE_NOT_VALID; return 0; }
//     inline const ValueType * vals() const { if (valid_) return &vals_[0]; else SPARSE_NOT_VALID; return 0; }
//     inline const ValueType & vals() const { if (valid_) return vals_[0]; else SPARSE_NOT_VALID; return vals_[0]; }
    inline const Vector < ValueType > & vecVals() const { return vals_; }
    inline Vector < ValueType > & vecVals() { return vals_; }

    //!!preferred
    inline const Vector < ValueType > & values() const { return vals_; }
    inline Vector < ValueType > values(const IndexArray & mask)
        const { return vals_[mask]; }
    inline Vector < ValueType > & values() { return vals_; }

    void update(const Vector < ValueType > & v) { vals_ = v; }

    // inline Index size() const { return rows(); }
    inline Index nVals() const { return vals_.size(); }
    // inline Index cols() const { return _cols; }
    // inline Index rows() const { return _rows; }
    inline Index nCols() const { return _cols; }
    inline Index nRows() const { return _rows; }

    /*! Return copy of the c-th column. */
    Vector < ValueType > col(Index c) const{
        ASSERT_RANGE(c, 0, cols())
        Vector < ValueType > b(cols(), 0.0);
        b[c] = 1.0;
        return this->mult(b);
    }
    /*! Return copy of the r-th row. */
    Vector < ValueType > row(Index r) const {
        ASSERT_RANGE(r, 0, rows())
        Vector < ValueType > ret(cols(), 0.0);

        for (int col = colPtr_[r]; col < colPtr_[r + 1]; col ++){
            ret[rowIdx_[col]] = vals_[col];
        }
        return ret;
    }

    void save(const std::string & fileName) const {
        if (!valid_) SPARSE_NOT_VALID;
        std::fstream file;
        openOutFile(fileName, & file);

        file.setf(std::ios::scientific, std::ios::floatfield);
        file.precision(14);

        for (Index i = 0; i < size(); i++){
            for (SIndex j = colPtr_[i]; j < colPtr_[i + 1]; j ++){
                file << i << "\t" << rowIdx_[j]
                          << "\t" << vals_[j] << std::endl;
            }
        }
        file.close();
    }

    bool valid() const { return valid_; }
    void setValid(bool v) { valid_ = v; }

    /*! Just set the matrix dimensions, not useful until proper fill of pattern. Only use for internal use.*/
    inline virtual void resize(Index rows, Index cols){
            // __MS(_rows, colPtr_.size() - 1, rows)
            // __MS(_cols, max(rowIdx_) + 1, cols)
        if (_rows != rows || _cols != cols){


            if (_cols != cols){
                // __M
                // colPtr_.resize(_rows + 1, 0);
            }
            if (_rows != rows){
                // __M
                colPtr_.resize(_rows + 1, 0);
            }
            _rows = rows;
            _cols = cols;
        }
    }

protected:

    // int to be cholmod compatible!!!!!!!!

    std::vector < int > colPtr_;
    std::vector < int > rowIdx_;
    Vector < ValueType > vals_;

    bool valid_;
    int stype_;
};


/*! SparseMatrix specialized type traits in sparsematrix.cpp */
template< typename ValueType >
void SparseMatrix< ValueType >::copy_(const SparseMapMatrix< double, Index > & S){THROW_TO_IMPL}
template< typename ValueType >
void SparseMatrix< ValueType >::copy_(const SparseMapMatrix< Complex, Index > & S){THROW_TO_IMPL}

template <> DLLEXPORT void SparseMatrix<double>::setVal(Index i, Index j, const double & val);
template <> DLLEXPORT void SparseMatrix<Complex>::setVal(Index i, Index j, const Complex & val);

template <> DLLEXPORT double SparseMatrix<double>::getVal(Index i, Index j, bool warn) const;
template <> DLLEXPORT Complex SparseMatrix<Complex>::getVal(Index i, Index j, bool warn) const;

template <> DLLEXPORT void SparseMatrix<double>::copy_(const SparseMapMatrix< double, Index > & S);
template <> DLLEXPORT void SparseMatrix<Complex>::copy_(const SparseMapMatrix< Complex, Index > & S);

template <> DLLEXPORT void SparseMatrix< double >::
    add(const ElementMatrix < double > & A, const double & f, const double & scale);
template <> DLLEXPORT void SparseMatrix< double >::
    add(const ElementMatrix < double > & A, const Pos & f, const double & scale);
template <> DLLEXPORT void SparseMatrix< double >::
    add(const ElementMatrix < double > & A, const RSmallMatrix & f, const double & scale);

template <> DLLEXPORT void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const Complex & f, const double & scale);
template <> DLLEXPORT void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const Pos & f, const double & scale);
template <> DLLEXPORT void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const CSmallMatrix & f, const double & scale);

inline CSparseMatrix operator + (const CSparseMatrix & A, const RSparseMatrix & B){
    CSparseMatrix ret(A);
    ret.vecVals() += toComplex(B.vecVals());
    return ret;
}
template < class ValueType >
SparseMatrix< ValueType > operator + (const SparseMatrix< ValueType > & A,
                                      const SparseMatrix< ValueType > & B){
    SparseMatrix< ValueType > ret(A);
    return ret += B;
}

template < class ValueType >
SparseMatrix< ValueType > operator - (const SparseMatrix< ValueType > & A,
                                      const SparseMatrix< ValueType > & B){
    SparseMatrix< ValueType > ret(A);
    return ret -= B;
}

template < class ValueType >
SparseMatrix < ValueType > operator * (const SparseMatrix < ValueType > & A,
                                       const ValueType & b){
    SparseMatrix< ValueType > ret(A);
    return ret *= b;
}

template < class ValueType >
SparseMatrix < ValueType > operator * (const ValueType & b,
                                       const SparseMatrix < ValueType > & A){
    SparseMatrix< ValueType > ret(A);
    return ret *= b;
}

inline RVector operator * (const RSparseMatrix & A, const RVector & b){
    return A.mult(b);
}
inline RVector transMult(const RSparseMatrix & A, const RVector & b){
    return A.transMult(b);
}

inline CVector operator * (const CSparseMatrix & A, const CVector & b){
    return A.mult(b);
}
inline CVector operator * (const CSparseMatrix & A, const RVector & b){
    return A.mult(toComplex(b));
}
inline CVector transMult(const CSparseMatrix & A, const CVector & b){
    return A.transMult(b);
}
inline CVector transMult(const CSparseMatrix & A, const RVector & b){
    return A.transMult(toComplex(b));}

inline RSparseMatrix real(const CSparseMatrix & A){
    return RSparseMatrix(A.vecColPtr(), A.vecRowIdx(),
                         real(A.vecVals()), A.stype());
}
inline RSparseMatrix imag(const CSparseMatrix & A){
    return RSparseMatrix(A.vecColPtr(), A.vecRowIdx(),
                         imag(A.vecVals()), A.stype());
}

template <> DLLEXPORT void SparseMatrix< double >::
fillMassMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern);
template <> DLLEXPORT void SparseMatrix< Complex >::
fillMassMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern);

template <> DLLEXPORT void SparseMatrix< double >::
fillStiffnessMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern);
template <> DLLEXPORT void SparseMatrix< Complex >::
fillStiffnessMatrix(const Mesh & mesh, const RVector & a, bool rebuildPattern);

template <> DLLEXPORT void SparseMatrix< double >::
buildSparsityPattern(const Mesh & mesh);
template <> DLLEXPORT void SparseMatrix< Complex >::
buildSparsityPattern(const Mesh & mesh);

template <> DLLEXPORT void SparseMatrix< double >::
buildSparsityPattern(const std::vector < std::set< Index > > & idxMap);
template <> DLLEXPORT void SparseMatrix< Complex >::
buildSparsityPattern(const std::vector < std::set< Index > > & idxMap);

template <> DLLEXPORT void SparseMatrix< double >::
addSparsityPattern(const std::vector < std::set< Index > > & idxMap);
template <> DLLEXPORT void SparseMatrix< Complex >::
addSparsityPattern(const std::vector < std::set< Index > > & idxMap);

template <> DLLEXPORT void SparseMatrix< double >::
reduce(const IVector & ids, bool keepDiag);
template <> DLLEXPORT void SparseMatrix< Complex >::
reduce(const IVector & ids, bool keepDiag);

template <> DLLEXPORT IndexArray SparseMatrix< double >::
createReduceMask(const IVector & ids);

template <> DLLEXPORT IndexArray SparseMatrix< Complex >::
createReduceMask(const IVector & ids);

template <> DLLEXPORT IndexArray SparseMatrix< double >::
createDiagonalMask( );

template <> DLLEXPORT IndexArray SparseMatrix< Complex >::
createDiagonalMask();

} // namespace GIMLI
#endif //GIMLI_SPARSEMATRIX__H
