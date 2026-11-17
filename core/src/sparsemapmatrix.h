/******************************************************************************
 *   Copyright (C) 2006-2021 by the GIMLi development team                    *
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

#ifndef GIMLI_SPARSEMAPMATRIX__H
#define GIMLI_SPARSEMAPMATRIX__H

#include "gimli.h"
#include "sparsematrix.h"
#include "elementmatrix.h"
#include "vector.h"
#include "vectortemplates.h"
#include "matrix.h"
#include "mesh.h"
#include "meshentities.h"
#include "node.h"
#include "stopwatch.h"

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <iostream>
#include <cmath>

namespace GIMLI{

//! based on: Ulrich Breymann, Addison Wesley Longman 2000 , revised edition ISBN 0-201-67488-2, Designing Components with the C++ STL
template< class ValueType, class IndexType, class ContainerType >
class MatrixElement {
public:
  typedef std::pair< IndexType, IndexType > IndexPair;
  typedef MatrixElement< ValueType, IndexType, ContainerType > & Reference;

  MatrixElement(ContainerType & Cont, IndexType r, IndexType c)
    : C(Cont), I(C.find(IndexPair(r, c))), row(r), column(c) {
  }

  /* An assignment operator is required which in turn requires a
     reference to an object of type MatrixElement. When both the
     left- and right-hand side are identical, nothing has to happen.
     Otherwise, as above, it has to be checked whether the value of
     the right-hand element is 0 or not. The resulting behavior is
     described together with the above assignment operator, so that
     here it is simply called: */

    Reference operator = (Reference rhs) {
        if (this != & rhs) {    // not identical?
            return operator = (rhs.asValue());  // see above
        } return *this;
    }

  /* The constructor initializes the private variables with all
     information that is needed. The container itself is located in
     the sparseMatrix class; here, the reference to it is entered.
     If the passed indices for row and column belong to an element
     not yet stored in the container, the iterator has the value
     C.end(). */

  ValueType asValue() const {
    if (I == C.end()) return ValueType(0); else return (*I).second;
  }

  //** type conversion operator;
  operator ValueType () const { return asValue();  }

  /* According to the definition of the sparse matrix, 0 is returned
     if the element is not present in the container. Otherwise, the
     result is the second part of the object of type value_type
     stored in the container. */

  Reference operator = (const ValueType & x) {
    // not equal 0?
    if (x != ValueType(0) || 1) { // we need the element to force some sought  matrix shape
      /* If the element does not yet exist, it is put,  together
	 with the indices, into an object of type value_type and
	 inserted with insert(): */

        if (I == C.end()) {
            assert(C.size() < C.max_size());
            I = (C.insert(typename ContainerType::value_type(IndexPair(row, column), x))).first;
        } else (*I).second = x;
    }

    /* insert() returns a pair whose first part is an iterator
       pointing to the inserted object. The second part is of type
       bool and indicates whether the insertion took place because
       no element with this key existed. This is, however, not
       evaluated here because, due to the precondition (I ==
       C.end()), the second part must always have the value true.
       If, instead, the element already exists, the value is
       entered into the second part of the value_type object. If
       the value is equal 0, in order to save space the element is
       deleted if it existed. */

    else {                    // x = 0
      if (I != C.end()) {
        C.erase(I);
        I = C.end();
      }
    } return *this;
  }

  Reference operator += (const ValueType & x) {
    if (x != ValueType(0) || 1 ) { // we need the element to force some sought  matrix shape
      if (I == C.end()) {
        assert(C.size() < C.max_size());
        I = (C.insert(typename ContainerType::value_type(IndexPair(row, column), x))).first;
      } else (*I).second += x;
    }
    return *this;
  }

  Reference operator -= (const ValueType & x) {
    if (x != ValueType(0) || 1) {
      if (I == C.end()) {
        assert(C.size() < C.max_size());
        I = (C.insert(typename ContainerType::value_type(IndexPair(row, column), -x))).first;
      } else (*I).second -= x;
    }
    return *this;
  }

private:
  ContainerType & C;
  typename ContainerType::iterator I;
  IndexType row, column;

};  // class MatrixElement


//! based on: Ulrich Breymann, Addison Wesley Longman 2000 , revised edition ISBN 0-201-67488-2, Designing Components with the C++ STL
template< class ValueType, class IndexType >
class SparseMapMatrix : public SparseMatrixBase {
public:
    typedef std::pair< IndexType, IndexType > IndexPair;
    typedef std::map< IndexPair, ValueType, std::less< IndexPair > > ContainerType;
    typedef typename ContainerType::iterator          iterator;
    typedef typename ContainerType::const_iterator    const_iterator;
    typedef MatrixElement< ValueType, IndexType, ContainerType > MatElement;

    /*!stype .. symmetric style. stype=0 (full), stype=1 (UpperRight), stype=2 (LowerLeft)*/
    SparseMapMatrix(IndexType r=0, IndexType c=0, int stype=0)
        : SparseMatrixBase(), stype_(stype) {
        _rows = r;
        _cols = c;
    }

    SparseMapMatrix(const std::string & filename)
        : SparseMatrixBase(), stype_(0) {
        this->load(filename);
    }

    SparseMapMatrix(const SparseMapMatrix< ValueType, IndexType > & S)
        : SparseMatrixBase(){
        clear();
        _rows = S.rows();
        _cols = S.cols();
        stype_ = S.stype();

        for (const_iterator it = S.begin(); it != S.end(); it ++){
            this->setVal(it->first.first, it->first.second, it->second);
        }
    }
    /*!Copy constructor from Matrix. Drops all absolute values lower than droptol. */
    SparseMapMatrix(const Matrix< ValueType > & S,
                    const ValueType & dropTol=0.0)
        : SparseMatrixBase(){
        clear();
        _rows = S.rows();
        _cols = S.cols();
        stype_ = 0;

        ValueType v;
        for (Index i = 0; i < S.rows(); i ++ ){
            for (Index j = 0; j < S.cols(); j ++ ){
                v = S[i][j];
                if (abs(v) > dropTol){
                    this->setVal(i, j, v);
                }
            }
        }
    }
    #ifndef PYGIMLI_CAST // disallow automatic type conversion from python
    SparseMapMatrix(const SparseMatrix< ValueType > & S)
        : SparseMatrixBase(){
        this->copy_(S);
    }
    #endif

    /*! Construct Map Matrix from 3 arrays of the same length.
     *Number of columns are max(j)+1 and Number of rows are max(i)+1.*/
    SparseMapMatrix(const IndexArray & i, const IndexArray & j, const RVector & v)
        : SparseMatrixBase(){
        ASSERT_EQUAL(i.size(), j.size())
        ASSERT_EQUAL(i.size(), v.size())
        stype_ = 0;
        _rows = max(i)+1;
        _cols = max(j)+1;
        for (Index n = 0; n < i.size(); n ++ ) (*this)[i[n]][j[n]] = v[n];
    }

    SparseMapMatrix< ValueType, IndexType > & operator = (const SparseMapMatrix< ValueType, IndexType > & S){
        if (this != &S){
            this->copy(S);
        } return *this;
    }

    SparseMapMatrix & operator = (const SparseMatrix< ValueType > & S){
        this->copy_(S);
        return *this;
    }

    virtual ~SparseMapMatrix() {}

    /*! Return entity rtti value. */
    virtual uint rtti() const { return GIMLI_SPARSE_MAP_MATRIX_RTTI; }

    void resize(Index rows, Index cols){
        _rows = rows;
        _cols = cols;
    }

    void copy_(const SparseMapMatrix< ValueType, IndexType > & S){
        this->clear();
        _rows = S.rows();
        _cols = S.cols();
        stype_ = S.stype();
        for (const_iterator it = S.begin(); it != S.end(); it ++){
            this->setVal(it->first.first, it->first.second, it->second);
        }
    }
    /*! Inline copy.  */
    void copy(const SparseMapMatrix< double, Index > & S){
        this->copy_(S);
    }
    /*! Inline copy.  */
    void copy(const SparseMapMatrix< Complex, Index > & S){
        this->copy_(S);
    }

    void copy_(const SparseMatrix< double > & S);
    void copy_(const SparseMatrix< Complex > & S);
    /*! Inline copy.  */
    void copy(const SparseMatrix< double > & S){ this->copy_(S); }
    /*! Inline copy. */
    void copy(const SparseMatrix< Complex > & S){ this->copy_(S); }


    /*! Add this values to the matrix. */
    inline void add(const IndexArray & rows, const IndexArray & cols,
                    const RVector & vals) {
        ASSERT_EQUAL(vals.size(), rows.size())
        ASSERT_EQUAL(vals.size(), cols.size())

        for (Index i = 0; i < vals.size(); i ++){
            (*this)[rows[i]][cols[i]] += vals[i];
        }
    }

    virtual void clear() {
        C_.clear();
        _cols = 0; _rows = 0; stype_ = 0;
    }

    void cleanRow(IndexType row){
        ASSERT_RANGE(row, 0, this->rows())

        for (auto it = begin(); it != end();){
            if (idx1(it) == row){
                it = C_.erase(it);
            } else {
                ++it;
            }
        }
    }

    void cleanCol(IndexType col){
        ASSERT_RANGE(col, 0, this->cols())
        for (auto it = begin(); it != end();){
            if (idx2(it) == col){
                it = C_.erase(it);
            } else {
                ++it;
            }
        }
    }
    void reduce(const IVector & ids, bool keepDiag=true);

    /*! symmetric type. 0 = unsymmetric, -1 symmetric lower part, 1 symmetric upper part.*/
    inline int stype() const {return stype_;}

    inline void setRows(IndexType r) { _rows = r ; }
    inline void setCols(IndexType c) { _cols = c ; }

    // virtual IndexType cols()     const { return _cols; }
    inline IndexType nRows()     const { return _rows; }
    inline IndexType nCols()     const { return _cols; }

    // inline IndexType size()     const { return C_.size(); }
    inline IndexType max_size() const { return C_.max_size(); }
    inline IndexType nVals()    const { return C_.size(); }

    inline iterator begin() { return C_.begin(); }
    inline iterator end() { return C_.end(); }

    inline const_iterator begin() const { return C_.begin(); }
    inline const_iterator end()   const { return C_.end(); }

    /*!Scale with scale */
    void add(const ElementMatrix < double > & A,
             const ValueType & f=1.0, const ValueType & scale=1.0);
    void add(const ElementMatrix < double > & A,
             const Pos & f, const ValueType & scale=1.0);
    void add(const ElementMatrix < double > & A,
             const SmallMatrix< ValueType> & f, const ValueType & scale=1.0);
    /*!Scale with values from vector scale. Take values from scale[A.ids()]. */
    void add(const ElementMatrix < double > & A,
             const Vector < ValueType > & f, const ValueType & scale=1.0);

    void addToCol(Index id, const ElementMatrix < double > & A,
                  ValueType scale=1.0, bool isDiag=false);
    void addToRow(Index id, const ElementMatrix < double > & A,
                  ValueType scale=1.0, bool isDiag=false);

#define DEFINE_SPARSEMAPMATRIX_UNARY_MOD_OPERATOR__(OP) \
    SparseMapMatrix< ValueType, IndexType > & operator OP##= (const ValueType & v){\
        for (iterator it = begin(); it != end(); it ++) (*it).second OP##= v; \
        return *this; \
    } \

    DEFINE_SPARSEMAPMATRIX_UNARY_MOD_OPERATOR__(+)
    DEFINE_SPARSEMAPMATRIX_UNARY_MOD_OPERATOR__(-)
    DEFINE_SPARSEMAPMATRIX_UNARY_MOD_OPERATOR__(*)
    DEFINE_SPARSEMAPMATRIX_UNARY_MOD_OPERATOR__(/)

#undef DEFINE_SPARSEMMAPATRIX_UNARY_MOD_OPERATOR__

    SparseMapMatrix< ValueType, IndexType > & operator += (const SparseMapMatrix< ValueType, IndexType > & A){
        for (const_iterator it = A.begin(); it != A.end(); it ++){
            this->addVal(it->first.first, it->first.second, it->second);
        }
        return *this;
    }
    SparseMapMatrix< ValueType, IndexType > & operator -= (const SparseMapMatrix< ValueType, IndexType > & A){
        for (const_iterator it = A.begin(); it != A.end(); it ++){
            this->addVal(it->first.first, it->first.second, -it->second);
        }
        return *this;
    }

    //SparseMatrix< T > & operator += (const ElementMatrix < T > & A){
    void operator += (const ElementMatrix < double > & A){
        for (Index i = 0, imax = A.size(); i < imax; i++){
            for (Index j = 0, jmax = A.size(); j < jmax; j++){
                (*this)[A.idx(i)][A.idx(j)] += A.getVal(i, j);
            }
        }
    }

    // ## fails python conversion
    // /*! Negation operator that returns a copy of this with negative values. */
    // SparseMapMatrix< ValueType, IndexType > operator - () {
    //     return *this * ValueType(-1.0); }

    class Aux {  // for index operator below
    public:
        Aux(IndexType r, IndexType maxs, ContainerType & Cont, int stype)
            : _row(r), _maxColumns(maxs), _C(Cont), stype_(stype) { }

        MatElement operator [] (IndexType c) {
            if ((c < 0 || c >= _maxColumns) || (stype_ < 0 && c < _row) || (stype_ > 0 && c > _row)) {
                throwLengthError(
                                  WHERE_AM_I + " idx = " + str(c) + ", " + str(_row) + " maxcol = "
                                  + str(_maxColumns) + " stype: " + str(stype_));
            }
            return MatElement(_C, _row, c);
        }
    protected:
        IndexType _row, _maxColumns;
        ContainerType & _C;
        int stype_;
    };

    Aux operator [] (IndexType r) {
        if (r < 0 || r >= _rows){
            throwLengthError(
                              WHERE_AM_I + " idx = " + str(r) + " maxrow = "
                              + str(_rows));
        }
        return Aux(r, cols(), C_, stype_);
    }

    inline IndexType idx1(const const_iterator & I) const { return (*I).first.first; }
    inline IndexType idx2(const const_iterator & I) const { return (*I).first.second; }

//     inline IndexType idx1(iterator & I) { return (*I).first.first; }
//     inline IndexType idx2(iterator & I) { return (*I).first.second; }

    inline const ValueType & val(const const_iterator & I) const { return (*I).second;  }

    inline ValueType & val(const iterator & I) { return (*I).second;  }

    inline ValueType getVal(IndexType i, IndexType j) { return (*this)[i][j]; }

    inline void setVal(IndexType i, IndexType j, const ValueType & val) {
        if ((stype_ < 0 && i > j) || (stype_ > 0 && i < j)) return;

        if (i >= _rows) _rows = i+1;
        if (j >= _cols) _cols = j+1;
        //__MS(i,j,val)
        (*this)[i][j] = val;
        // if ((i >= 0 && i < _rows) && (j >=0 && j < _cols)) {
        // } else {
        //     throwLengthError(
        //                       WHERE_AM_I +
        //                       " i = " + str(i) + " max_row = " + str(_rows) +
        //                       " j = " + str(j) + " max_col = " + str(_cols)
        //                      );
        // }
    }
    inline void addVal(IndexType i, IndexType j, const ValueType & val) {
        if ((stype_ < 0 && i > j) || (stype_ > 0 && i < j)) return;
        if (i >= _rows) _rows = i+1;
        if (j >= _cols) _cols = j+1;
        (*this)[i][j] += val;

        // if ((i >= 0 && i < _rows) && (j >=0 && j < _cols)) {
        //     (*this)[i][j] += val;
        // } else {
        //     throwLengthError(
        //                       WHERE_AM_I +
        //                       " i = " + str(i) + " max_row = " + str(_rows) +
        //                       " j = " + str(j) + " max_col = " + str(_cols)
        //                      );
        // }
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
    /*! ret = [this * ax, this * ay, this * az. with ret = r3 and a is sqeezed pos vector. a = [ax, ay, az] */
    virtual void mult(const Vector < ValueType > & a,
                      Vector < Pos > & ret) const;
    /*! ret = [this.T * ax, this.T * ay, this.T * az. with ret = r3 and a is sqeezed pos vector. a = [ax, ay, az] */
    virtual void transMult(const Vector < ValueType > & a,
                           Vector < Pos > & ret) const;

    Vector < ValueType > col(const Index c) const {
        ASSERT_RANGE(c, 0, cols())
        Vector < ValueType > b(rows(), 0.0);
        for (const_iterator it = this->begin(); it != this->end(); it ++){
            if (idx2(it) == c){
                b[idx1(it)] = val(it);
            }
        }
        return b;
        // Vector < ValueType > null(this->cols(), 0.0);
        // null[i] = 1.0;
        // return this->mult(null);
    }

    Vector < ValueType > row(const Index r) const {
        ASSERT_RANGE(r, 0, rows())
        Vector < ValueType > b(cols(), 0.0);
        for (const_iterator it = this->begin(); it != this->end(); it ++){
            if (idx1(it) == r){
                b[idx2(it)] = val(it);
            }
        }

        // for (int col = colPtr_[r]; col < colPtr_[r + 1]; col ++){
        //     b[rowIdx_[col]] = vals_[col];
        // }
        return b;

        // Vector < ValueType > null(this->rows(), 0.0);
        // null[i] = 1.0;
        // return this->transMult(null);
    }

    void save(const std::string & filename) const {
        std::fstream file; openOutFile(filename, &file);

        for (const_iterator it = begin(); it != end(); it ++){
            file <<  idx1(it) << " " << idx2(it) << " " << val(it) << std::endl;
        }

        file.close();
    }

    void load(const std::string & filename){
        std::fstream file; openInFile(filename, &file);
        std::vector < IndexType > vi, vj;
        std::vector < ValueType > vval;
        IndexType i, j;
        ValueType val;
        while (file >> i >> j >> val){
            vi.push_back(i);
            vj.push_back(j);
            vval.push_back(val);
        }
        file.close();
        setRows(IndexType(max(vi) + 1));
        setCols(IndexType(max(vj) + 1));

        for (Index i = 0; i < vi.size(); i ++){
            (*this)[vi[i]][vj[i]] = vval[i];
        }
    }

    /*! Import columnwise from bmat starting at offset */
    void importCol(const std::string & filename, double dropTol, Index colOffset){
    //std::cout << "rows: " << Jcluster.rows() << " cols: " << Jcluster.cols()<< std::endl;

        FILE *file; file = fopen(filename.c_str(), "r+b");
        if (!file) {
            throwError(WHERE_AM_I + " " + filename + ": " + strerror(errno));
        }
        Index ret = 0;


        uint32 rows = 0;
        ret = fread(&rows, sizeof(uint32), 1, file);
        if (ret == 0) throwError("fail reading file " + filename);
        uint32 cols = 0;
        ret = fread(&cols, sizeof(uint32), 1, file);
        if (ret == 0) throwError("fail reading file " + filename);
        ValueType val;
        for (uint i = 0; i < rows; i ++){
            for (uint j = 0; j < cols; j ++){
                ret = fread(&val, sizeof(ValueType), 1, file);
                if (ret == 0) throwError("fail reading file " + filename);
                if (abs(val) > dropTol) this->setVal(i, j + colOffset, val);
            }
        }

        fclose(file);
    }
    // no default arg here .. pygimli@win64 linker bug
    void importCol(const std::string & filename, double dropTol=1e-3){
        importCol(filename, dropTol, 0);
    }

    /*! Fill existing arrays with values, row and column indieces of this
     * SparseMapMatrix*/
    void fillArrays(Vector < ValueType > & vals, IndexArray & rows, IndexArray & cols){
        vals.resize(C_.size());
        rows.resize(C_.size());
        cols.resize(C_.size());
        Index i = 0;

        for (const_iterator it = this->begin(); it != this->end(); it ++, i ++){
            rows[i] = idx1(it);
            cols[i] = idx2(it);
            vals[i] = val(it);
        }
    }

    Vector< ValueType > values() const {
        Vector< ValueType > ret(C_.size());
        Index i = 0;
        for (const_iterator it = this->begin(); it != this->end(); it ++, i ++){
            ret[i] = val(it);
        }
        return ret;
    }
    IVector rowIDs() const {
        IVector ret(C_.size());
        Index i = 0;
        for (const_iterator it = this->begin(); it != this->end(); it ++, i ++){
            ret[i] = idx1(it);
        }
        return ret;
    }
    IVector colIDs() const {
        IVector ret(C_.size());
        Index i = 0;
        for (const_iterator it = this->begin(); it != this->end(); it ++, i ++){
            ret[i] = idx2(it);
        }
        return ret;
    }

protected:

//   IndexType _rows, _cols;
  ContainerType C_;
  // 0 .. nonsymmetric, -1 symmetric lower part, 1 symmetric upper part
  int stype_;
};// class SparseMapMatrix


/*! Scales a matrix A from left and right vectors such that A -> diag(l) * A * diag(r) */
template< class Vec >
void scaleMatrix(SparseMapMatrix< double, Index > & S,
                 const Vec & l, const Vec & r) {

    if (S.cols() != r.size())
        throwLengthError(WHERE_AM_I + " " + str(S.cols())
                                            + " != " + str(r.size()));
    if (S.rows() != l.size())
        throwLengthError(WHERE_AM_I + " " + str(S.rows())
                                            + " != " + str(l.size()));

    for (SparseMapMatrix< double, Index >::iterator it = S.begin(); it != S.end(); it ++){
                S.val(it) *= l[S.idx1(it)] * r[S.idx2(it)];
//       int i = S.idx1(it);
//     int j = S.idx2(it);
//     S[i][j] *= (l[i] * r[j]);
    }
    return;
}

/*! Performs a rank 1 update of a matrix such that A -> A + u * v^T */
template< class Vec >
void rank1Update(SparseMapMatrix< double, Index > & S, const Vec & u, const Vec & v) {

    if (S.cols() != v.size())
        throwLengthError(WHERE_AM_I + " " + str(S.cols())
                                + " != " + str(v.size()));
    if (S.rows() != u.size())
        throwLengthError(WHERE_AM_I + " " + str(S.rows())
                                + " != " + str(u.size()));

    for (SparseMapMatrix< double, Index >::iterator it = S.begin(); it != S.end(); it ++){
        S.val(it) += u[S.idx1(it)] * v[S.idx2(it)];
    }
    return;
}

template < class ValueType, class IndexType >
void save(const SparseMapMatrix< ValueType, IndexType > & S,
         const std::string & fname){
    S.save(fname);
}

template < class ValueType, class IndexType >
int load(SparseMapMatrix< ValueType, IndexType > & S,
         const std::string & fname){
    return S.load(fname);
}

inline RSparseMapMatrix real(const CSparseMapMatrix & A){
    RSparseMapMatrix R(A.rows(), A.cols());
    for (CSparseMapMatrix::const_iterator it = A.begin(); it != A.end(); it ++){
        R.setVal(it->first.first, it->first.second, it->second.real());
    }
    return R;
}

inline RSparseMapMatrix imag(const CSparseMapMatrix & A){
    RSparseMapMatrix R(A.rows(), A.cols());
    for (CSparseMapMatrix::const_iterator it = A.begin(); it != A.end(); it ++){
        R.setVal(it->first.first, it->first.second, it->second.imag());
    }
    return R;
}

inline RVector operator * (const RSparseMapMatrix & A, const RVector & b){
    return A.mult(b);
}
inline CVector operator * (const CSparseMapMatrix & A, const CVector & b){
    return A.mult(b);
}
inline CVector operator * (const CSparseMapMatrix & A, const RVector & b){
    return A.mult(toComplex(b));
}
inline RSparseMapMatrix operator + (const RSparseMapMatrix & A,
                                    const RSparseMapMatrix & B){
    RSparseMapMatrix tmp(A);
    return tmp += B;
}

inline RSparseMapMatrix operator - (const RSparseMapMatrix & A,
                                    const RSparseMapMatrix & B){
    RSparseMapMatrix tmp(A);
    return tmp -= B;
}

#define DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(OP) \
    inline DLLEXPORT RSparseMapMatrix operator OP (const RSparseMapMatrix & A, \
                                            const double & v){\
        return RSparseMapMatrix(A) OP##= v; \
    } \
    inline DLLEXPORT RSparseMapMatrix operator OP (const RSparseMapMatrix & A, \
                                                   long v){\
        return RSparseMapMatrix(A) OP##= v; \
    } \

    DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(+)
    DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(-)
    DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(*)
    DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(/)

#undef DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__

#define DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(OP) \
inline DLLEXPORT RSparseMapMatrix operator OP (const double & v, \
                                     const RSparseMapMatrix & A){\
    return RSparseMapMatrix(A) OP##= v; \
} \
inline DLLEXPORT RSparseMapMatrix operator OP (long v, \
                                     const RSparseMapMatrix & A){\
    return RSparseMapMatrix(A) OP##= v; \
} \

DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(+)
DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__(*)

#undef DEFINE_SPARSEMAPMATRIX_EXPR_OPERATOR__

/*! SparseMapMatrix specialized type traits in sparsematrix.cpp */
template< typename ValueType, typename Index >
void SparseMapMatrix< ValueType, Index >::
    copy_(const SparseMatrix< double > & S){THROW_TO_IMPL}

template< typename ValueType, typename Index >
void SparseMapMatrix< ValueType, Index >::
    copy_(const SparseMatrix< Complex > & S){THROW_TO_IMPL}

template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    copy_(const SparseMatrix< double > & S);

// template <> DLLEXPORT void SparseMapMatrix< double, Index >::
//     mult(const Vector < double > & a, Vector < double > & ret) const;
// template <> DLLEXPORT void SparseMapMatrix< double, Index >::
//     transMult(const Vector < double > & a, Vector < double > & ret) const;

// template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
//     mult(const Vector < Complex > & a, Vector < Complex > & ret) const;
// template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
//     transMult(const Vector < Complex > & a, Vector < Complex > & ret) const;

template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    mult(const RVector & a, Vector < Pos > & ret) const;
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    mult(const CVector & a, Vector < Pos > & ret) const;
template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    transMult(const RVector & a, Vector < Pos > & ret) const;
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    transMult(const CVector & a, Vector < Pos > & ret) const;

template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const double & f,
        const double & scale);
template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const Pos & f,
        const double & scale);
template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const RSmallMatrix & f,
        const double & scale);
template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const Vector < double > & f,
        const double & scale);

template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const Complex & f,
        const Complex & scale);
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const Pos & f,
        const Complex & scale);
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const CSmallMatrix & f,
        const Complex & scale);
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const Vector < Complex > & f,
        const Complex & scale);

template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    addToCol(Index id, const ElementMatrix < double > & A, double scale,
             bool isDiag);
template <> DLLEXPORT void SparseMapMatrix< double, Index >::
    addToRow(Index id, const ElementMatrix < double > & A, double scale,
             bool isDiag);
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    addToCol(Index id, const ElementMatrix < double > & A, Complex scale,
             bool isDiag);
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >::
    addToRow(Index id, const ElementMatrix < double > & A, Complex scale,
             bool isDiag);

template <> DLLEXPORT void SparseMapMatrix< double, Index >
    ::reduce(const IVector & ids, bool keepDiag);
template <> DLLEXPORT void SparseMapMatrix< Complex, Index >
    ::reduce(const IVector & ids, bool keepDiag);

/*!Convenience functions.*/

inline RVector transMult(const RSparseMapMatrix & A, const RVector & b){
    return A.transMult(b);
}
inline CVector transMult(const CSparseMapMatrix & A, const CVector & b){
    return A.transMult(b);
}
inline CVector transMult(const CSparseMapMatrix & A, const RVector & b){
    return A.transMult(toComplex(b));
}

/*! RVector[i] = A[i]*ret[i] */
DLLEXPORT void
    mult(const std::vector < RSparseMapMatrix > & A,
         const RVector & b, std::vector< RVector > & ret);

/*! PosVector[i] = [A[i]*ret[i][0..dof], A[i]*ret[i][dof..2*dof], ..*/
DLLEXPORT void
    mult(const std::vector < RSparseMapMatrix > & A,
         const RVector & b, std::vector< PosVector > & ret);

DLLEXPORT std::vector< RVector >
    mult(const std::vector < RSparseMapMatrix > & A,
         const RVector & b);


} // namespace GIMLI

#endif //GIMLI_SPARSEMAPMATRIX__H
