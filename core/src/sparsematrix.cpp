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
#include "sparsematrix.h"

namespace GIMLI{

template<>
void SparseMatrix< double >::copy_(const SparseMapMatrix< double, Index > & S){
    this->clear();
    Index col = 0, row = 0;
    double val;
    cols_ = S.cols();
    rows_ = S.rows();
    std::vector < std::map < Index, double > > rowColMap(S.rows());

    for (typename SparseMapMatrix< double, Index>::const_iterator
        it = S.begin(); it != S.end(); it ++){
        row = S.idx1(it);
        col = S.idx2(it);
        val = S.val(it);
        rowColMap[row].insert(std::pair< Index, double >(col, val));
    }
    colPtr_.resize(S.rows() + 1);
    rowIdx_.resize(S.nVals());
    vals_.resize(S.nVals());
    stype_  = S.stype();

    this->colPtr_[0] = 0;

    Index colPtr = 0;
    row = 0;
    for (typename std::vector < std::map < Index, double > >::iterator
         it = rowColMap.begin(); it != rowColMap.end(); it++){
        for (typename std::map < Index, double >::iterator
             itR = (*it).begin(); itR != (*it).end(); itR++){
            rowIdx_[colPtr] = itR->first;
            vals_[colPtr] = (double)itR->second;
            colPtr ++;
        }
        row ++;
        this->colPtr_[row] = colPtr;
    }
    valid_ = true;
}

template<>
void SparseMatrix< Complex >::copy_(const SparseMapMatrix< Complex, Index > & S){
    this->clear();
    Index col = 0, row = 0;
    Complex val;
    cols_ = S.cols();
    rows_ = S.rows();
    std::vector < std::map < Index, Complex > > rowColMap(S.rows());

    for (typename SparseMapMatrix< Complex, Index>::const_iterator
        it = S.begin(); it != S.end(); it ++){
        row = S.idx1(it);
        col = S.idx2(it);
        val = S.val(it);
        rowColMap[row].insert(std::pair< Index, Complex >(col, val));
    }
    colPtr_.resize(S.rows() + 1);
    rowIdx_.resize(S.nVals());
    vals_.resize(S.nVals());
    stype_  = S.stype();

    this->colPtr_[0] = 0;

    Index colPtr = 0;
    row = 0;
    for (typename std::vector < std::map < Index, Complex > >::iterator
         it = rowColMap.begin(); it != rowColMap.end(); it++){
        for (typename std::map < Index, Complex >::iterator
             itR = (*it).begin(); itR != (*it).end(); itR++){
            rowIdx_[colPtr] = itR->first;
            vals_[colPtr] = (Complex)itR->second;
            colPtr ++;
        }
        row ++;
        this->colPtr_[row] = colPtr;
    }
    valid_ = true;
}

template<>
void SparseMapMatrix< double, Index >::copy_(const SparseMatrix< double > & S){
    clear();
    cols_ = S.cols();
    rows_ = S.rows();
    stype_ = S.stype();

    // std::vector < int > colPtr(S.vecColPtr());
    // std::vector < int > rowIdx(S.vecRowIdx());
    // Vector < ValueType > vals(S.vecVals());

    for (Index i = 0; i < S.rows(); i++){
        for (int j = S.vecColPtr()[i]; j < S.vecColPtr()[i + 1]; j ++){
            (*this)[i][S.vecRowIdx()[j]] = S.vecVals()[j];
        }
    }
}
template <> void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, double scale){
    A.integrate();
    double tol = 1e-25;
    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j) * scale;
                this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
            if (::fabs(v) > tol){
                // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << " " << v * scale)
            }
        }
    }
}
template <> void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const Pos & scale){
    A.integrate();
    THROW_TO_IMPL
    // double tol = 1e-25;
    // for (Index i = 0, imax = A.rows(); i < imax; i++){
    //     for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
    //         double v = A.getVal(i, j) * scale;
    //             this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
    //         if (::fabs(v) > tol){
    //             // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << " " << v * scale)
    //         }
    //     }
    // }
}
template <> void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const RMatrix & scale){
    A.integrate();
    THROW_TO_IMPL
    // double tol = 1e-25;
    // for (Index i = 0, imax = A.rows(); i < imax; i++){
    //     for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
    //         double v = A.getVal(i, j) * scale;
    //             this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
    //         if (::fabs(v) > tol){
    //             // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << " " << v * scale)
    //         }
    //     }
    // }
}

template <> void SparseMapMatrix< double, Index >::
    add(const ElementMatrix < double > & A, const Vector < double > & scale){
    A.integrate();
    __MS("inuse?")
    double tol = 1e-25;
    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j) * scale[i];
                this->addVal(A.rowIDs()[i], A.colIDs()[j], v);
            if (::fabs(v) > tol){
            }
        }
    }
}
template <> void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, Complex scale){
    A.integrate();

    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j);
            (*this)[A.rowIDs()[i]][A.colIDs()[j]] += v * scale;
        }
    }
}
template <> void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const Vector < Complex > & scale){
    A.integrate();
    __MS("inuse?")
    for (Index i = 0, imax = A.rows(); i < imax; i++){
        for (Index j = 0, jmax = A.mat().cols(); j < jmax; j++){
            double v = A.getVal(i, j);
            (*this)[A.rowIDs()[i]][A.colIDs()[j]] += v * scale[i];
        }
    }
}
template <> void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const CMatrix & scale){
    A.integrate();
    THROW_TO_IMPL
}
template <> void SparseMapMatrix< Complex, Index >::
    add(const ElementMatrix < double > & A, const Pos & scale){
    A.integrate();
    THROW_TO_IMPL
}

template <> void SparseMapMatrix< double, Index >::
    addToCol(Index id, const ElementMatrix < double > & A, double scale, bool isDiag){
    A.integrate();
    for (Index i = 0, imax = A.size(); i < imax; i++){
            if (isDiag){
                (*this)[A.idx(i)][id] += A.getVal(i, i);
            } else {
                (*this)[A.idx(i)][id] += A.getVal(0, i);
            }
    }
}

template <> void SparseMapMatrix< double, Index >::
    addToRow(Index id, const ElementMatrix < double > & A, double scale, bool isDiag){
    A.integrate();

    for (Index i = 0, imax = A.size(); i < imax; i++){
            if (isDiag){
                (*this)[id][A.idx(i)] += A.getVal(i, i);
            } else {
                (*this)[id][A.idx(i)] += A.getVal(0, i);
            }
        }
}

template <> void SparseMapMatrix< Complex, Index >::
    addToCol(Index id, const ElementMatrix < double > & A, Complex scale, bool isDiag){
        A.integrate();
THROW_TO_IMPL
}
template <> void SparseMapMatrix< Complex, Index >::
    addToRow(Index id, const ElementMatrix < double > & A, Complex scale, bool isDiag){
        A.integrate();
THROW_TO_IMPL
}

template <> void SparseMatrix< double >::add(const ElementMatrix< double > & A,
                                            double scale){
    if (A.oldStyle()){
        if (!valid_) SPARSE_NOT_VALID;
        for (Index i = 0, imax = A.size(); i < imax; i++){
            for (Index j = 0, jmax = A.size(); j < jmax; j++){
                addVal(A.idx(i), A.idx(j), scale * A.getVal(i, j));
            }
        }
    } else {
        A.integrate();
        for (Index i = 0, imax = A.rows(); i < imax; i++){
            for (Index j = 0, jmax = A.cols(); j < jmax; j++){
                // __MS(A.rowIDs()[i] << " " << A.colIDs()[j] << "  "
                //       << scale << " " << A.getVal(i, j))
                addVal(A.rowIDs()[i], A.colIDs()[j], scale * A.getVal(i, j));
            }
        }
    }
}
template <> void SparseMatrix< double >::add(const ElementMatrix< double > & A,
                                            const Pos & scale){
    THROW_TO_IMPL
}
template <> void SparseMatrix< double >::add(const ElementMatrix< double > & A,
                                            const Matrix < double > & scale){
    THROW_TO_IMPL
}

template <> void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, Complex scale){
    if (!valid_) SPARSE_NOT_VALID;
    for (Index i = 0, imax = A.size(); i < imax; i++){
        for (Index j = 0, jmax = A.size(); j < jmax; j++){
            addVal(A.idx(i), A.idx(j), scale * A.getVal(i, j));
        }
    }
}
template <> void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const Pos & scale){
    THROW_TO_IMPL
}
template <> void SparseMatrix< Complex >::
    add(const ElementMatrix < double > & A, const CMatrix & scale){
    THROW_TO_IMPL
}

} // namespace GIMLI{


