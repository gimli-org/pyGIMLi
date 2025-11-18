/******************************************************************************
 *   Copyright (C) 2005-2024 by the GIMLi development team                    *
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
#include "blockmatrix.h"

namespace GIMLI{

template <> Vector < double > BlockMatrix< double >::mult(const Vector < double >  & b) const{
    // no need to check here .. let the matrices itself check
    // if (b.size() != this->cols()){
    //     throwLengthError(WHERE_AM_I + " wrong size of vector b (" +
    //     str(b.size()) + ") needed: " + str(this->cols()));
    // }

    RVector ret(rows_);

    for (Index i = 0; i < entries_.size(); i ++ ){
        BlockMatrixEntry entry = entries_[i];

        MatrixBase *mat = matrices_[entry.matrixID];

        ret.addVal(mat->mult(b.getVal(entry.colStart,
                                        entry.colStart + mat->cols())) *
                    entry.scale,
                    entry.rowStart, entry.rowStart + mat->rows());
    }

    return ret;
}
template <> Vector < double >
BlockMatrix< double >::transMult(const Vector < double > & b) const{
        // no need to check here .. let the matrices itself check
    // if (b.size() != this->rows()){
    //     throwLengthError(WHERE_AM_I + " wrong size of vector b (" +
    //     str(b.size()) + ") needed: " + str(this->rows()));
    // }

    RVector ret(cols_);
        for (Index i = 0; i < entries_.size(); i++){
        BlockMatrixEntry entry = entries_[i];

        MatrixBase *mat = matrices_[entry.matrixID];

        ret.addVal(mat->transMult(b.getVal(entry.rowStart,
                                            entry.rowStart + mat->rows())) *
                    entry.scale,
                    entry.colStart, entry.colStart + mat->cols());
    }
    return ret;
}
template <>
RSparseMapMatrix BlockMatrix< double >::sparseMapMatrix() const {

    RSparseMapMatrix ret(this->rows(), this->cols());

    for (Index i = 0; i < entries_.size(); i++){
        BlockMatrixEntry entry(entries_[i]);

        MatrixBase *mat = matrices_[entry.matrixID];

        RVector vals(0);
        IndexArray rows(0);
        IndexArray cols(0);

        switch (mat->rtti()){
            case GIMLI_SPARSE_CRS_MATRIX_RTTI:{
                RSparseMapMatrix S(*dynamic_cast< RSparseMatrix * >(mat));
                S.fillArrays(vals, rows, cols);
            }  break;
            case GIMLI_SPARSE_MAP_MATRIX_RTTI:
                dynamic_cast< RSparseMapMatrix * >(mat)->fillArrays(vals, rows, cols);
                break;
            default:
                log(Critical, "Matrix type need to be either SparseMatrix or SparseMapMatrix");
                return ret;
        }
        ret.add(rows + entry.rowStart, cols + entry.colStart, vals * entry.scale);
    }
    return ret;
}

} // namespace GIMLI

