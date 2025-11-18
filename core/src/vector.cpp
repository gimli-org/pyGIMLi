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
#include "vector.h"
#include "elementmatrix.h"

namespace GIMLI{


template<>
void Vector< double >::add(const ElementMatrix < double > & A){
    return this->add(A, 1.0);
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const double & scale){
    if (A.oldStyle()){
        if (A.cols() == 1){
            addVal(A.col(0) * scale, A.rowIDs());
        } else {
            addVal(A.row(0) * scale, A.ids());
        }
    } else {
        A.integrate();
        for (Index i = 0; i < A.cols(); i++){
            for (Index j = 0; j < A.rows(); j++){
                data_[A.rowIDs()[j]] += A.mat()[j][i] * scale;
            }
        }
    }
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const RVector3 & scale){
    if (A.oldStyle()){
        THROW_TO_IMPL
    } else {
        A.integrate();
        for (Index i = 0; i < A.cols(); i++){
            for (Index j = 0; j < A.rows(); j++){
                data_[A.rowIDs()[j]] += A.mat()[j][i] * scale[i];
            }
        }
    }
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const RMatrix & scale){
    if (A.oldStyle()){
        THROW_TO_IMPL
    } else {
        THROW_TO_IMPL
    }
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const RVector & scale){
    if (!A.oldStyle()){
        THROW_TO_IMPL
    }
    A.integrate();
    //!! warning this will lead to incorrect results with non constant scale
    //!! use new fea style for correct integration
    if (A.cols() == 1){
        addVal(A.col(0) * scale.get_(A.rowIDs()), A.rowIDs());
    } else {
        addVal(A.row(0) * scale.get_(A.ids()), A.ids());
    }
}

template <>
void Vector< RVector3 >::clean(){
     if (size_ > 0) {
         for (Index i = 0; i < size_; i ++) {
             data_[i] = RVector3();
         }
     }
}

IndexArray range(Index start, Index stop, Index step){
    IndexArray ret(0);
    for (Index i = start; i < stop; i += step){
        ret.push_back(i);
    }
    return ret;
}

IndexArray range(Index stop){
    return range(0, stop, 1);
}

} // namespace GIMLI{


