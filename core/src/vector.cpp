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
#include "meshentities.h"

namespace GIMLI{

#if USE_EIGEN3
template <> Vector< double > &
Vector< double >::operator = (const Eigen::VectorXd & v) {
    this->resize(v.size());
    for (Index i=0; i < v.size(); i++){
        this->data_[i] = v(i);
    }
    return *this;
}
template <> Vector< double > &
Vector< double >::operator += (const Eigen::VectorXd & v) {
    this->resize(v.size());
    for (Index i=0; i < v.size(); i++){
        this->data_[i] += v(i);
    }
    return *this;
}
template <> Vector< double > &
Vector< double >::addVal(const Eigen::VectorXd & v, const IVector & ids){
    ASSERT_EQUAL_SIZE(v, ids)
    for (Index i=0; i < v.size(); i++){
        this->data_[ids[i]] += v(i);
    }
    return *this;
}
template <> Vector< double > &
Vector< double >::setVal(const Eigen::VectorXd & v, const IVector & ids){
    ASSERT_EQUAL_SIZE(v, ids)
    for (Index i=0; i < v.size(); i++){
        this->data_[ids[i]] = v(i);
    }
    return *this;
}

#endif

Index _getColSTep(const ElementMatrix < double > & A){
    //** rows are indices
    //** cols are components
    Index colStep = 1;
    if (A.nCoeff() == 2 and A.cols() == 4){
        // xx, xy, yy
        // A is 2d grad so we ignore d_ij
        colStep = 3; // [0, 3]
    } else if (A.nCoeff() == 2 and A.cols() == 3){
        // xx, yy, xy // 2d mapping
        colStep = 1; // [0, 1]
    } else if (A.nCoeff() == 3 and A.cols() == 9){
        // xx, xy, xz, yx, yy, zy, zx, zy, zz
        // A is 3d grad so we ignore d_ij
        colStep = 4; // [0, 4, 8]
    } else if (A.nCoeff() == 3 and A.cols() == 6){
        // xx, yy, zz, xy, yz, zx // 3d mapping
        colStep = 1;
    }
    return colStep;
}

template<>
void Vector< double >::add(const ElementMatrix < double > & A, bool neg){
    return this->add(A, 1.0, neg);
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const double & scale, bool neg){
    // inuse?
    if (neg == true) THROW_TO_IMPL
    A.integrate();

    //__MS(A.oldStyle(), scale, A)
    if (A.oldStyle()){
        if (A.cols() == 1){
            addVal(A.col(0) * scale, A.rowIDs());
        } else {
            addVal(A.mat().row(0) * scale, A.ids());
        }
    } else {
        // __MS("A.cols()", A.cols(), "A.rows()", A.rows(), "A.nCoeff()", A.nCoeff(), "A.dofPerCoeff()", A.dofPerCoeff())
        // __MS(scale)
        // __MS(A.mat())
        Index maxCols = A.cols();
        Index colStep = _getColSTep(A);
        // switch to A.mat() transpose
        if (A.elastic()) {
            maxCols = A.entity()->dim();
        }
        for (Index i = 0; i < maxCols; i+= colStep){
            for (Index j = 0; j < A.rows(); j++){
                // __MS(i, j, A.rowIDs()[j])
                data_[A.rowIDs()[j]] += A.mat()(j,i) * scale;
            }
        }
    }
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const RVector3 & scale, bool neg){
    // inuse?
    if (neg == true) THROW_TO_IMPL

    // __MS("inuse?")
    if (A.oldStyle()){
        THROW_TO_IMPL
    } else {
        // switch to A.mat() transpose
        A.integrate();
        if (A.elastic()) {
            __MS("elasticity", A.mat().rows(), A.mat().cols())
            THROW_TO_IMPL
        }

        Index colStep = _getColSTep(A);
        for (Index i = 0; i < A.cols(); i+=colStep){
            for (Index j = 0; j < A.rows(); j++){
                data_[A.rowIDs()[j]] += A.mat()(j,i) * scale[i/colStep];
            }
        }
    }
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const RVector & scale){
    // inuse?
    // if (neg == true) THROW_TO_IMPL
    A.integrate();

    if (A.oldStyle()){
        //!! warning this will lead to incorrect results with non constant scale
        //!! use new fea style for correct integration
        // __M
        if (A.cols() == 1){
            // addVal(A.mat().col(0) * scale.get_(A.rowIDs()), A.rowIDs());
            addVal(A.mat().col(0), A.rowIDs(), scale);
        } else {
            //** This[ids] += vals[:] * scale[ids]
            // __MS(A.mat().row(0))
            // addVal(A.mat().row(0) * scale.get_(A.ids()), A.ids());
            addVal(A.mat().row(0), A.ids(), scale);
        }
    } else {
        Index jID = 0;
        if (A.elastic()) {
            __MS("elasticity", A.mat().rows(), A.mat().cols())
            THROW_TO_IMPL
        }

        Index colStep = _getColSTep(A);
        for (Index j = 0; j < A.rows(); j++){
            jID = A.rowIDs()[j];

            if (A.nCoeff() == 1 or scale.size() == A.dofPerCoeff()*A.nCoeff()){
                // scale unsqueezed PosList
                for (Index i = 0; i < A.cols(); i+= colStep){
                    // __MS(i, j, jID, scale[jID], A.mat()(j, i))
                    data_[jID] += A.mat()(j, i) * scale[jID];
                }
            } else if (scale.size() == A.dofPerCoeff()){
                // scale is RVector(nodeCount()) and nCoeff > 1

                for (Index i = 0; i < A.cols(); i+= colStep){
                    data_[jID] += A.mat()(j, i) * scale[jID%A.dofPerCoeff()];
                }
            }
        }
        return;
    }
}
template <>
void Vector< double >::add(const ElementMatrix < double > & A,
                           const RSmallMatrix & scale, bool neg){
    if (A.oldStyle()){
        THROW_TO_IMPL
    } else {
        THROW_TO_IMPL
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


