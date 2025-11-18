/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Thomas Günther thomas@resistivity.net                                    *
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

#include "trans.h"

namespace GIMLI{

template <> DLLEXPORT RVector TransLogLU<RVector>::rangify(const RVector & a) const {
    RVector v(a);
    double minA = min(v);
    double maxA = max(v);
    double lb1 = this->lowerBound() * (1.0 + TRANSTOL);
    if (minA < lb1){
        // log(Verbose, "TransLogLU cap min: ",
        //                 minA, "< lower bound(", this->lowerBound(),")");
        capMin(v, lb1);
    }

    double ub1 = upperbound_ * (1.0 - TRANSTOL);
    if (maxA > ub1 ){
        // log(Verbose, "TransLogLU cap max: ",
        //                 maxA, "> upper bound(", upperbound_,")");
        capMax(v, ub1);
    }
    return v;
}

template <> RVector TransLogLU<RVector>::trans(const RVector & a) const {
    if (std::fabs(upperbound_) < TOLERANCE) return TransLog< RVector >::trans(a);

    RVector tmp = rangify(a);
    return (log(tmp - this->lowerBound()) - log(upperbound_ - tmp));
}

template <> RVector TransLogLU<RVector>::invTrans(const RVector & a) const {
    if (std::fabs(upperbound_) < TOLERANCE) return TransLog< RVector >::invTrans(a);
    RVector expm(a);
    capMax(expm, 50.0);
    expm = exp(expm);
    return (expm * this->upperBound() + this->lowerBound()) / (expm + 1.0);
}

template <> RVector TransLogLU<RVector>::deriv(const RVector & a) const {
    if (std::fabs(upperbound_) < TOLERANCE) return TransLog< RVector >::deriv(a);
    RVector tmp = rangify(a);
    return (1.0 / (tmp - this->lowerBound()) + 1.0 / (this->upperBound() - tmp));
}

typedef RVector Vec;

template <> Vec TransCumulative < Vec >::trans(const Vec & a) const {
    Vec tmp(a.size());
    if (indices_.size() > 0){
        for (Index i = 0; i < transVec_.size(); i ++){
            tmp.setVal(transVec_[i]->trans(a(indices_[i])),
                        indices_[i]);
        }
    } else {
        for (Index i = 0; i < transVec_.size(); i ++){
            tmp.setVal(transVec_[i]->trans(a(slice_[i].first,
                                            slice_[i].second)),
                        slice_[i]);
        }
    }
    return tmp;
}

template <> Vec TransCumulative < Vec >::invTrans(const Vec & a) const {
    Vec tmp(a.size());
    if (indices_.size() > 0){
        for (Index i = 0; i < transVec_.size(); i ++){
            tmp.setVal(transVec_[i]->invTrans(a(indices_[i])),
                        indices_[i]);
        }
    } else {
        for (Index i = 0; i < transVec_.size(); i ++){
            tmp.setVal(transVec_[i]->invTrans(a(slice_[i].first, slice_[i].second)),
                        slice_[i]);
        }
    }
    return tmp;
}

template <> Vec TransCumulative < Vec >::deriv(const Vec & a) const {
    Vec tmp(a.size());
    if (indices_.size() > 0){
        for (Index i = 0; i < transVec_.size(); i ++){
            tmp.setVal(transVec_[i]->deriv(a(indices_[i])),
                        indices_[i]);
        }
    } else {
        for (Index i = 0; i < transVec_.size(); i ++){
            tmp.setVal(transVec_[i]->deriv(a(slice_[i].first, slice_[i].second)),
                        slice_[i]);
        }
    }
    return tmp;
}

template <> void TransCumulative < Vec >::add(Trans< Vec > & trans, Index size) {
    Index start = 0;
    if (!slice_.empty()) start = slice_.back().second;
    this->add(trans, start, start + size);
}
template <> void TransCumulative < Vec >::add(Trans< Vec > & trans, Index start, Index end) {
    // __MS(typeid(&trans).name() << " "<< start <<""<< end)
    transVec_.push_back(&trans);
    slice_.push_back(std::pair< Index, Index >(start, end));
}
template <> void TransCumulative < Vec >::add(Trans< Vec > & trans, const IVector & indices) {
    transVec_.push_back(&trans);
    indices_.push_back(indices);
}


} // namespace GIMLI
