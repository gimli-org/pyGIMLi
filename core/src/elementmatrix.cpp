/******************************************************************************
 *   Copyright (C) 2006-2025 by the GIMLi development team                    *
 *   Carsten Rücker carsten@gimli.org
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

#include "elementmatrix.h"
#include "shape.h"
#include "meshentities.h"
#include "node.h"
#include "pos.h"
#include "sparsematrix.h"

#include "integration.h"


namespace GIMLI{

template < >
std::ostream & operator << (std::ostream & str,
                            const ElementMatrix< double > & e){
    for (Index i = 0; i < e.colIDs().size(); i ++) str << e.colIDs()[i] << " " ;

    str << std::endl;
    for (Index i = 0; i < e.size(); i ++){
        str << e.rowIDs()[i] << "\t: ";
        for (Index j = 0; j < e.colIDs().size(); j ++){
            str << e.getVal(i, j) << " ";
        }
        str << std::endl;
    }
    return str;
}

/*
// #if USE_EIGEN3
// #define DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(OP) \
// template < > ElementMatrix < double > & \
// ElementMatrix < double >::operator OP##= (double val) { \
//     if (this->_newStyle){ \
//         if (this->_integrated){ \
//             mat_.array() OP##= val; \
//         } \
//         for (auto & m: _matX){ \
//             m.array() OP##= val; \
//         } \
//         return *this;\
//     } \
//     mat_.array() OP##= val; \
//     return *this; \
// } \
*/
// DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(+)
// DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(-)
// DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(/)
// DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(*)
// #else

#define DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(OP) \
template < > ElementMatrix < double > & \
ElementMatrix < double >::operator OP##= (double val) { \
    if (this->_newStyle){ \
        if (this->_integrated){ \
            mat_ OP##= val; \
        } \
        for (auto & m: _matX){ \
            m OP##= val; \
        } \
        return *this;\
    } \
    mat_ OP##= val; \
    return *this;\
} \

DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(+)
DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(-)
DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(/)
DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(*)

#undef DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__
// #endif

template < >  ElementMatrix < double > &
ElementMatrix < double >::operator += (const ElementMatrix < double > & E){
    ASSERT_EQUAL_SIZE(this->_matX, E.matX())

    for (Index i=0; i < this->_matX.size(); i ++){
        if (this->_matX[i].rows() == E.matX()[i].rows() &&
            this->_matX[i].cols() == E.matX()[i].cols()){
            this->_matX[i] += E.matX()[i];
        } else {
            THROW_TO_IMPL
        }
    }
    if (1 || this->isIntegrated()){ // force allways to stay compatible to oldfea
        if (this->rows() == E.rows() && this->cols() == E.cols()){
            this->mat_ += E.mat();
        } else {
            THROW_TO_IMPL
        }
        // we could cross check if manual integration works too
    }
    return *this;
}
// template < > DLLEXPORT const Vector< double > &
// ElementMatrix < double >::operator[](Index row) const {
//     __MS("inuse?")

//     return mat_[row];
//     // return mat_.row(row);
// }

/*! Set specific columns to zero. */
template < > void
ElementMatrix < double >::cleanCols(const IndexArray & c){
    this->mat_.cleanCols(c);
}

template < > DLLEXPORT Vector< double >
ElementMatrix < double >::row_RM(Index i) const {
    // #if USE_EIGEN3
    //     Vector< double > ret(this->cols());
    //     for (Index j = 0; j < ret.size(); j ++){
    //         ret[j] = mat_(i,j);
    //     }
    //     return ret;
    //     //__MS("eigen3")
    // #else
        return this->mat_.row(i);
    // #endif
}
template < > DLLEXPORT Vector< double >
ElementMatrix < double >::col(Index i) const {
    Vector< double > ret(this->rows());
    for (Index j = 0; j < ret.size(); j ++){
        ret[j] = mat_(j,i);
    }
    return ret;
}

template < > DLLEXPORT void
ElementMatrix < double >::setMatXI(Index i, const RSmallMatrix & m){
    ASSERT_RANGE(i, 0, this->_matX.size())
    //this->_integrated = false; ## better here but lead to problems until integration is marked correctly for ref imnpl.
    this->_matX[i] = m;
}
// template < > DLLEXPORT void
// ElementMatrix < double >::setMatXI_RM(Index i, const RSmallMatrix  & m){
// // #if USE_EIGEN3
// //     toEigenMatrix(m, this->_matX[i]);
// // #else
//     return this->setMatXI(i, m);
// // #endif
// }
// template < > DLLEXPORT void
// ElementMatrix < double >::setMat_RM(const RSmallMatrix  & m){
// // #if USE_EIGEN3
// // __MS(m)
//     // toEigenMatrix(m, this->mat_);
// // __MS(this->mat_)
// // #else
//     return this->setMat(m);
// // #endif
// }
// template < > DLLEXPORT RSmallMatrix
// ElementMatrix < double >::mat_RM() const {
//     // __MS(this->mat().rows(), this->mat().cols() )
// // #if USE_EIGEN3
// //     RSmallMatrix  r;
// //     toRSmallMatrix (this->mat(), r);
// //     return r;
// // #else
//     return this->mat();
// #endif
// }


template < > DLLEXPORT void
ElementMatrix < double >::setX(const PosVector & p) {
    _x = &p;

    // __MS('x', _matX.size(),_x->size())
    if (_matX.size() != _x->size()){
        _matX.resize(_x->size());
    }

}
template < > DLLEXPORT void
ElementMatrix < double >::setW(const RVector & w) {
    _w = &w;
    // __MS('w', _matX.size(),_w->size())
    if (_matX.size() != _w->size()){
        _matX.resize(_w->size());
    }
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::u(const MeshEntity & ent,
                            const RVector & w,
                            const PosVector & x,
                            bool verbose){

    uint nVerts = ent.nodeCount();
    std::map< uint, RVector >::const_iterator it = uCache_.find(ent.rtti());

    if (it == uCache_.end()) {
        uint nRules = w.size();

        RVector u(nVerts);
        RSmallMatrix  N(nVerts, nRules);

        RVector tmp;
        for (uint i = 0; i < nRules; i ++){
            tmp = ent.N(x[i]);
            N.setCol(i, tmp);
        }
        for (uint i = 0; i < nVerts; i ++){
            u[i] = sum(w * N[i]);
        }
        uCache_[ent.rtti()] = u;
        it = uCache_.find(ent.rtti());
    }

    double A = ent.shape().domainSize();
    // double J = det(ent.shape().createJacobian());
    // __MS(A, J)
    for (Index i = 0; i < nVerts; i ++){
        // __MS(i, A, it->second[i])
        mat_(0,i) = A * it->second[i];
        if (this->_nDof > 0){
            if (ent.dim() == 2){
                // #if USE_EIGEN3
                //     __MS("EIGENCHECK")
                //     //mat_(nVerts) = setVal(mat_[0][i], nVerts + i);
                // #else
                mat_[nVerts].setVal(mat_[0][i], nVerts + i);
                // #endif
            }
            if (ent.dim() == 3){
                // #if USE_EIGEN3
                //     __MS("EIGENCHECK")
                // #else
                mat_[2*nVerts].setVal(mat_[0][i], 2*nVerts + i);
                // #endif
            }
        }
    }

    if (verbose) std::cout << "int u " << *this << std::endl;
    return *this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::u2(const MeshEntity & ent,
                                                        const RVector & w,
                                                        const PosVector & x,
                                                        bool verbose){

    uint nVerts = ent.nodeCount();
    std::map< uint, RSmallMatrix >::const_iterator it = u2Cache_.find(ent.rtti());
    // const auto &u2Cache = u2Cache_;
    // auto it = u2Cache.find(ent.rtti());

    if (it == u2Cache_.end()) {
        // __M
        Index nRules = w.size();

        RSmallMatrix  u2(nVerts, nVerts);
        RSmallMatrix  N(nVerts, nRules);

        RVector tmp;
        for (Index i = 0; i < nRules; i ++){
            tmp = ent.N(x[i]);
            N.setCol(i, tmp);
        }

        // print(N);

        double t;
        for (uint i = 0; i < nVerts; i ++){
            for (uint j = i; j < nVerts; j ++){

                // print(i, N[i]);
                // print(j, N[j]);
                // print('*', N[i],N[j]);
                // print('*', N[i]*N[j]);
                // print("*********************************************");

                t = sum(w * N[j] * N[i]);

                u2[i][j] = t;
                u2[j][i] = t;

                // print(i, j, t, u2[i][j]);
            }
                // exit(1);
        }
        // print(u2);

        u2Cache_[ent.rtti()] = u2;
        it = u2Cache_.find(ent.rtti());
    }

    double A = ent.shape().domainSize();
//** very slow yet, pimp this with expressions (matrix::operator = (matExpression &))
//    mat_ = it->second * A;
    // mat_ = it->second * A;

    for (Index i = 0; i < nVerts; i ++){
        for (Index j = 0; j < nVerts; j ++){
            // print(i, j, A * it->second(i,j));
            mat_(i, j) = A * it->second(i,j);
        }
    }

    if (verbose) std::cout << "int u2 " << *this << std::endl;

    return *this;
}

template <> DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::dudi(
                                                        const MeshEntity & ent,
                                                        const RVector & w,
                                                        const PosVector & x,
                                                        Index dim,
                                                        bool verbose){

    this->fillIds(ent);
    Index nVerts = ent.nodeCount();
    Index nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        for (Index i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
            if (ent.dim() > 1){
                dNds_.resize(nVerts, nRules);
                dNds_.setCol(i, ent.dNdL(x[i], 1));
            }
            if (ent.dim() > 2){
                dNdt_.resize(nVerts, nRules);
                dNdt_.setCol(i, ent.dNdL(x[i], 2));
            }
        }
        dNdx_.resize(nVerts, nRules);
    }

    double drdi = ent.shape().drstdxyz(0, dim);
    double dsdi = ent.shape().drstdxyz(1, dim);
    double dtdi = ent.shape().drstdxyz(2, dim);

    // double A = ent.shape().domainSize();

    for (Index i = 0; i < nVerts; i ++){
        switch (ent.dim()){
            case 1: dNdx_[i].assign(drdi * dNdr_[i]); break;
            case 2: dNdx_[i].assign(drdi * dNdr_[i] + dsdi * dNds_[i]); break;
            case 3: dNdx_[i].assign(drdi * dNdr_[i] + dsdi * dNds_[i] + dtdi * dNdt_[i]); break;
        }

        mat_(i,i) = sum(w * dNdx_[i]);
    }
    if (verbose) std::cout << "int dudx " << *this << std::endl;
    return *this;
}


template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::ux2(const MeshEntity & ent,
                                                         const RVector & w,
                                                         const PosVector & x,
                                                         bool verbose){

    uint nVerts = ent.nodeCount();
    uint nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        for (uint i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
        }
        dNdx_.resize(nVerts, nRules);
    }

    double drdx = ent.shape().drstdxyz(0, 0);
//     double drdx = ent.shape().invJacobian()[0];

    double A = ent.shape().domainSize();

    for (uint i = 0; i < nVerts; i ++){
        dNdx_[i].assign(drdx * dNdr_[i]);
    }
    for (uint i = 0; i < nVerts; i ++){
        for (uint j = i; j < nVerts; j ++){
            mat_(i,j) = A * sum(w * (dNdx_[i] * dNdx_[j]));
            mat_(j,i) = mat_(i,j);
        }
    }
    if (verbose) std::cout << "int ux2uy2 " << *this << std::endl;
    return *this;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::ux2uy2(const MeshEntity & ent,
                                 const RVector & w,
                                 const PosVector & x,
                                 bool verbose){

//     __MS(w)

    Index nVerts = ent.nodeCount();
    Index nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        dNds_.resize(nVerts, nRules);

        for (uint i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
            dNds_.setCol(i, ent.dNdL(x[i], 1));

//             __MS(i << " " << dNdr_[i])
        }

        dNdx_.resize(nVerts, nRules);
        dNdy_.resize(nVerts, nRules);
    }

//     double drdx = ent.shape().invJacobian()[0]; //ent.shape().drstdxyz(0,0)
//     double drdy = ent.shape().invJacobian()[1];
//     double dsdx = ent.shape().invJacobian()[2];
//     double dsdy = ent.shape().invJacobian()[3];

    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);

    double A = ent.shape().domainSize();
    for (Index i = 0; i < nVerts; i ++){
        dNdx_[i].assign(drdx * dNdr_[i] + dsdx * dNds_[i]);
        dNdy_[i].assign(drdy * dNdr_[i] + dsdy * dNds_[i]);
    }

    for (Index i = 0; i < nVerts; i ++){
        // dNidx_.assign(drdx * dNdr_[i] + dsdx * dNds_[i]);
        // dNidy_.assign(drdy * dNdr_[i] + dsdy * dNds_[i]);

        for (Index j = i; j < nVerts; j ++){
            mat_(i,j) = A * sum(w * (dNdx_[i] * dNdx_[j] + dNdy_[i] * dNdy_[j]));

//             mat_(i,j) = A * sum(w * (dNidx_ * (drdx * dNdr_[j] + dsdx * dNds_[j]) +
//                                       dNidy_ * (drdy * dNdr_[j] + dsdy * dNds_[j])));
//             mat_(i,j) = A * sum(w * ((drdx * dNdr_[i] + dsdx * dNds_[i]) * (drdx * dNdr_[j] + dsdx * dNds_[j]) +
//                                       (drdy * dNdr_[i] + dsdy * dNds_[i]) * (drdy * dNdr_[j] + dsdy * dNds_[j])));
            mat_(j,i) = mat_(i,j);
        }
    }

    if (verbose) std::cout << "int ux2uy2 " << *this << std::endl;
    return *this;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::ux2uy2uz2(const MeshEntity & ent,
                                    const RVector & w,
                                    const PosVector & x,
                                    bool verbose){

    Index nVerts = ent.nodeCount();
    Index nRules = w.size();

    if (dNdr_.rows() != nVerts){
        dNdr_.resize(nVerts, nRules);
        dNds_.resize(nVerts, nRules);
        dNdt_.resize(nVerts, nRules);

        for (Index i = 0; i < nRules; i ++){
            dNdr_.setCol(i, ent.dNdL(x[i], 0));
            dNds_.setCol(i, ent.dNdL(x[i], 1));
            dNdt_.setCol(i, ent.dNdL(x[i], 2));
        }

        dNdx_.resize(nVerts, nRules);
        dNdy_.resize(nVerts, nRules);
        dNdz_.resize(nVerts, nRules);
    }

    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double drdz = ent.shape().drstdxyz(0, 2);

    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);
    double dsdz = ent.shape().drstdxyz(1, 2);

    double dtdx = ent.shape().drstdxyz(2, 0);
    double dtdy = ent.shape().drstdxyz(2, 1);
    double dtdz = ent.shape().drstdxyz(2, 2);

    // double drdx = ent.shape().invJacobian()[0];
    // double drdy = ent.shape().invJacobian()[1];
    // double drdz = ent.shape().invJacobian()[2];
    //
    // double dsdx = ent.shape().invJacobian()[3];
    // double dsdy = ent.shape().invJacobian()[4];
    // double dsdz = ent.shape().invJacobian()[5];
    //
    // double dtdx = ent.shape().invJacobian()[6];
    // double dtdy = ent.shape().invJacobian()[7];
    // double dtdz = ent.shape().invJacobian()[8];

    double A = ent.shape().domainSize();
    for (Index i = 0; i < nVerts; i ++){
        dNdx_[i].assign(drdx * dNdr_[i] + dsdx * dNds_[i] + dtdx * dNdt_[i]);
        dNdy_[i].assign(drdy * dNdr_[i] + dsdy * dNds_[i] + dtdy * dNdt_[i]);
        dNdz_[i].assign(drdz * dNdr_[i] + dsdz * dNds_[i] + dtdz * dNdt_[i]);
    }

    for (Index i = 0; i < nVerts; i ++){
        for (Index j = i; j < nVerts; j ++){
            mat_(i,j) = A * sum(w * (dNdx_[i] * dNdx_[j] + dNdy_[i] * dNdy_[j] + dNdz_[i] * dNdz_[j]));

//             mat_(i,j) = A * sum(w * ((drdx * dNdr_[i] + dsdx * dNds_[i] + dtdx * dNdt_[i]) *
//                                       (drdx * dNdr_[j] + dsdx * dNds_[j] + dtdx * dNdt_[j]) +
//                                       (drdy * dNdr_[i] + dsdy * dNds_[i] + dtdy * dNdt_[i]) *
//                                       (drdy * dNdr_[j] + dsdy * dNds_[j] + dtdy * dNdt_[j]) +
//                                       (drdz * dNdr_[i] + dsdz * dNds_[i] + dtdz * dNdt_[i]) *
//                                       (drdz * dNdr_[j] + dsdz * dNds_[j] + dtdz * dNdt_[j])));

            mat_(j,i) = mat_(i,j);
        }
    }
    if (verbose) std::cout << "int ux2uy2uz2 " << *this << std::endl;
    return *this;
}

template < > DLLEXPORT void
ElementMatrix < double >::findWeightsAndPoints(const MeshEntity & ent,
                                               const RVector * &w,
                                               const PosVector * &x,
                                               int order){
    switch (ent.rtti()) {
        case MESH_EDGE_CELL_RTTI:
        case MESH_EDGE3_CELL_RTTI: {
            w = &IntegrationRules::instance().edgWeights(2);
            x = &IntegrationRules::instance().edgAbscissa(2);
        } break;
        case MESH_TRIANGLE_RTTI: {
            w = &IntegrationRules::instance().triWeights(1);
            x = &IntegrationRules::instance().triAbscissa(1);
        } break;
        case MESH_TRIANGLE6_RTTI: {
            w = &IntegrationRules::instance().triWeights(2);
            x = &IntegrationRules::instance().triAbscissa(2);
        } break;
        case MESH_QUADRANGLE_RTTI: {
            w = &IntegrationRules::instance().quaWeights(2);
            x = &IntegrationRules::instance().quaAbscissa(2);
        } break;
        case MESH_QUADRANGLE8_RTTI: {
            w = &IntegrationRules::instance().quaWeights(3);
            x = &IntegrationRules::instance().quaAbscissa(3);
        } break;
        case MESH_TETRAHEDRON_RTTI: {
            w = &IntegrationRules::instance().tetWeights(1);
            x = &IntegrationRules::instance().tetAbscissa(1);
        } break;
        case MESH_TETRAHEDRON10_RTTI: {
            w = &IntegrationRules::instance().tetWeights(2);
            x = &IntegrationRules::instance().tetAbscissa(2);
        } break;
        case MESH_HEXAHEDRON_RTTI: {
            w = &IntegrationRules::instance().hexWeights(2);
            x = &IntegrationRules::instance().hexAbscissa(2);
        } break;
        case MESH_HEXAHEDRON20_RTTI: {
            w = &IntegrationRules::instance().hexWeights(4);
            x = &IntegrationRules::instance().hexAbscissa(4);
        } break;
        case MESH_TRIPRISM_RTTI: {
            w = &IntegrationRules::instance().priWeights(2);
            x = &IntegrationRules::instance().priAbscissa(2);
        } break;
        case MESH_TRIPRISM15_RTTI: {
            w = &IntegrationRules::instance().priWeights(4);
            x = &IntegrationRules::instance().priAbscissa(4);
        } break;
        default:
            std::cerr << ent.rtti() << std::endl;
            THROW_TO_IMPL
            break;
    }
}

template < > DLLEXPORT
void ElementMatrix < double >::fillGradientBase(
                                    const MeshEntity & ent,
                                    const RVector & w,
                                    const PosVector & x,
                                    Index nC,
                                    bool voigtNotation){

    Index nRules = x.size();
    Index nDof = this->mat_.cols();
    Index nVerts = ent.nodeCount();
    if (_B.size() != nRules){
        _B.resize(nRules);
        for (Index i = 0; i < nRules; i ++ ){
            _B[i].resize(nC, nDof);
        }
    }
    if (dNdr_.rows() != nRules){
        if (ent.dim() > 0) dNdr_.resize(nRules, nVerts);
        // if (ent.dim() > 1) dNds_.resize(nRules, nVerts);
        if (ent.dim() > 2) dNdt_.resize(nRules, nVerts);

        for (Index i = 0; i < nRules; i ++){
            if (ent.dim() > 0) dNdr_[i] = ent.dNdL(x[i], 0);
            if (ent.dim() > 1) dNds_[i] = ent.dNdL(x[i], 1);
            if (ent.dim() > 2) dNdt_[i] = ent.dNdL(x[i], 2);
        }

        if (ent.dim() > 0) dNdx_.resize(nRules, nVerts);
        if (ent.dim() > 1) dNdy_.resize(nRules, nVerts);
        if (ent.dim() > 2) dNdz_.resize(nRules, nVerts);
    }
    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double drdz = ent.shape().drstdxyz(0, 2);
    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);
    double dsdz = ent.shape().drstdxyz(1, 2);
    double dtdx = ent.shape().drstdxyz(2, 0);
    double dtdy = ent.shape().drstdxyz(2, 1);
    double dtdz = ent.shape().drstdxyz(2, 2);

    for (Index i = 0; i < nRules; i ++){
        if (ent.dim() == 1){
            dNdx_[i].assign(dNdr_[i] * drdx);
        } else if (ent.dim() == 2){
            dNdx_[i].assign(dNdr_[i] * drdx + dNds_[i] * dsdx);
            dNdy_[i].assign(dNdr_[i] * drdy + dNds_[i] * dsdy);
        } else if (ent.dim() == 3){
            dNdx_[i].assign(dNdr_[i] * drdx + dNds_[i] * dsdx + dNdt_[i] * dtdx);
            dNdy_[i].assign(dNdr_[i] * drdy + dNds_[i] * dsdy + dNdt_[i] * dtdy);
            dNdz_[i].assign(dNdr_[i] * drdz + dNds_[i] * dsdz + dNdt_[i] * dtdz);
        }
    }

    double a = 1./std::sqrt(2.);
    if (voigtNotation){
        a = 1.0;
    }

    for (Index i = 0; i < nRules; i ++){
        // __MS(i, x[i])
        //dNdx_[i] = (dNdr_[i] * drdx);
        //dNdx_.setCol(i, (dNdr_[i] * drdx));

        if (this->_nDof == 0){
            // scalar field
            if (ent.dim() > 0){
                _B[i][0].setVal(dNdx_[i], 0, nVerts);
            }
            if (ent.dim() > 1){
                _B[i][1].setVal(dNdy_[i], 0, nVerts);
            }
            if (ent.dim() > 2){
                _B[i][2].setVal(dNdz_[i], 0, nVerts);
            }
        } else {
            // vector field
            if (ent.dim() == 1){
                _B[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts);
            } else if (ent.dim() == 2){
                _B[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts);
                _B[i][1].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts);

                if (nC > ent.dim()){ // elastic material
                    _B[i][2].setVal(dNdy_[i] * a, 0 * nVerts, 1 * nVerts);
                    _B[i][2].setVal(dNdx_[i] * a, 1 * nVerts, 2 * nVerts);
                }
            } else if (ent.dim() == 3){
                _B[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts);
                _B[i][1].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts);
                _B[i][2].setVal(dNdz_[i], 2 * nVerts, 3 * nVerts);

                if (nC > ent.dim()){ // elastic material
                    _B[i][3].setVal(dNdy_[i] * a, 0 * nVerts, 1 * nVerts);
                    _B[i][3].setVal(dNdx_[i] * a, 1 * nVerts, 2 * nVerts);

                    _B[i][4].setVal(dNdz_[i] * a, 1 * nVerts, 2 * nVerts);
                    _B[i][4].setVal(dNdy_[i] * a, 2 * nVerts, 3 * nVerts);

                    _B[i][5].setVal(dNdz_[i] * a, 0 * nVerts, 1 * nVerts);
                    _B[i][5].setVal(dNdx_[i] * a, 2 * nVerts, 3 * nVerts);
                }
            }
        }
    }
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU(const MeshEntity & ent,
                                                           const RVector & w,
                                                           const PosVector & x,
                                                           Index nC,
                                                           bool voigtNotation){
    this->fillIds(ent, nC); // also cleans
    this->fillGradientBase(ent, w, x, nC, voigtNotation);

    for (Index i = 0; i < w.size(); i ++ ){
        for (Index j = 0; j < nC; j ++) {
            // __MS(*this)
            // __MS(_B[i])
            // __MS(_B[i][j])
// #if USE_EIGEN3
// __MS("EIGENNEEDFIX")
// #else

            mat_[j * ent.nodeCount()] += _B[i][j] * w[i] * ent.size();
// #endif

        }
        // check performance if this works
        // iterator over weights in fill Gradient
        // this *= ent.size();
    }
    //mat_ *= ent.size();

    return * this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU(const Cell & cell,
                                 Index nC,
                                 bool voigtNotation){
    const RVector * w = 0;
    const PosVector * x = 0;
    this->findWeightsAndPoints(cell, w, x, 1);
    return this->gradU(cell, *w, *x, nC, voigtNotation);
}

template < > DLLEXPORT
RVector ElementMatrix < double >::stress(const MeshEntity & ent,
                                         const RSmallMatrix & C,
                                         const RVector & u, bool voigtNotation){
    const RVector * w = 0;
    const PosVector * x = 0;

    this->findWeightsAndPoints(ent, w, x, 1);
    this->fillIds(ent, C.size()); // also cleans
    this->fillGradientBase(ent, *w, *x,
                           max(C.size(), ent.dim()),
                           voigtNotation);

    RVector ret(C.rows(), 0.0);
    for (Index i = 0; i < w->size(); i ++ ){
        // C * B * u
        ret += C * (_B[i] * u[_ids]) * (*w)[i];
    }

    // performance optimization using _grad from fill GradientBase
    // _grad = sum(_B*w) // sum in GradientBase
    //ret =  C * (_grad * u[_ids])

    return ret;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU2(const MeshEntity & ent,
                                                            const Matrix< double > & C,
                                                            const RVector & w,
                                                            const PosVector & x,
                                                            bool voigtNotation){
    // __MS("inuse?") // somewhere in pgtest
    this->fillIds(ent, C.size()); // also cleans
    this->fillGradientBase(ent, w, x,
                           max(C.size(), ent.dim()),
                           voigtNotation);
    double beta = 0.0;
    if (C.size() == 1){

        for (Index i = 0; i < w.size(); i ++ ){
            // B.T * B
            if (i > 0) beta = 1.0;
            // #if USE_EIGEN3
            //     __MS("eigen3")
            // #else
            _B[i].transMult(_B[i], mat_, w[i] * ent.size() * C[0][0], beta);
            // #endif
        }
    } else {
        for (Index i = 0; i < w.size(); i ++ ){
            // B.T * C * B
            if (i > 0) beta = 1.0;
            // #if USE_EIGEN3
            //     __MS("eigen3")
            // #else
            matMultABA(_B[i], C, mat_, _abaTmp, w[i] * ent.size(), beta);
            // #endif
        }
        // check performance if this works
        // iterator over weights in fill Gradient
        // matMultABA(_B, C, mat_, _abaTmp, 1.0);
        // this *= ent.size();
    }
    // M += w[i] * B.T @ C @ B
    //mat_ *= ent.size();

    return * this;
}


template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::gradU2(const Cell & cell,
                                                    const Matrix< double > & C,
                                                    bool voigtNotation){
    // __MS("inuse?") // somewhere in pgtest

    const RVector * w = 0;
    const PosVector * x = 0;
    this->findWeightsAndPoints(cell, w, x, 1);
    return this->gradU2(cell, C, *w, *x, voigtNotation);
}


template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::u(const MeshEntity & ent){
    this->fillIds(ent);

    switch(ent.rtti()){
        case MESH_BOUNDARY_NODE_RTTI:
            mat_(0,0) = 1.0;
        break;
        case MESH_EDGE_CELL_RTTI:
        case MESH_EDGE_RTTI:
        case MESH_EDGE3_CELL_RTTI:
        case MESH_EDGE3_RTTI:
            return u(ent, IntegrationRules::instance().edgWeights(2),
                     IntegrationRules::instance().edgAbscissa(2), false); //ch
        case MESH_TRIANGLE_RTTI:
        case MESH_TRIANGLEFACE_RTTI:
            return u(ent, IntegrationRules::instance().triWeights(2),
                     IntegrationRules::instance().triAbscissa(2), false); //ch
        case MESH_TRIANGLE6_RTTI:
        case MESH_TRIANGLEFACE6_RTTI:
            return u(ent, IntegrationRules::instance().triWeights(2),
                     IntegrationRules::instance().triAbscissa(2), false); //ch
        case MESH_QUADRANGLE_RTTI:
        case MESH_QUADRANGLE8_RTTI:
            return u(ent, IntegrationRules::instance().quaWeights(2),
                     IntegrationRules::instance().quaAbscissa(2), false); //ch
        case MESH_QUADRANGLEFACE_RTTI:
        case MESH_QUADRANGLEFACE8_RTTI:
            return u(ent, IntegrationRules::instance().quaWeights(2),
                     IntegrationRules::instance().quaAbscissa(2), false); //ch
        case MESH_TETRAHEDRON_RTTI:
        case MESH_TETRAHEDRON10_RTTI:
            return u(ent, IntegrationRules::instance().tetWeights(2),
                     IntegrationRules::instance().tetAbscissa(2), false);
        case MESH_HEXAHEDRON_RTTI:
        case MESH_HEXAHEDRON20_RTTI:
            return u(ent, IntegrationRules::instance().hexWeights(2),
                     IntegrationRules::instance().hexAbscissa(2), false);
        case MESH_TRIPRISM_RTTI:
        case MESH_TRIPRISM15_RTTI:
            return u(ent, IntegrationRules::instance().priWeights(2),
                     IntegrationRules::instance().priAbscissa(2), false);

        default: std::cerr << WHERE_AM_I << " celltype not specified " << ent.rtti() << std::endl;
    }
    return *this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::u2(const MeshEntity & ent){
    this->fillIds(ent);

    switch(ent.rtti()){
    case MESH_BOUNDARY_NODE_RTTI:
        mat_(0,0) = 1.0;
    break;
    case MESH_EDGE_CELL_RTTI:
    case MESH_EDGE_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         mat_[0][0] = J / 3.0;
//         mat_[1][0] = J / 6.0;
//         mat_[0][1] = J / 6.0;
//         mat_[1][1] = J / 3.0;
//         std::cout << "2 " << *this << std::endl;
/*}
    break;*/
        u2(ent, IntegrationRules::instance().edgWeights(2), IntegrationRules::instance().edgAbscissa(2), false);
        break;
    case MESH_EDGE3_CELL_RTTI:
    case MESH_EDGE3_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         mat_[0][0] =   J / 30.0 * 4.0;
//         mat_[1][0] = - J / 30.0;
//         mat_[2][0] =   J / 15.0;
//         mat_[0][1] = - J / 30.0;
//         mat_[1][1] =   J / 30.0 * 4.0;
//         mat_[2][1] =   J / 15.0;
//         mat_[0][2] =   J / 15.0;
//         mat_[1][2] =   J / 15.0;
//         mat_[2][2] =   J / 30.0 * 16.0;
//          std::cout << "2 " << *this << std::endl;
    //} break;
        u2(ent, IntegrationRules::instance().edgWeights(3), IntegrationRules::instance().edgAbscissa(3), false);
        break;
    case MESH_TRIANGLE_RTTI:
    case MESH_TRIANGLEFACE_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         double Jl = J / 24.0;
//         mat_[0][0] = Jl * 2.0;
//         mat_[1][0] = Jl;
//         mat_[2][0] = Jl;
//         mat_[0][1] = Jl;
//         mat_[1][1] = Jl * 2.0;
//         mat_[2][1] = Jl;
//         mat_[0][2] = Jl;
//         mat_[1][2] = Jl;
//         mat_[2][2] = Jl * 2.0;
//         std::cout << "2 " << *this << std::endl;
//}break;
        u2(ent, IntegrationRules::instance().triWeights(2), IntegrationRules::instance().triAbscissa(2), false);
        break;
    case MESH_QUADRANGLE_RTTI:
    case MESH_QUADRANGLEFACE_RTTI:
        u2(ent, IntegrationRules::instance().quaWeights(2), IntegrationRules::instance().quaAbscissa(2), false);
        break;
    case MESH_QUADRANGLE8_RTTI:
    case MESH_QUADRANGLEFACE8_RTTI:
        u2(ent, IntegrationRules::instance().quaWeights(3), IntegrationRules::instance().quaAbscissa(3), false);
        break;
    case MESH_TRIANGLE6_RTTI:
    case MESH_TRIANGLEFACE6_RTTI:
    //{
//         double J = ent.shape().jacobianDeterminant();
//         if (J < 0) {
//             std::cerr << WHERE_AM_I << "JacobianDeterminant < 0 (" << J << ")" << std::endl;
//             std::cerr << ent.shape() << std::endl;
//             J = std::fabs(J);
//         }
//         for (uint i = 0; i < nVerts; i++){
//             for (uint j = 0; j < nVerts; j++){
//                 mat_(i,j) = J * (*Tri6_u2)[i][j];//compound[i][j];
//             }
//         }
//         std::cout << "2 " << *this << std::endl;
    //} break;
        return u2(ent, IntegrationRules::instance().triWeights(4), IntegrationRules::instance().triAbscissa(4), false); //ch
    case MESH_TETRAHEDRON_RTTI:
        return u2(ent, IntegrationRules::instance().tetWeights(2), IntegrationRules::instance().tetAbscissa(2), false);
    case MESH_TETRAHEDRON10_RTTI:
        return u2(ent, IntegrationRules::instance().tetWeights(4), IntegrationRules::instance().tetAbscissa(4), false);
    case MESH_HEXAHEDRON_RTTI:
        return u2(ent, IntegrationRules::instance().hexWeights(2), IntegrationRules::instance().hexAbscissa(2), false);
    case MESH_HEXAHEDRON20_RTTI:
        return u2(ent, IntegrationRules::instance().hexWeights(4), IntegrationRules::instance().hexAbscissa(4), false);
    case MESH_TRIPRISM_RTTI:
        return u2(ent, IntegrationRules::instance().priWeights(2), IntegrationRules::instance().priAbscissa(2), false);
    case MESH_TRIPRISM15_RTTI:
        return u2(ent, IntegrationRules::instance().priWeights(4), IntegrationRules::instance().priAbscissa(4), false);
    default:
      std::cerr << ent.rtti() << std::endl;
      THROW_TO_IMPL
    }

    return *this;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::ux2uy2uz2(const Cell & cell, bool useCache){

    fillIds(cell);

//     double J = cell.jacobianDeterminant();
//     if (J <= 0) std::cerr << WHERE_AM_I << " JacobianDeterminant < 0 (" << J << ") " << cell << std::endl;
//      std::cout << J << std::endl;

    switch (cell.rtti()) {
    case MESH_EDGE_CELL_RTTI:
    case MESH_EDGE3_CELL_RTTI:
        ux2(cell, IntegrationRules::instance().edgWeights(2),
            IntegrationRules::instance().edgAbscissa(2), false); break;
    case MESH_TRIANGLE_RTTI: {
        double J = cell.size()*2.;
    ////////////////////////////////////////////////////////////////////
//         double dN1dx = cell.shape().deriveCoordinates(0, 0);
//         double dN2dx = cell.shape().deriveCoordinates(1, 0);
//         double dN3dx = cell.shape().deriveCoordinates(2, 0);
//         double dN1dy = cell.shape().deriveCoordinates(0, 1);
//         double dN2dy = cell.shape().deriveCoordinates(1, 1);
//         double dN3dy = cell.shape().deriveCoordinates(2, 1);
//         mat_[0][0] = J / 2.0 * (dN1dx * dN1dx + dN1dy * dN1dy);
//         mat_[1][0] = J / 2.0 * (dN2dx * dN1dx + dN2dy * dN1dy);
//         mat_[2][0] = J / 2.0 * (dN3dx * dN1dx + dN3dy * dN1dy);
//         mat_[1][1] = J / 2.0 * (dN2dx * dN2dx + dN2dy * dN2dy);
//         mat_[2][1] = J / 2.0 * (dN3dx * dN2dx + dN3dy * dN2dy);
//         mat_[2][2] = J / 2.0 * (dN3dx * dN3dx + dN3dy * dN3dy);
//         mat_[0][1] = mat_[1][0];
//         mat_[0][2] = mat_[2][0];
//         mat_[1][2] = mat_[2][1];
////////////////////////////////////////////////////////////////////
        // this is much faster than numerical integration

//         double x21 = cell.node(1).x()-cell.node(0).x();
//         double x31 = cell.node(2).x()-cell.node(0).x();
//         double y21 = cell.node(1).y()-cell.node(0).y();
//         double y31 = cell.node(2).y()-cell.node(0).y();
//
//         double a =   ((x31) * (x31) + (y31) * (y31)) / J;
//         double b = - ((x31) * (x21) + (y31) * (y21)) / J;
//         double c =   ((x21) * (x21) + (y21) * (y21)) / J;

        double x1 = cell.node(0).x();
        double x2 = cell.node(1).x();
        double x3 = cell.node(2).x();
        double y1 = cell.node(0).y();
        double y2 = cell.node(1).y();
        double y3 = cell.node(2).y();

        double a =   ((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)) / J;
        double b = - ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / J;
        double c =   ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) / J;

        mat_(0,0) = a *  0.5 + b        + c *  0.5 ;
        mat_(1,0) = a * -0.5 + b * -0.5            ;
        mat_(2,0) =            b * -0.5 + c * -0.5 ;
        mat_(1,1) = a *  0.5                       ;
        mat_(2,1) =            b *  0.5            ;
        mat_(2,2) =                       c *  0.5 ;

        mat_(0,1) = mat_(1,0);
        mat_(0,2) = mat_(2,0);
        mat_(1,2) = mat_(2,1);

        //std::cout << "2" << *this << std::endl;*/
        //ux2uy2(cell, IntegrationRules::instance().triWeights(1), IntegrationRules::instance().triAbscissa(1), false);
    } break;
    case MESH_TRIANGLE6_RTTI: {
///////////////////////////////////////////////////////////////////////////////////
//         double x1 = cell.node(0).x();
//         double x2 = cell.node(1).x();
//         double x3 = cell.node(2).x();
//         double y1 = cell.node(0).y();
//         double y2 = cell.node(1).y();
//         double y3 = cell.node(2).y();
//
//         double b1 = y2-y3;
//         double b2 = y3-y1;
//         double b3 = y1-y2;
//         double c1 = x3-x2;
//         double c2 = x1-x3;
//         double c3 = x2-x1;
//         double a = (c3*b2-c2*b3)/2.0;
//         double b1sqr=b1*b1;
//         double b2sqr=b2*b2;
//         double b3sqr=b3*b3;
//         double c1sqr=c1*c1;
//         double c2sqr=c2*c2;
//         double c3sqr=c3*c3;
//         double b1b2=b1*b2;
//         double b1b3=b1*b3;
//         double b2b3=b2*b3;
//         double c1c2=c1*c2;
//         double c1c3=c1*c3;
//         double c2c3=c2*c3;
//
//         mat_[0][0]=(b1sqr+c1sqr) / (4.0*a);
//         mat_[0][1]=(-b1b2-c1c2)/(12.0*a);
//         mat_[0][2]=(-b1b3-c1c3)/(12.0*a);
//         mat_[0][3]=(b1b2+c1c2)/(3.0*a);
//         mat_[0][4]=0.0;
//         mat_[0][5]=(b1b3+c1c3)/(3.0*a);
//         mat_[1][1]=(b2sqr+c2sqr)/(4.0*a);
//         mat_[1][2]=(-b2b3-c2c3)/(12.0*a);
//         mat_[1][3]=(b1b2+c1c2)/(3.0*a);
//         mat_[1][4]=(b2b3+c2c3)/(3.0*a);
//         mat_[1][5]=0.0;
//         mat_[2][2]=(b3sqr+c3sqr)/(4.0*a);
//         mat_[2][3]=0.0;
//         mat_[2][4]=(b2b3+c2c3)/(3.0*a);
//         mat_[2][5]=(b1b3+c1c3)/(3.0*a);
//         mat_[3][3]=2.0*((b1sqr+b1b2+b2sqr)+(c1sqr+c1c2+c2sqr))/(3.0*a);
//         mat_[3][4]=((b1b2+b2sqr+2.0*b1b3+b2b3)+(c1c2+c2sqr+2.0*c1c3+c2c3))/(3.0*a);
//         mat_[3][5]=((b1sqr+b1b3+b1b2+2.0*b2b3)+(c1sqr+c1c3+c1c2+2.0*c2c3))/(3.0*a);
//         mat_[4][4]=2.0*((b2sqr+b2b3+b3sqr)+(c2sqr+c2c3+c3sqr))/(3.0*a);
//         mat_[4][5]=((2.0*b1b2+b2b3+b1b3+b3sqr)+(2.0*c1c2+c2c3+c1c3+c3sqr))/(3.0*a);
//         mat_[5][5]=2.0*((b1sqr+b1b3+b3sqr)+(c1sqr+c1c3+c3sqr))/(3.0*a);
//
//         for (int i = 1, imax = 6; i < imax; i ++){
//              for (int j = 0, jmax = i; j < jmax; j ++){
//                 mat_(i,j)=mat_(j,i);
//              }
//         }

//working version
//         double x1 = cell.node(0).x();
//         double x2 = cell.node(1).x();
//         double x3 = cell.node(2).x();
//         double y1 = cell.node(0).y();
//         double y2 = cell.node(1).y();
//         double y3 = cell.node(2).y();
//
//         double a = ((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)) / J / 6.0;
//         double b = - ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / J / 6.0;
//         double c = ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) / J / 6.0;
//         for (int i = 0, imax = 6; i < imax; i ++){
//             for (int j = 0, jmax = 6; j < jmax; j ++){
//                 mat_(i,j) = Triangle6_S1[i][j] * a +
//                                  Triangle6_S2[i][j] * b +
//                                  Triangle6_S3[i][j] * c;
//             }
//         }
//         std::cout << "2" << *this << std::endl;
        ux2uy2(cell,
               IntegrationRules::instance().triWeights(2),
               IntegrationRules::instance().triAbscissa(2), false); //ch
    } break;
    case MESH_QUADRANGLE_RTTI:
        ux2uy2(cell,
               IntegrationRules::instance().quaWeights(2),
               IntegrationRules::instance().quaAbscissa(2), false); break;
    case MESH_QUADRANGLE8_RTTI:
        ux2uy2(cell,
               IntegrationRules::instance().quaWeights(3),
               IntegrationRules::instance().quaAbscissa(3), false); break;
    case MESH_TETRAHEDRON_RTTI:
    //{
//         double x_xi = cell.shape().partDerivationRealToUnity(0, 1);
//         double x_eta = cell.shape().partDerivationRealToUnity(0, 2);
//         double x_zeta = cell.shape().partDerivationRealToUnity(0, 3);
//
//         double y_xi = cell.shape().partDerivationRealToUnity(1, 1);
//         double y_eta = cell.shape().partDerivationRealToUnity(1, 2);
//         double y_zeta = cell.shape().partDerivationRealToUnity(1, 3);
//
//         double z_xi = cell.shape().partDerivationRealToUnity(2, 1);
//         double z_eta = cell.shape().partDerivationRealToUnity(2, 2);
//         double z_zeta = cell.shape().partDerivationRealToUnity(2, 3);
//
//         double xi_x =  1.0 / J * det(y_eta, z_eta, y_zeta, z_zeta);
//         double xi_y = -1.0 / J * det(x_eta, z_eta, x_zeta, z_zeta);
//         double xi_z =  1.0 / J * det(x_eta, y_eta, x_zeta, y_zeta);
//
//         double eta_x = -1.0 / J * det(y_xi, z_xi, y_zeta, z_zeta);
//         double eta_y =  1.0 / J * det(x_xi, z_xi, x_zeta, z_zeta);
//         double eta_z = -1.0 / J * det(x_xi, y_xi, x_zeta, y_zeta);
//
//         double zeta_x =  1.0 / J * det(y_xi, z_xi, y_eta, z_eta);
//         double zeta_y = -1.0 / J * det(x_xi, z_xi, x_eta, z_eta);
//         double zeta_z =  1.0 / J * det(x_xi, y_xi, x_eta, y_eta);
//
//         double a = J / 6.0 * (xi_x * xi_x + xi_y * xi_y + xi_z * xi_z);
//         double b = J / 6.0 * (eta_x * eta_x + eta_y * eta_y + eta_z * eta_z);
//         double c = J / 6.0 * (zeta_x * zeta_x + zeta_y * zeta_y + zeta_z * zeta_z);
//
//         double d = J / 6.0 * (xi_x * eta_x + xi_y * eta_y + xi_z * eta_z);
//         double e = J / 6.0 * (xi_x * zeta_x + xi_y * zeta_y + xi_z * zeta_z);
//         double f = J / 6.0 * (eta_x * zeta_x + eta_y * zeta_y + eta_z * zeta_z);
//
//     double u_xi2[4][4] = {
//       {   a,  -a, 0.0, 0.0 },
//       {  -a,   a, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_eta2[4][4] = {
//       {   b, 0.0,  -b, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       {  -b, 0.0,   b, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_zeta2[4][4] = {
//       {   c, 0.0, 0.0,  -c },
//       { 0.0, 0.0, 0.0, 0.0 },
//       { 0.0, 0.0, 0.0, 0.0 },
//       {  -c, 0.0, 0.0,   c }
//     };
//     double u_xi__u_eta[4][4] = {
//       { 2.0 * d,  -d,  -d, 0.0 },
//       {      -d, 0.0,   d, 0.0 },
//       {      -d,   d, 0.0, 0.0 },
//       {     0.0, 0.0, 0.0, 0.0 }
//     };
//     double u_xi__u_zeta[4][4] = {
//       { 2.0 * e,  -e, 0.0,  -e },
//       {      -e, 0.0, 0.0,   e },
//       {     0.0, 0.0, 0.0, 0.0 },
//       {      -e,   e, 0.0, 0.0 }
//     };
//     double u_eta__u_zeta[4][4] = {
//       { 2.0 * f, 0.0,  -f,  -f },
//       {     0.0, 0.0, 0.0, 0.0 },
//       {      -f, 0.0, 0.0,   f },
//       {      -f, 0.0,   f, 0.0 }
//     };
//
//     for (uint i = 0; i < dim; i++){
//       for (uint j = 0; j < dim; j++){
//  mat_(i,j) = u_xi2[i][j] + u_eta2[i][j] +  u_zeta2[i][j] +
//    u_xi__u_eta[i][j] + u_xi__u_zeta[i][j] + u_eta__u_zeta[i][j];
//       }
//     }
//     std::cout << "0 " << *this << std::endl;
//} break;
        ux2uy2uz2(cell,
                  IntegrationRules::instance().tetWeights(1),
                  IntegrationRules::instance().tetAbscissa(1), false); //ch
        break;
    case MESH_TETRAHEDRON10_RTTI:
//     {
// ////////////////////////////////////////////////////////////////////
//     double x_xi = cell.shape().partDerivationRealToUnity(0, 1);
//     double x_eta = cell.shape().partDerivationRealToUnity(0, 2);
//     double x_zeta = cell.shape().partDerivationRealToUnity(0, 3);
//
//     double y_xi = cell.shape().partDerivationRealToUnity(1, 1);
//     double y_eta = cell.shape().partDerivationRealToUnity(1, 2);
//     double y_zeta = cell.shape().partDerivationRealToUnity(1, 3);
//
//     double z_xi = cell.shape().partDerivationRealToUnity(2, 1);
//     double z_eta = cell.shape().partDerivationRealToUnity(2, 2);
//     double z_zeta = cell.shape().partDerivationRealToUnity(2, 3);
//
//     double xi_x =  1.0 / J * det(y_eta, z_eta, y_zeta, z_zeta);
//     double xi_y = -1.0 / J * det(x_eta, z_eta, x_zeta, z_zeta);
//     double xi_z =  1.0 / J * det(x_eta, y_eta, x_zeta, y_zeta);
//
//     double eta_x = -1.0 / J * det(y_xi, z_xi, y_zeta, z_zeta);
//     double eta_y =  1.0 / J * det(x_xi, z_xi, x_zeta, z_zeta);
//     double eta_z = -1.0 / J * det(x_xi, y_xi, x_zeta, y_zeta);
//
//     double zeta_x =  1.0 / J * det(y_xi, z_xi, y_eta, z_eta);
//     double zeta_y = -1.0 / J * det(x_xi, z_xi, x_eta, z_eta);
//     double zeta_z =  1.0 / J * det(x_xi, y_xi, x_eta, y_eta);
//
//     double a = J / 6.0 * (xi_x * xi_x + xi_y * xi_y + xi_z * xi_z);
//     double b = J / 6.0 * (eta_x * eta_x + eta_y * eta_y + eta_z * eta_z);
//     double c = J / 6.0 * (zeta_x * zeta_x + zeta_y * zeta_y + zeta_z * zeta_z);
//
//     double d = J / 6.0 * (xi_x * eta_x + xi_y * eta_y + xi_z * eta_z);
//     double e = J / 6.0 * (xi_x * zeta_x + xi_y * zeta_y + xi_z * zeta_z);
//     double f = J / 6.0 * (eta_x * zeta_x + eta_y * zeta_y + eta_z * zeta_z);
//
//     RSTLMatrix compound(10, 10);
//     compound = a * (*Tet10_u_xi2) + b * (*Tet10_u_eta2) + c * (*Tet10_u_zeta2)
//         + (2.0*d) * (*Tet10_u_xi_u_eta)
//         + (2.0*e) * (*Tet10_u_xi_u_zeta)
//         + (2.0*f) * (*Tet10_u_eta_u_zeta);
//     for (uint i = 0; i < dim; i++){
//       for (uint j = 0; j < dim; j++){
//         //** * 6.0 because a b c d e f / 6
//         mat_(i,j) = compound[i][j] * 6.0;
//       }
//     }
// //    std::cout << " 2 " << *this << std::endl;
//
//     break;
//   }
        ux2uy2uz2(cell,
                  IntegrationRules::instance().tetWeights(2),
                  IntegrationRules::instance().tetAbscissa(2), false);   break;
    case MESH_HEXAHEDRON_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().hexWeights(2),
                  IntegrationRules::instance().hexAbscissa(2), false);   break;
    case MESH_HEXAHEDRON20_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().hexWeights(4),
                  IntegrationRules::instance().hexAbscissa(4), false);   break;
    case MESH_TRIPRISM_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().priWeights(2),
                  IntegrationRules::instance().priAbscissa(2), false);   break;
    case MESH_TRIPRISM15_RTTI:
        ux2uy2uz2(cell,
                  IntegrationRules::instance().priWeights(4),
                  IntegrationRules::instance().priAbscissa(4), false);   break;

    default:
        std::cerr << cell.rtti() << std::endl;
        THROW_TO_IMPL
        break;
    }

    return *this;
}

template < >
void ElementMatrix < double >::fillIds(const MeshEntity & ent, Index nC){
    Index nDims = 1;
    Index nNodes = ent.nodeCount();

    // __MS(nC)
    if (this->_nDof > 0){
        nDims = ent.dim();
        if (nC > ent.dim()){
            //THROW_TO_IMPL
        }

        if (size() != nNodes * nDims) resize(nNodes * nDims);

        for (Index dim = 0; dim < nDims; dim++){
            for (Index i = 0; i < nNodes; i ++) {
                _ids[dim * nNodes + i]  = dim * this->_nDof + ent.node(i).id();
                _idsC[dim * nNodes + i] = _ids[i + dim * nNodes];
                _idsR[dim * nNodes + i] = _ids[i + dim * nNodes];
            }
        }
    } else {
        nDims = 1; //nC;
        this->resize(nNodes*nC, nNodes);
        // mat_.resize(nNodes*nC, nNodes);
        // _ids.resize(nNodes*nC);

        for (Index dim = 0; dim < nDims; dim++){
            for (Index i = 0; i < nNodes; i ++) {
                _ids[dim * nNodes + i] = ent.node(i).id();
                _idsC[dim * nNodes + i] = _ids[dim * nNodes + i];
                _idsR[i] = _ids[dim * nNodes + i];
            }
        }
    }

    // __MS(_ids)
    // __MS(_idsR)
    // __MS(_idsC)

    *this *= 0.0;
}

//** old interface constructor
template < > DLLEXPORT
ElementMatrix < double >::ElementMatrix(Index dof) {
    // __MS("EM()", this);
    this->init(0, 0, 0);
    this->_nDof = dof;
    this->_newStyle = false;

    // this->_nCoeff = 0;
    // this->_dofPerCoeff = 0;
    // this->_dofOffset = 0;
    // _ent = 0;
    // _w = 0;
    // _x = 0;
}

//** new interface starts here
template < > DLLEXPORT
ElementMatrix < double >::ElementMatrix(Index nCoeff, Index dofPerCoeff,
                                        Index dofOffset){
    // __MS("EM(a,b,c)", this);
    this->init(nCoeff, dofPerCoeff, dofOffset);
}

template < > DLLEXPORT
ElementMatrix < double >::ElementMatrix(const ElementMatrix < double > & E){
    // __MS("EM(a,b,c)", this);
    this->copyFrom(E, true);
}

template < > DLLEXPORT
ElementMatrix < double >::ElementMatrix(const ElementMatrix < double > & E,
                                        bool withMat){
    this->copyFrom(E, withMat);
}

template < > DLLEXPORT
void ElementMatrix < double >::init(Index nCoeff, Index dofPerCoeff,
                                    Index dofOffset){
    // __MS(this, nCoeff, dofPerCoeff, dofOffset)
    if (nCoeff > 1 && dofPerCoeff == 0){
        __MS(nCoeff, dofPerCoeff, dofOffset)
        log(Error, "Number of coefficents > 1 but no dofPerCoefficent given");
    }
    this->_newStyle = true;
    this->_order = 0;
    this->_nCoeff = nCoeff;
    this->_dofPerCoeff = dofPerCoeff;
    this->_dofOffset = dofOffset;

    _ent = 0;
    _w = 0;
    _x = 0;

    _div = false;
    _integrated = false;
    _elastic = false;

    _valid = false;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::add(const ElementMatrix < double > & B,
                              Index dim, double b) {
    this->_integrated = false;

    if (this->rows() == B.rows() && this->cols() == B.cols()){
        return (*this) += B;
    }

    if (dim == 1){
        this->resize(this->rows(), dim);

        for (Index r=0; r < this->_matX.size(); r++){
            auto & Ar = this->_matX[r];
            const auto & Br = B.matX()[r];
            // print("A:", Ar);
            // print("B:", Br);
// #if USE_EIGEN3
//             // print("As:", Ar.colwise().sum());
//             // print("Bs:", Br.colwise().sum());

//             Ar = Ar.colwise().sum().eval() + Br.colwise().sum();
// #else
            for (Index i=0; i < Br.size(); i ++){
                Ar[0] += Br[i];
            }
            for (Index i=1; i < Ar.size(); i ++){
                Ar[0] += Ar[i];
            }
            Ar.resize(1, Ar.cols());

// #endif
            // print("A2:", Ar);
        }
    } else {
        this->resize(this->rows(), max(this->cols(), B.cols()));

        for (Index r=0; r < this->_matX.size(); r++){
            auto & Ar = this->_matX[r];
            const auto & Br = B.matX()[r];
// #if USE_EIGEN3

//             if (Ar.rows() == 1 && Ar.cols() == Br.cols()){
//                 auto Ar0 = Ar.row(0).eval();

//                 Ar.resize(Br.rows(), Ar.cols());

//                 for (Index i = 0; i < Br.rows(); i ++){
//                     Ar.row(i) = Ar0 + Br.row(i);
//                 }
//             } else if (Br.rows() == 1 && Ar.cols() == Br.cols()){
//                 for (Index i = 0; i < Ar.rows(); i ++){
//                     Ar.row(i) += Br.row(0);
//                 }
//             } else {
//                 THROW_TO_IMPL
//             }
// #else
            Ar += Br;
// #endif
            // print("A2:", Ar);
        }

    }
    return *this;
}


template < > DLLEXPORT
void ElementMatrix < double >::copyFrom(const ElementMatrix < double > & E,
                                        bool withMat){
// __M
    this->_newStyle = true;
    this->_order = E.order(); //quadrature order
    this->_nCoeff = E.nCoeff(); //number of cofficients e.g, dim
    this->_dofPerCoeff = E.dofPerCoeff();
    this->_dofOffset = E.dofOffset();

    this->_ent = E.entity();
    this->_w = E.w();
    this->_x = E.x();
    this->_matX = E.matX();

    this->_idsC = E.colIDs();
    this->_idsR = E.rowIDs();
    this->_div = E.isDiv();
    this->_elastic = E.elastic();

    if (withMat == true) {
        this->_integrated = E.isIntegrated();
        this->mat_ = E.mat();
    } else {
        this->_integrated = false;
        this->mat_.resize(E.mat().rows(), E.mat().cols());
    }
    this->_valid = E.valid();
}


template < >
void ElementMatrix < double >::resize(Index rows, Index cols, bool setIds) {
    if (cols == 0) cols = rows;

    _idsR.resize(rows, 0);
    _idsC.resize(cols, 0);
    _ids.resize(rows);
    mat_.resize(rows, cols);

    if (this->_ent && setIds){
        Index nVerts = this->_ent->nodeCount();

        for (Index i = 0; i < this->_nCoeff; i++){
            this->_idsR.setVal(this->_ent->ids() + i * this->_dofPerCoeff + this->_dofOffset,
                            i * nVerts, (i+1) * nVerts);
        }
    }
}


template < > DLLEXPORT
void ElementMatrix < double >::integrate() const {
    // __MS(_newStyle, this->_integrated, this->_valid, this)

    if (_newStyle && !this->_integrated && this->_valid){
        // __M
        const RVector &w = *this->_w;
        Index nRules(w.size());

        this->mat_.setZero();

        for (Index i = 0; i < nRules; i ++){
            // improve either with matrix expressions or shift w as scale
            // __MS(_matX[i] * (w[i] * this->_ent->size()))

            this->mat_.transAdd((_matX[i] * (w[i] * this->_ent->size())));

            // MAT_TRANS_ADD(this->mat_, (_matX[i] * (w[i] * this->_ent->size())))

            // __MS(this->mat_)
        }

        this->_integrated = true;
    }
}


template < > DLLEXPORT
void ElementMatrix < double >::fillEntityAndOrder_(const MeshEntity & ent, Index order){
    this->_order = order;
    this->_ent = &ent;
    this->_integrated = false;

    this->_x = &IntegrationRules::instance().abscissa(ent.shape(), this->_order);
    this->_w = &IntegrationRules::instance().weights(ent.shape(), this->_order);
}


template < > DLLEXPORT
void ElementMatrix < double >::resizeMatX_U_(){
    Index nRules(_x->size());
    Index nVerts(_ent->nodeCount());
    Index nCoeff(_nCoeff); //components
    Index nCols(nCoeff);

    if (nCols == 0){
        log(Critical, "ElementMatrix need to be initialized");
    }
    this->resize(nVerts*nCoeff, nCols);

    _matX.resize(nRules);

    for (Index i = 0; i < nRules; i ++ ){
        // transpose might be better?? check
        // fill per row is cheaper
        _matX[i].resize(nCoeff, nVerts*nCoeff);
        _matX[i].setZero();
    }
}


template < > DLLEXPORT
void ElementMatrix < double >::fillMatX_U_(bool sum){
    Index nRules(_x->size());
    Index nVerts(_ent->nodeCount());
    Index nCoeff(this->_nCoeff);

    const PosVector &x = *this->_x;

    if (0 && nCoeff == 1){
        for (Index i = 0; i < nRules; i ++ ){
            _matX[i][0].setVal(_ent->N(x[i]), 0, nVerts);
        }
    } else {

        RSmallMatrix N(nRules, nVerts);

        for (Index i = 0; i < nRules; i ++ ){
            N[i] = _ent->N(x[i]);
            for (Index n = 0; n < nCoeff; n ++ ){
                _matX[i][n].setVal(N[i], n*nVerts, (n+1)*nVerts);
            }
        }
    }

    this->setValid(true);
    if (sum) this->integrate();
}


template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::pot(
                        const MeshEntity & ent, Index order, bool sum){

    if (this->_valid && this->order() == order && this->_ent == &ent){
        return *this;
    }

    this->fillEntityAndOrder_(ent, order);
    this->resizeMatX_U_();
    this->fillMatX_U_(sum);

    return *this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::pot(
                        const MeshEntity & ent, Index order, bool sum,
                        Index nCoeff, Index dofPerCoeff, Index dofOffset){

    // __M
    if (disableCacheForDBG() || !this->valid() ||
        this->order() != order ||
        this->_ent != &ent ||
        this->_nCoeff != nCoeff){

        this->init(nCoeff, dofPerCoeff, dofOffset);
        this->pot(ent, order, sum);
    }
    if (sum == true) this->integrate();

    // __MS(this->mat_)
    return *this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::grad(
                        const MeshEntity & ent, Index order,
                        bool elastic, bool sum, bool div, bool kelvin){

    if (this->valid() && \
        this->order() == order && \
        this->_ent == & ent && \
        this->elastic() == elastic){
        return *this;
    }

    this->_order = order;
    this->_ent = &ent;
    this->_div = div;
    this->_integrated = false;
    this->_elastic = elastic;

    //this->findWeightsAndPoints(ent, this->_w, this->_x, this->_order);

    this->_w = &IntegrationRules::instance().weights(ent.shape(), this->_order);
    this->_x = &IntegrationRules::instance().abscissa(ent.shape(), this->_order);

    const PosVector &x = *this->_x;
    // const RVector &w = *this->_w;

    Index nRules(_x->size());
    Index nVerts(_ent->nodeCount());
    Index nCoeff(this->_nCoeff); //components
    Index nCols(ent.dim() * nCoeff);

    if (nCols == 0){
        log(Critical, "ElementMatrix need to be initialized");
    }

    double a = 1.0;
    if (elastic == true){
        //** special case for constitutive matrix
        nCols = ent.dim();
        if (kelvin){
            a = 1./std::sqrt(2.);
        }

        if (ent.dim() == 2){
            nCols += 1;
        } else if (ent.dim() == 3){
            nCols += 3;
        }
    }

    this->resize(nVerts*nCoeff, nCols);
    // this->_idsR.resize(nVerts*nCoeff, 0);
    // this->_idsC.resize(nCols, 0);

    // for (Index i = 0; i < nCoeff; i++){
    //     this->_idsR.setVal(ent.ids() + i * _dofPerCoeff + _dofOffset,
    //                        i * nVerts, (i+1) * nVerts);
    // }

    // __MS(ent.dim(), nCols, elastic)
    // matrices per quadrature point
    //_matX

    _matX.resize(nRules);
    for (Index i = 0; i < nRules; i ++ ){
        //** _matX stored transposed
        // transpose might be better?? check
        // fill per row is cheaper
        _matX[i].resize(nCols, nVerts*nCoeff);
        _matX[i].setZero();
    }

    // __MS(this->rows(), this->cols())
    // __MS(_matX[0].rows(), _matX[0].cols())

    if (dNdr_.rows() != nRules){
        if (ent.dim() > 0) dNdr_.resize(nRules, nVerts);
        if (ent.dim() > 1) dNds_.resize(nRules, nVerts);
        if (ent.dim() > 2) dNdt_.resize(nRules, nVerts);

        for (Index i = 0; i < nRules; i ++){
            if (ent.dim() > 0) dNdr_[i] = ent.dNdL(x[i], 0);
            if (ent.dim() > 1) dNds_[i] = ent.dNdL(x[i], 1);
            if (ent.dim() > 2) dNdt_[i] = ent.dNdL(x[i], 2);
        }

        if (ent.dim() > 0) dNdx_.resize(nRules, nVerts);
        if (ent.dim() > 1) dNdy_.resize(nRules, nVerts);
        if (ent.dim() > 2) dNdz_.resize(nRules, nVerts);
    }
    double drdx = ent.shape().drstdxyz(0, 0);
    double drdy = ent.shape().drstdxyz(0, 1);
    double drdz = ent.shape().drstdxyz(0, 2);
    double dsdx = ent.shape().drstdxyz(1, 0);
    double dsdy = ent.shape().drstdxyz(1, 1);
    double dsdz = ent.shape().drstdxyz(1, 2);
    double dtdx = ent.shape().drstdxyz(2, 0);
    double dtdy = ent.shape().drstdxyz(2, 1);
    double dtdz = ent.shape().drstdxyz(2, 2);

    for (Index i = 0; i < nRules; i ++){
        if (ent.dim() == 1){
            dNdx_[i].assign(dNdr_[i] * drdx);
        } else if (ent.dim() == 2){
            dNdx_[i].assign(dNdr_[i] * drdx + dNds_[i] * dsdx);
            dNdy_[i].assign(dNdr_[i] * drdy + dNds_[i] * dsdy);
        } else if (ent.dim() == 3){
            dNdx_[i].assign(dNdr_[i] * drdx + dNds_[i] * dsdx + dNdt_[i] * dtdx);
            dNdy_[i].assign(dNdr_[i] * drdy + dNds_[i] * dsdy + dNdt_[i] * dtdy);
            dNdz_[i].assign(dNdr_[i] * drdz + dNds_[i] * dsdz + dNdt_[i] * dtdz);
        }
    }
    for (Index i = 0; i < nRules; i ++){

        if (nCoeff == 1){
            //** scalar field
            if (ent.dim() > 0){
                _matX[i][0].setVal(dNdx_[i], 0, nVerts); //dNx/dx
// SET_MAT_ROW_SLICE(_matX[i],0,dNdx_[i], 1, 0, nVerts);
            }
            if (ent.dim() > 1){
                _matX[i][1].setVal(dNdy_[i], 0, nVerts); //dNy/dy
// SET_MAT_ROW_SLICE(_matX[i],1,dNdy_[i], 1, 0, nVerts);
            }
            if (ent.dim() > 2){
                _matX[i][2].setVal(dNdz_[i], 0, nVerts); //dNz/dz
// SET_MAT_ROW_SLICE(_matX[i],2,dNdz_[i], 1, 0, nVerts);
            }
        } else {
            //** vector field
            if (ent.dim() == 1){
                _matX[i][0].setVal(dNdx_[i], 0, nVerts); //dNx/dx
// SET_MAT_ROW_SLICE(_matX[i],0,dNdx_[i], 1, 0 * nVerts, 1 * nVerts);
            } else if (ent.dim() == 2){
                if (elastic == true){
                    // special case for constitutive matrix (2x3)
        _matX[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts); //dNx/dx
        _matX[i][1].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts); //dNy/dy
        _matX[i][2].setVal(dNdy_[i] * a, 0 * nVerts, 1 * nVerts); //dNy/dx
        _matX[i][2].setVal(dNdx_[i] * a, 1 * nVerts, 2 * nVerts); //dNx/dy
// SET_MAT_ROW_SLICE(_matX[i], 0,dNdx_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dx
// SET_MAT_ROW_SLICE(_matX[i], 1,dNdy_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dy
// SET_MAT_ROW_SLICE(_matX[i], 2,dNdy_[i], a, 0 * nVerts, 1 * nVerts); //dNy/dx
// SET_MAT_ROW_SLICE(_matX[i], 2,dNdx_[i], a, 1 * nVerts, 2 * nVerts); //dNx/dy
                } else{
                // full matrix still unsure how to add constitutive for it
            _matX[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts); //dNx/dx
            _matX[i][1].setVal(dNdy_[i], 0 * nVerts, 1 * nVerts); //dNx/dy
            _matX[i][2].setVal(dNdx_[i], 1 * nVerts, 2 * nVerts); //dNy/dx
            _matX[i][3].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts); //dNy/dy
// SET_MAT_ROW_SLICE(_matX[i], 0, dNdx_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dx
// SET_MAT_ROW_SLICE(_matX[i], 1, dNdy_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dy
// SET_MAT_ROW_SLICE(_matX[i], 2, dNdx_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dx
// SET_MAT_ROW_SLICE(_matX[i], 3, dNdy_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dy
                   }
            } else if (ent.dim() == 3){


                if (elastic == true){
                    // special case for constitutive matrix (3x6)
_matX[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts); //dNx/dx
_matX[i][1].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts); //dNy/dy
_matX[i][2].setVal(dNdz_[i], 2 * nVerts, 3 * nVerts); //dNz/dz
_matX[i][3].setVal(dNdy_[i] * a, 0 * nVerts, 1 * nVerts);//dNy/dx
_matX[i][3].setVal(dNdx_[i] * a, 1 * nVerts, 2 * nVerts);//dNx/dy
_matX[i][4].setVal(dNdz_[i] * a, 1 * nVerts, 2 * nVerts);//dNz/dy
_matX[i][4].setVal(dNdy_[i] * a, 2 * nVerts, 3 * nVerts);//dNy/dz
_matX[i][5].setVal(dNdz_[i] * a, 0 * nVerts, 1 * nVerts);//dNz/dx
_matX[i][5].setVal(dNdx_[i] * a, 2 * nVerts, 3 * nVerts);//dNx/dz
// SET_MAT_ROW_SLICE(_matX[i],0,dNdx_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dx
// SET_MAT_ROW_SLICE(_matX[i],1,dNdy_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dy
// SET_MAT_ROW_SLICE(_matX[i],2,dNdz_[i], 1, 2 * nVerts, 3 * nVerts); //dNz/dz

// SET_MAT_ROW_SLICE(_matX[i],3,dNdy_[i], a, 0 * nVerts, 1 * nVerts);//dNy/dx
// SET_MAT_ROW_SLICE(_matX[i],3,dNdx_[i], a, 1 * nVerts, 2 * nVerts);//dNx/dy
// SET_MAT_ROW_SLICE(_matX[i],4,dNdz_[i], a, 1 * nVerts, 2 * nVerts);//dNz/dy
// SET_MAT_ROW_SLICE(_matX[i],4,dNdy_[i], a, 2 * nVerts, 3 * nVerts);//dNy/dz
// SET_MAT_ROW_SLICE(_matX[i],5,dNdz_[i], a, 0 * nVerts, 1 * nVerts);//dNz/dx
// SET_MAT_ROW_SLICE(_matX[i],5,dNdx_[i], a, 2 * nVerts, 3 * nVerts);//dNx/dz
                } else {
                // full matrix still unsure how to add constitutive for it
            _matX[i][0].setVal(dNdx_[i], 0 * nVerts, 1 * nVerts); //dNx/dx
            _matX[i][1].setVal(dNdy_[i], 0 * nVerts, 1 * nVerts); //dNx/dy
            _matX[i][2].setVal(dNdz_[i], 0 * nVerts, 1 * nVerts); //dNx/dz
            _matX[i][3].setVal(dNdx_[i], 1 * nVerts, 2 * nVerts); //dNy/dx
            _matX[i][4].setVal(dNdy_[i], 1 * nVerts, 2 * nVerts); //dNy/dy
            _matX[i][5].setVal(dNdz_[i], 1 * nVerts, 2 * nVerts); //dNy/dz
            _matX[i][6].setVal(dNdx_[i], 2 * nVerts, 3 * nVerts); //dNz/dx
            _matX[i][7].setVal(dNdy_[i], 2 * nVerts, 3 * nVerts); //dNz/dy
            _matX[i][8].setVal(dNdz_[i], 2 * nVerts, 3 * nVerts); //dNz/dz

// SET_MAT_ROW_SLICE(_matX[i],0, dNdx_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dx
// SET_MAT_ROW_SLICE(_matX[i],1, dNdy_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dy
// SET_MAT_ROW_SLICE(_matX[i],2, dNdz_[i], 1, 0 * nVerts, 1 * nVerts); //dNx/dz
// SET_MAT_ROW_SLICE(_matX[i],3, dNdx_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dx
// SET_MAT_ROW_SLICE(_matX[i],4, dNdy_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dy
// SET_MAT_ROW_SLICE(_matX[i],5, dNdz_[i], 1, 1 * nVerts, 2 * nVerts); //dNy/dz
// SET_MAT_ROW_SLICE(_matX[i],6, dNdx_[i], 1, 2 * nVerts, 3 * nVerts); //dNz/dx
// SET_MAT_ROW_SLICE(_matX[i],7, dNdy_[i], 1, 2 * nVerts, 3 * nVerts); //dNz/dy
// SET_MAT_ROW_SLICE(_matX[i],8, dNdz_[i], 1, 2 * nVerts, 3 * nVerts); //dNz/dz
                }
            }
        }
    }
    this->setValid(true);

    if (sum){
        this->integrate();
    }
    return *this;
}

template < > DLLEXPORT
ElementMatrix < double > & ElementMatrix < double >::grad(
                                    const MeshEntity & ent, Index order,
                                    bool elastic, bool sum, bool div,
                                    Index nCoeff, Index dof, Index dofOffset,
                                    bool kelvin){
    if (!this->valid() ||
        this->order() != order ||
        this->elastic() != elastic ||
        this->_ent != &ent ||
        this->_nCoeff != nCoeff){

        this->init(nCoeff, dof, dofOffset);
        this->grad(ent, order, elastic, sum, div, kelvin);
    }
    if (sum == true) this->integrate();
    return *this;
}

template < > DLLEXPORT ElementMatrix < double > &
ElementMatrix < double >::identity(const MeshEntity & ent, Index order,
                                   Index nCoeff, Index dofPerCoeff, Index dofOffset){


    this->init(nCoeff, dofPerCoeff, dofOffset);
    this->_ent = & ent;
    this->_order = order;

    this->_w = &IntegrationRules::instance().weights(ent.shape(), this->_order);
    this->_x = &IntegrationRules::instance().abscissa(ent.shape(), this->_order);

    const RVector &w = *this->_w;
    Index nVerts = ent.nodeCount();
    Index nRules(w.size());
    Index nCols = ent.dim() * this->_nCoeff;

    this->resize(nVerts*nCoeff, nCols);
    _matX.resize(nRules);

    for (Index i = 0; i < nRules; i ++ ){
        //** _matX stored transposed
        // transpose might be better?? check
        // fill per row is cheaper
        _matX[i].resize(nCols, nVerts*nCoeff);
        _matX[i].setZero();
        if (this->_nCoeff == ent.dim()){
            if (ent.dim() == 2){
                _matX[i].row(0).fill(1.0); // du_x/dx
                _matX[i].row(3).fill(1.0); // du_y/dy
            } else if (ent.dim() == 3){
                _matX[i].row(0).fill(1.0); // du_x/dx
                _matX[i].row(4).fill(1.0); // du_y/dy
                _matX[i].row(8).fill(1.0); // du_z/dz
            } else {
                THROW_TO_IMPL
            }
        } else {
            THROW_TO_IMPL
        }
    }
    this->setValid(true);
    // test if integration is needed for all cases
    // bool sum = true;
    // if (sum == true) this->integrate();

    return *this;
}

template < > DLLEXPORT RVector
ElementMatrix < double >::trace() const {
    THROW_TO_IMPL
    return 0;
}

template < > DLLEXPORT RSmallMatrix
ElementMatrix < double >::traceX() const {
    // __M
    RSmallMatrix  ret(_matX.size(), _matX[0].cols());
    // __MS(ret.rows(), ret.cols())
    Index i = 0;
    for (auto & m: _matX){
        if (m.size() == 2){
            // 1D vecspace no elastic
            ret[i] = m[0];
        } else if (m.size() == 4){
            // 2D vecspace no elastic
            ret[i] = m[0] + m[3];
        } else if (m.size() == 9){
            // 3D vecspace no elastic
            ret[i] = m[0] + m[4] + m[8];
        } else {
            __MS(m)
            THROW_TO_IMPL
        }
        i++;
    }
    // __MS(ret)
    return ret;
}


void _prepDot(const ElementMatrix < double > & A,
              const ElementMatrix < double > & B,
              ElementMatrix < double > & C){

    C.copyFrom(A, false);
    C.resize(A.rowIDs().size(), B.rowIDs().size());
    C.setIds(A.rowIDs(), B.rowIDs());

    if (A.order() != B.order()){
        __M
        log(Critical, "_prepDot. Elementmatrices need the same integration order",
            A.order(), ", ", B.order());
    }
}

void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         double b,
         ElementMatrix < double > & C){


    if (A.isIntegrated() && B.isIntegrated()){
        // shortcut for const per cell cases
        //THROW_TO_IMPL
    }

    bool verbose = false;
    _prepDot(A, B, C);

    const RVector &w = *A.w();
    Index nRules(w.size());

    double beta = 0.0;
    if (verbose){
        __MS("A: ", A.rows(), A.cols())
        __MS("B: ", B.rows(), B.cols())
        __MS("C: ", C.rows(), C.cols())
    }

    //** C = sum( A_i * B_i * b * w_i * dA)
    for (Index r = 0; r < nRules; r++){
        if (r > 0) beta = 1.0;

        const RSmallMatrix & Ai = A.matX()[r];
        const RSmallMatrix & Bi = B.matX()[r];
        RSmallMatrix & Ci = (*C.pMatX())[r];

        double wS = w[r] * A.entity()->size();

        if (verbose){
            __MS("Ci:\n ", Ci.rows(), Ci.cols())
            __MS("rule#: ", r, "b: " , b, "wS: ", wS, "w:", w[r],
                    "S:", A.entity()->size())
            __MS("Ai:\n ", Ai.rows(), Ai.cols(), "\n ", Ai)
            __MS("Bi:\n ", Bi.rows(), Bi.cols(), "\n ", Bi)
        }

        if (A.isDiv() || (Ai.rows() > 1 && Bi.rows()==1)){
            //
            // !! compare the following with the working
            //
            // RSmallMatrix  sA(1, Ai.cols());
            // sA[0] += Ai[0]; // v_x/dx
            // if (A.entity()->dim() == 2){
            //     ASSERT_VEC_SIZE(Ai, 2*2)
            //     sA[0] += Ai[3]; // v_y/dy
            // } else if (A.entity()->dim() == 3){
            //     ASSERT_VEC_SIZE(Ai, 3*3)
            //     sA[0] += Ai[4]; // v_y/dy
            //     sA[0] += Ai[8]; // v_z/dz
            // }
            //matTransMult(sA, B.matX()[r], *C.pMat(), c, beta);
            IndexArray ids;
            if (A.isDiv()){
                if (A.entity()->dim() == 1){
                    ids = IndexArray(std::vector < Index >({0}));
                } else if (A.entity()->dim() == 2){
                    ids = IndexArray(std::vector < Index >({0,3}));
                } else if  (A.entity()->dim() == 3){
                    ids = IndexArray(std::vector < Index >({0,4,8}));
                }
            } else {
                ids = range(Ai.rows());
            }

            // __MS(ids)
            // RSmallMatrix  Aii(Ai[0]);
            // for (Index i = 1; i < ids.size(); i ++){
            //     Aii += Ai[ids[i]];
            // }
            // double bi = 0;
            RSmallMatrix Aii(1, Ai.cols());
            Aii.row(0) = Ai.row(0);

            for (Index i = 1; i < ids.size(); i ++){
                Aii.row(0) += Ai.row(ids[i]);
            }

            if (verbose) __MS("Aii", Aii)

            Aii.transMult(Bi, Ci, b, 0.0);

            // __MS(Ci.rows(), Ci.cols());
                // THROW_TO_IMPL

        } else if (B.isDiv() || (Bi.rows() > 1 && Ai.rows() == 1)){
            // RSmallMatrix  sB(1, Bi.cols());
            // sB[0] += Bi[0]; // v_x/dx
            // if (B.entity()->dim() == 2){
            //     ASSERT_VEC_SIZE(Bi, 2*2)
            //     sB[0] += Bi[3]; // v_y/dy
            // } else if (B.entity()->dim() == 3){
            //     ASSERT_VEC_SIZE(Bi, 3*3)
            //     sB[0] += Bi[4]; // v_y/dy
            //     sB[0] += Bi[8]; // v_z/dz
            // }

            IndexArray ids;
            if (B.isDiv()){
                if (B.entity()->dim() == 1){
                    ids = IndexArray(std::vector < Index >({0}));
                } else if (B.entity()->dim() == 2){
                    ids = IndexArray(std::vector < Index >({0,3}));
                } else if (B.entity()->dim() == 3){
                    ids = IndexArray(std::vector < Index >({0,4,8}));
                }
            } else {
                ids = range(Bi.rows());
            }

            // __MS(ids)
            // RSmallMatrix Bii(Bi.row(0));
            RSmallMatrix Bii(1, Bi.cols());
            Bii.row(0) = Bi.row(0);

            for (Index i = 1; i < ids.size(); i ++){
                Bii.row(0) += Bi.row(ids[i]);
            }
            Ai.transMult(Bii, Ci, b, 0.0);

            // double bi = 0;
            // for (auto i: ids){
            //     // __MS(i)
            //     if (i > 0) bi = 1.0;
            //     RSmallMatrix  Bii(1, Bi[i].size());
            //     Bii[0] = Bi[i];
            //     matTransMult(Ai, Bii, Ci, b, bi);
            // }

            // matTransMult(Ai, sB, Ci, c, 0.0);
            //matTransMult(A.matX()[r], sB, *C.pMat(), c, beta);

        } else {
            if (verbose){
                __MS("no div transMult")
            }
            Ai.transMult(Bi, Ci, b, 0.0);
        }
        if (verbose){
            __MS("Ci:\n ", Ci.rows(), " ", Ci.cols(), "\n ", Ci)
        }

        // integration of pretransposed submatrices .. needs check! and sync with integrate
        if (beta == 1.0){
            *C.pMat() += Ci * wS;
        } else {
            *C.pMat() = Ci * wS;
        }
        if (verbose){
            __MS("C:", *C.pMat())
        }
    }
    // __MS("C:\n ", C.rows(), " ", C.cols(), "\n ", C)


    //# do we need submatrices after dot? yes e.g. for: (u*v+u*v)*u
    C.setValid(true);
    C.integrated(true);

    // __MS(C);
    // C*=0;

    // __MS(A.isIntegrated(),'\n', A);
    // __MS(A.isIntegrated(),'\n', A);
    // __MS(b);

    // matMult(A.mat(), B.mat(), *C.pMat(), b, 0.0);
    // __MS(C);
    // exit(1);
}
void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         const RSmallMatrix  & c,
         ElementMatrix < double > & C){
// __M
    _prepDot(A, B, C);

    bool reduceCe = false;
    if (c.rows() != A.cols() || c.cols() != B.cols()){
        if (A.elastic() and B.elastic() and A.entity()->dim() == 2 and c.cols()== 6){
            reduceCe = true;
            //## special case for 2D elastic with 3D elasticity tensor
            // __MS("special case for 2D elastic with 3D elasticity tensor")

        } else {

            __MS(c.rows(), c.cols())
            log(Critical, "Parameter matrix need to match Elementmatrix shapes: "
                "A:(", A.rows(), ",", A.cols(), ")",
                "B:(", B.rows(), ",", B.cols(), ")");

            return;
        }
    }

    const RVector &w = *A.w();
    (*C.pMat()) *= 0.0; // needed because matMult allways adds

    RSmallMatrix AtC;

    //** mat = sum A_i.T * C_i * B_i
    // log(Info, "A:(", A.rows(), ",", A.cols(), ")",
    //           "B:(", B.rows(), ",", B.cols(), ")");
// #if USE_EIGEN3
//     RSmallMatrix ce;
//     toEigenMatrix(c, ce);
// #else

// #endif

    // double beta = 0.0;

    for (Index i = 0; i < w.size(); i ++ ){
        // if (i > 0) beta = 1.0;

        const RSmallMatrix & Ai = A.matX()[i];
        const RSmallMatrix & Bi = B.matX()[i];
        RSmallMatrix & Ci = (*C.pMatX())[i];
        double wS = w[i] * A.entity()->size();

        AtC *= 0.0; // needed because matMult allways adds

        // __MS ("Ai:(", Ai.rows(), ",", Ai.cols(), ")",
        //       "Bi:(", Bi.rows(), ",", Bi.cols(), ")");

    // #if USE_EIGEN3
    // #else
    //     matTransMult(Ai, c, AtC, 1.0, 0.0);
    // #endif
        if (reduceCe){
            RSmallMatrix ce(3,3);
            for (Index i = 0; i < 2; i ++){
                for (Index j = 0; j < 2; j ++){
                    ce(i,j) = c(i,j);
                }
            }
            ce(2,2) = c(4,4);

            Ai.transMult(ce, AtC, 1.0, 0.0);
        } else {
            Ai.transMult(c, AtC, 1.0, 0.0);
        }


        AtC.mult(Bi, Ci, 1.0, 0.0);

        *C.pMat() += Ci * wS;
    }
    C.setValid(true);
    C.integrated(true);
}

void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         const Pos & c,
         ElementMatrix < double > & C){
THROW_TO_IMPL
}
void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         const RVector & c,
         ElementMatrix < double > & C){
THROW_TO_IMPL
}
void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         const PosVector & c,
         ElementMatrix < double > & C){
THROW_TO_IMPL
}
void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         const std::vector < RSmallMatrix  > & c,
         ElementMatrix < double > & C){
THROW_TO_IMPL
}
void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
         const FEAFunction & c,
         ElementMatrix < double > & C){
THROW_TO_IMPL
}

const ElementMatrix < double > dot(const ElementMatrix < double > & A,
                                   const ElementMatrix < double > & B){
    return dot(A, B, 1.0);
}
void dot(const ElementMatrix < double > & A,
         const ElementMatrix < double > & B,
               ElementMatrix < double > & ret){
    return dot(A, B, 1.0, ret);
}

#define DEFINE_DOT_MULT_WITH_RETURN(A_TYPE) \
ElementMatrix < double > dot(const ElementMatrix < double > & A, \
                             const ElementMatrix < double > & B, \
                             A_TYPE c){ \
    ElementMatrix < double > C; \
    dot(A, B, c, C); \
    return C; \
} \
ElementMatrix < double > mult(const ElementMatrix < double > & A, \
                                    A_TYPE b){ \
    ElementMatrix < double > C; \
    mult(A, b, C); \
    return C; \
}

DEFINE_DOT_MULT_WITH_RETURN(double)
DEFINE_DOT_MULT_WITH_RETURN(const RVector &)
DEFINE_DOT_MULT_WITH_RETURN(const Pos &)
DEFINE_DOT_MULT_WITH_RETURN(const PosVector &)
DEFINE_DOT_MULT_WITH_RETURN(const RSmallMatrix  &)
DEFINE_DOT_MULT_WITH_RETURN(const std::vector < RSmallMatrix  > &)
DEFINE_DOT_MULT_WITH_RETURN(const FEAFunction &)

#undef DEFINE_DOT_MULT_WITH_RETURN


// this *= const. Scalar
void mult(const ElementMatrix < double > & A,
          double f,
          ElementMatrix < double > & C){
    // __MS("** mult(A, d)")
    C.copyFrom(A, false);
    const RVector &w = *A.w();
    Index nRules(w.size());

    for (Index r = 0; r < nRules; r++){
        RSmallMatrix & iC = (*C.pMatX())[r];

        for (Index k = 0; k < iC.rows(); k ++){
            iC.row(k) *= f;

        }
        // __MS("r", r, rt)
    }
    C.setValid(true);
    C.integrate();
}
// this *= constant Pos
void mult(const ElementMatrix < double > & A,
          const Pos & f,
          ElementMatrix < double > & C){
    // __MS("** mult(A, p)")
// __MS(f)
    C.copyFrom(A, false);

    const PosVector &x = *A.x();
    // const RVector &w = *A.w();

    Index nRules(x.size());

    for (Index r = 0; r < nRules; r++){
        RSmallMatrix & iC = (*C.pMatX())[r];

        if (0 and iC.rows() == 1){
            __M
            iC.row(0) *= f.sum();
        } else {
            for (Index k = 0; k < iC.rows(); k ++){
                iC.row(k) *= f[k];
            }
        }
    }
    C.setValid(true);
    C.integrate();
}
// this *= scalar per quadrature
void mult(const ElementMatrix < double > & A,
          const RVector & b,
          ElementMatrix < double > & C){
    // __MS("** mult(A, rv)", b.size(), A.cols(), A.nCoeff(), A.elastic())
    // __MS(b)

    // if ((b.size() == A.cols() && A.cols() == A.nCoeff()) ||
    //     (b.size() == A.cols() && A.elastic())){
    // scale per component only for non gradients
    // only else case would be grad(v)* E0(ndarray(dim x dim)) -> ambiguous
    // if (b.size() == A.cols() && A.cols() == A.nCoeff()){
    if (b.size() == A.cols()){
        //** const scalar scale of matrix components
        // __MS("** mult(A, rv) with b.size() == A.cols()")
        C.copyFrom(A, false);
        //const PosVector &x = *A.x();
        Index nRules(C.w()->size());
        // __MS(A.rows(), A.cols(), b.size())

        for (Index r = 0; r < nRules; r++){
            RSmallMatrix & iC = (*C.pMatX())[r];
            // print("iC:", iC);
            for (Index k = 0; k < iC.rows(); k ++){
                // print('r', r, 'k', k, mr(k), b[r]);
                iC.row(k) *= b[k];
            }
        }
        C.integrate();
        return;
    }

    mult_d_q(A, b, C);
}
// vector per quadrature
void mult(const ElementMatrix < double > & A,
          const PosVector & f,
          ElementMatrix < double > & C){
    // result is no bilinear form, so keep it a rowMatrix check!!
    // __MS("** mult(A, pv)")
    C.copyFrom(A, false);
    const PosVector &x = *A.x();

    Index nRules(x.size());

    ASSERT_VEC_SIZE(f, nRules)
    ASSERT_VEC_SIZE(C.matX(), nRules)

    for (Index r = 0; r < nRules; r++){
        RSmallMatrix & iC = (*C.pMatX())[r];
        if (0 and iC.rows() == 1){
            __M
            iC.row(0) *= f[r].sum();
        } else {
            for (Index k = 0; k < iC.rows(); k ++){
                // __MS(r ," ", k ," ", f[r][k])
                iC.row(k) *= f[r][k];
                // iC.row(k) *= f[k][r];
            }
        }
        // __MS(iC)
    }

    C.setValid(true);
    C.integrate();
}
// constant Matrix
void mult(const ElementMatrix < double > & A,
          const RSmallMatrix  &  b,
          ElementMatrix < double > & C){
    // __MS("** mult(A, rm)")
    // __MS("b:\n", b)

    if (b.rows()*b.cols() == A.cols()){
        //** const scalar scale of each matrix components
        THROW_TO_IMPL
        //return mult(A, b.flatten(), C);
    }

    // result is no bilinear form, so keep it a rowMatrix

    C.copyFrom(A, false);

    const PosVector &x = *A.x();
    // const RVector &w = *A.w();

    Index nRules(x.size());


    if (b.rows() == A.matX().size() && b[0].size() == A.matX()[0].cols()){
        //** RVector per quadrature

        for (Index i = 0; i < nRules; i++){
            RSmallMatrix & Ci = (*C.pMatX())[i];
            const RSmallMatrix & Ai = A.matX()[i];
            Ci = Ai;
            Ci *= b[i];
        }
        C.integrate();
    } else {
        if (b.rows() != A.matX()[0].rows()){
            __MS("b:\n", b)
            __MS("A.matX()[0]:\n", A.matX()[0])
            log(Error, "Parameter matrix rows need to match Element sub matrix rows: ",
                A.matX()[0].rows());
            return;
        }

    // __MS(A.rows(), A.cols())
    // __MS(b.rows(), b.cols())
    // __MS(C.rows(), C.cols())
// #if USE_EIGEN3
//         RSmallMatrix be;
//         toEigenMatrix(b, be);
// #endif

        double beta = 0.0;
        for (Index i = 0; i < nRules; i++){
            if (i > 0) beta = 1.0;

            RSmallMatrix & Ci = (*C.pMatX())[i];
            const RSmallMatrix & Ai = A.matX()[i];
            // A.T * C
            Ci *= 0.0; // test and optimize me with C creation
            //matTransMult(Ai, b, Ci, 1.0, beta);
            // result is no bilinear form, so keep it a rowMatrix
// #if USE_EIGEN3
//             be.mult(Ai, Ci, 1.0, beta);
// #else
            b.mult(Ai, Ci, 1.0, beta);
// #endif
        }

        C.setValid(true);
        C.integrate(); // check if necessary
    }
}

// matrix per quadrature
void mult(const ElementMatrix < double > & A,
          const std::vector < RSmallMatrix  > & b,
          ElementMatrix < double > & C){
    // __MS("** mult(A, vmd)")
    // __MS("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    C.copyFrom(A, false);
    const PosVector &x = *A.x();

    Index nRules(x.size());

    ASSERT_VEC_SIZE(b, nRules)
    ASSERT_VEC_SIZE(C.matX(), nRules)

    double beta = 0.0;
    for (Index i = 0; i < nRules; i++){
        // if (i > 0) beta = 1.0; ## would be C +=

        RSmallMatrix & Ci = (*C.pMatX())[i];
        const RSmallMatrix & Ai = A.matX()[i];
        // A.T * C
        Ci *= 0.0; // test and optimize me with C creation

        // __MS(i, "/", nRules)
        // __MS("Ai:\n", Ai)
        // __MS("b[i]:\n", b[i])

        // Ai.transMult(b[i], Ci, 1.0, beta);
        b[i].mult(Ai, Ci, 1.0, beta);
// #endif
        // __MS(Ci)
    }

    C.setValid(true);

    //C.integrate();
    // __MS(C)
    // __MS("------------------------------------------------------------")
}

void mult(const ElementMatrix < double > & A, const FEAFunction & b,
          ElementMatrix < double > & C){
    // refactor with above
    // __MS(b.valueSize())

    if (b.valueSize() == 1){
        switch (b.getEvalOrder()){
            case 0: // cell center
                return mult(A, b.evalR1(A.entity()->center(), A.entity()), C);
                break;
            case 1:  // nodes
                THROW_TO_IMPL;
                break;
            case 2: { // quads
                RVector e;
                evaluateQuadraturePoints(*A.entity(), *A.x(), b, e);
                mult(A, e, C);
            } break;
            default:
                __M
                log(Error, "Eval order = ", b.getEvalOrder(), " is not defined.");
        }
    } else if (b.valueSize() == 3){
        switch (b.getEvalOrder()){
            case 0:
                return mult(A, b.evalR3(A.entity()->center(), A.entity()), C);
            case 1:
                THROW_TO_IMPL;
                break;
            case 2: {
                PosVector e;
                evaluateQuadraturePoints(*A.entity(), *A.x(), b, e);
                mult(A, e, C);
            } break;
            default:
                __M
                log(Error, "Eval order = ", b.getEvalOrder(), " is not defined.");
        }
    } else {
        switch (b.getEvalOrder()){
            case 0:
                return mult(A, b.evalRM(A.entity()->center(), A.entity()), C);
            case 1:
                THROW_TO_IMPL;
                break;
            case 2:{
                std::vector < RSmallMatrix  > e;
                evaluateQuadraturePoints(*A.entity(), *A.x(), b, e);
                mult(A, e, C);
            } break;
            default:
                __M
                log(Error, "Eval order = ", b.getEvalOrder(), " is not defined.");
        }
    }

    return;
}

#define DEFINE_DOT_MULT(A_TYPE) \
void mult_n(const ElementMatrix < double > & A, \
            A_TYPE b, ElementMatrix < double > & C){ \
    THROW_TO_IMPL \
}
DEFINE_DOT_MULT(const PosVector &)
DEFINE_DOT_MULT(const std::vector < RVector > &)
DEFINE_DOT_MULT(const std::vector < PosVector > &)
DEFINE_DOT_MULT(const std::vector < RSmallMatrix > &)
DEFINE_DOT_MULT(const std::vector < std::vector < RSmallMatrix  > > &)
#undef DEFINE_DOT_MULT

void mult_n(const ElementMatrix < double > & A, \
            const RVector & b, ElementMatrix < double > & C){
    // __MS("** mult_n(A, rv)")
    mult_d_n(A, b, C);
}
// this *= scalar per node
void mult_d_n(const ElementMatrix < double > & A,
              const RVector & b,
              ElementMatrix < double > & C){

    // __MS("** mult_s_n(A, rv)")
    C.copyFrom(A, false);
    Index nRules(C.w()->size());

    if (b.size() == A.rows()){
        //** b already sub of whole b
        for (Index r = 0; r < nRules; r++){
            RSmallMatrix & iC = (*C.pMatX())[r];
            // print("iC:", iC);
            for (Index k = 0; k < iC.rows(); k ++){
                // print('r', r, 'k', k, mr(k), b[r]);
                iC.row(k) *= b;
            }
        }
        C.integrate();
        return;
    } else if (b.size() == A.rows() / A.nCoeff()){
        //## b is already sub and scalar but A is on R3
        // optimize me!!
        Index nodeCount(A.rows() / A.nCoeff());
        for (Index r = 0; r < nRules; r++){
            RSmallMatrix & iC = (*C.pMatX())[r];
            // print("iC:", iC);
            for (Index k = 0; k < iC.rows(); k ++){
                // print('r', r, 'k', k, mr(k), b[r]);
                for (Index c = 0; c <= iC.cols(); c ++ ){
                    iC[k][c] *= b[c%nodeCount];
                }
            }
        }
        C.integrate();
        return;

    } else if (b.size() == A.nCoeff()*A.dofPerCoeff()){
        //## b is chained for all nodes and dims
        // optimize me!!
        return mult_d_n(A, b[A.rowIDs()], C);
        // //** b is whole b
        // for (Index r = 0; r < nRules; r++){
        //     RSmallMatrix & iC = (*C.pMatX())[r];
        //     // print("iC:", iC);
        //     for (Index k = 0; k < iC.rows(); k ++){
        //         // print('r', r, 'k', k, mr(k), b[r]);
        //         iC.row(k) *= b[A.rowIDs()];
        //     }
        // }
        // C.integrate();
    } else if (b.size() == A.dofPerCoeff()){
        //## b is for all nodes
        // optimize me!!
        for (Index r = 0; r < nRules; r++){
            RSmallMatrix & iC = (*C.pMatX())[r];
            // print("iC:", iC);
            for (Index k = 0; k < iC.rows(); k ++){
                // print('r', r, 'k', k, mr(k), b[r]);
                for (Index c = 0; c <= iC.cols(); c ++ ){
                    iC[k][c] *= b[A.rowIDs()[c]%A.dofPerCoeff()];
                }
            }
        }
        C.integrate();
        return;

    } else {

        print("b:", b.size(), "A:", A.cols(), "x",  A.rows());
        print("nCoeff:", A.nCoeff(), "dofperCoeff:", A.dofPerCoeff());
        print("A:", A);
        print("b:", b);
        THROW_TO_IMPL
    }
}


void mult_d_q(const ElementMatrix < double > & A,
              const RVector & b,
              ElementMatrix < double > & C){
    //** scalar per quadrature
    // __MS("** mult_s_q(A, rv)")
    C.copyFrom(A, false);
    Index nRules(C.w()->size());


    ASSERT_VEC_SIZE(b, nRules)
    ASSERT_VEC_SIZE(C.matX(), nRules)

    for (Index r = 0; r < nRules; r++){
        RSmallMatrix & mr = (*C.pMatX())[r];
        // print("iC:", mr);
        for (Index k = 0; k < mr.rows(); k ++){
            // print('r', r, 'k', k, mr(k), b[r]);
            // MAT_ROW_IMUL(mr, k, b[r])
            mr[k] *= b[r];
        }
    }

    C.setValid(true);
    C.integrate();
}


//******************************************************************************
// LINEAR-FORM -- integration template -- per CELL || QUADRATURE
//******************************************************************************
#define INTEGRATE_LINFORM(_F)                              \
    const RVector &w = *this->_w;                          \
    Index nRules(w.size());                                \
    ASSERT_VEC_SIZE(this->matX(), nRules)                  \
    RVector rt(this->rows());                              \
    Index rowStep = 1;                                     \
    Index maxRows = this->cols();                          \
    if (this->nCoeff() == 2 and this->cols() == 4){        \
        rowStep = 3;                                       \
    } else if (this->nCoeff() == 2 and this->cols() == 3){ \
        maxRows = entity()->dim();  /*elastic mapping */   \
    } else if (this->nCoeff() == 3 and this->cols() == 9){ \
        rowStep = 4;                                       \
    } else if (this->nCoeff() == 3 and this->cols() == 6){ \
        maxRows = entity()->dim(); /*elastic mapping */    \
    }                                                      \
    for (Index q = 0; q < nRules; q++){                    \
        const RSmallMatrix &mr(this->_matX[q]);            \
        for (Index k = 0; k < maxRows; k += rowStep){      \
            if (q == 0 and k == 0){                        \
                rt = mr.row(k) * w[q] _F;                  \
            } else {                                       \
                rt += mr.row(k) * w[q] _F;                 \
            }                                              \
        }                                                  \
    }                                                      \

//******************************************************************************
// LINEAR-FORM -- integration -- const Scalar -- f in R evalOnCells
//******************************************************************************
template < > void
ElementMatrix < double >::integrate(double f,
                                    RVector & R, double scale) const {
    this->integrate();
    R.add(*this, f*scale);
    // INTEGRATE_LINFORM()
    // rt *= this->_ent->size() * f * scale;
    // R.addVal(rt, this->rowIDs());
}
//******************************************************************************
// LINEAR-FORM -- integration -- const Vector -- ??
//******************************************************************************
template < > void
ElementMatrix < double >::integrate(const Pos & f,
                                    RVector & R, double scale) const {
    INTEGRATE_LINFORM(*f[k])
    rt *= this->_ent->size() * scale;
    R.addVal(rt, this->rowIDs());
}
//******************************************************************************
// LINEAR-FORM -- integration -- scalar per QUADRATURE -- f in R evalOnQuads
//******************************************************************************
template < > void
ElementMatrix < double >::integrate(const RVector & f,
                                    RVector & R, double scale) const {
    // only integrate if f fits size
    if (f.size() == this->_w->size()){
        INTEGRATE_LINFORM(*f[q])
        rt *= this->_ent->size() * scale;
        R.addVal(rt, this->rowIDs());
    }
}
//******************************************************************************
// LINEAR-FORM -- integration -- vector per QUADRATURE -- f in R^n evalOnQuads
//******************************************************************************
template < >
void ElementMatrix < double >::integrate(const PosVector & f,
                                         RVector & R, double scale) const {
    // __MS("** M.integrate(rv, ->R)", this->_ent->size(), " ", scale)
    // __MS("** M.integrate(rv, ->R): ", R)
    ASSERT_VEC_SIZE(f, this->_w->size())
    INTEGRATE_LINFORM(*f[q][k]) // #orig
    rt *= this->_ent->size() * scale;
    R.addVal(rt, this->rowIDs());
}
//******************************************************************************
// LINEAR-FORM -- integration -- DensMatrix per QUADRATURE -- f in R^n evalOnQuads
//******************************************************************************
template < > void
ElementMatrix < double >::integrate(const std::vector< RSmallMatrix > & f,
                                    RVector & R, double scale) const {
    ASSERT_VEC_SIZE(f, this->_w->size())

    INTEGRATE_LINFORM(*sum(f[q][k])) // #orig
    rt *= this->_ent->size() * scale;
    R.addVal(rt, this->rowIDs());
}

//******************************************************************************
// LINEAR-FORM -- integration -- FALLBACK per CELL || QUADRATURE
//******************************************************************************
#define DEFINE_INTEGRATOR_LF(A_TYPE) \
template < > \
void ElementMatrix < double >::integrate(A_TYPE f, \
                                         RVector & R, double scale) const { \
    THROW_TO_IMPL \
}
DEFINE_INTEGRATOR_LF(const RSmallMatrix & )  // const Matrix
DEFINE_INTEGRATOR_LF(const FEAFunction &)
#undef DEFINE_INTEGRATOR_LF

//*****************************************************************************
// LINEAR-FORM -- integration  -- scalar per NODE
//*****************************************************************************
template < > void
ElementMatrix < double >::integrate_n(const RVector & f,
                                      RVector & R, double scale) const {

    ASSERT_VEC_SIZE(f, this->dofPerCoeff())
    this->integrate();

    if (scale != 1.0) {
        R.add(*this, f*scale);
    } else {
        R.add(*this, f);
    }
}
//*****************************************************************************
// LINEAR-FORM -- integration  -- FALLBACK per NODE
//*****************************************************************************
#define DEFINE_INTEGRATOR_LF_N(A_TYPE) \
template < > \
void ElementMatrix < double >::integrate_n(A_TYPE f, \
                                           RVector & R, double scale) const { \
    THROW_TO_IMPL \
}
DEFINE_INTEGRATOR_LF_N(const PosVector & )                                   // unused
DEFINE_INTEGRATOR_LF_N(const std::vector< RSmallMatrix > &)                  // unused
DEFINE_INTEGRATOR_LF_N(const std::vector< RVector > &)                       // unused
DEFINE_INTEGRATOR_LF_N(const std::vector< PosVector > & )                    // unused
DEFINE_INTEGRATOR_LF_N(const std::vector< std::vector< RSmallMatrix  > > &)  // unused
#undef DEFINE_INTEGRATOR_LF_N


//*****************************************************************************
// BILINEAR-FORM -- Per CELL and QAUD integration
//*****************************************************************************
#define DEFINE_INTEGRATOR(A_TYPE) \
template < > \
void ElementMatrix < double >::integrate(const ElementMatrix < double > & B, \
                                A_TYPE f, SparseMatrixBase & A, double scale) \
                       const { \
    ElementMatrix < double > dAB; \
    dot(*this, B, f, dAB); \
    A.add(dAB, scale); \
}

DEFINE_INTEGRATOR(double)   // const scalar
DEFINE_INTEGRATOR(const RSmallMatrix  &)  // const Matrix
#undef DEFINE_INTEGRATOR

#define DEFINE_INTEGRATOR(A_TYPE) \
template < > \
void ElementMatrix < double >::integrate(const ElementMatrix < double > & B, \
                                A_TYPE f, SparseMatrixBase & A, double scale) \
                       const { \
    ElementMatrix < double > dAB; \
    ElementMatrix < double > mf; \
    GIMLI::mult(*this, f, mf); \
    dot(mf, B, 1, dAB); \
    A.add(dAB, scale); \
}

DEFINE_INTEGRATOR(const RVector &)  // scalar for each quadr
DEFINE_INTEGRATOR(const std::vector< RSmallMatrix  > &)// matrix for each quadrs
#undef DEFINE_INTEGRATOR

template < >
void ElementMatrix < double >::integrate(const ElementMatrix < double > & B, \
                       const Pos & v, SparseMatrixBase & A, double scale) const {
    ElementMatrix < double > dAB;
    ElementMatrix < double > vm;
    GIMLI::mult(B, v, vm);
    dot(*this, vm, 1, dAB);
    A.add(dAB, scale);
}

template < >
void ElementMatrix < double >::integrate(const ElementMatrix < double > & R, \
                       const PosVector & p, SparseMatrixBase & A, double scale) const {
THROW_TO_IMPL
}


template < >
void ElementMatrix < double >::integrate(const ElementMatrix < double > & B,
                                         const FEAFunction & f,
                                         SparseMatrixBase & A, double scale) const {
THROW_TO_IMPL
}

//
void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                              const FEAFunction & f,
                              RVector & ret){
THROW_TO_IMPL
    // ret.clear();
    // const PosVector *x;
    // for (auto &cell: mesh.cells()){
    //     x = &IntegrationRules::instance().abscissa(cell->shape(), order);
    //     for (Index i = 0; i < x->size(); i ++){

    //         if (f.valueSize() == 1){
    //             ret.push_back(f.evalR1(cell->shape().xyz((*x)[i]), cell));
    //         } else {
    //             __M
    //             log(Critical, "expecting FEAFunction with valueSize==1.",
    //                 f.valueSize());
    //         }
    //     }
    // }
}
void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                              const FEAFunction & f,
                              PosVector & ret){
    THROW_TO_IMPL
    // ret.clear();
    // const PosVector *x;
    // for (auto &cell: mesh.cells()){
    //     x = &IntegrationRules::instance().abscissa(cell->shape(), order);
    //     for (Index i = 0; i < x->size(); i ++){

    //         if (f.valueSize() > 1){
    //             ret.push_back(f.evalR3(cell->shape().xyz((*x)[i]), cell));
    //         } else {
    //             __M
    //             log(Critical, "expecting FEAFunction with valueSize==2 or 3",
    //                 f.valueSize());
    //         }
    //     }
    // }
}

void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                              const FEAFunction & f,
                              std::vector< RSmallMatrix  > & ret){
    THROW_TO_IMPL
}

void evaluateQuadraturePoints(const MeshEntity & ent, const PosVector & x,
                              const FEAFunction & f, RVector & ret){
    ret.resize(x.size());
    for (Index i = 0; i < x.size(); i ++){
        ret[i] = f.evalR1(ent.shape().xyz(x[i]), &ent);
    }
}

void evaluateQuadraturePoints(const MeshEntity & ent, const PosVector & x,
                              const FEAFunction & f, PosVector & ret){
    ret.resize(x.size());
    for (Index i = 0; i < x.size(); i ++){
        ret[i] = f.evalR3(ent.shape().xyz(x[i]), &ent);
    }
}

void evaluateQuadraturePoints(const MeshEntity & ent, const PosVector & x,
                              const FEAFunction & f,
                              std::vector < RSmallMatrix  > & ret){
    ret.resize(x.size());
    for (Index i = 0; i < x.size(); i ++){
        ret[i] = f.evalRM(ent.shape().xyz(x[i]), &ent);
    }
}

template < class ReturnType >
void evaluateQuadraturePoints_(const Mesh & mesh, Index order,
                               const FEAFunction & f, ReturnType & ret){
    ret.resize(mesh.cellCount());

    const PosVector *x;
    for (auto &cell: mesh.cells()){
        x = &IntegrationRules::instance().abscissa(cell->shape(), order);
        evaluateQuadraturePoints(*cell, *x, f, ret[cell->id()]);
    }
}

void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                              const FEAFunction & f,
                              std::vector< RVector > & ret){
    evaluateQuadraturePoints_(mesh, order, f, ret);
}

void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                              const FEAFunction & f,
                              std::vector< PosVector > & ret){
    evaluateQuadraturePoints_(mesh, order, f, ret);
}

void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                              const FEAFunction & f,
                              std::vector< std::vector< RSmallMatrix  > > & ret){
    evaluateQuadraturePoints_(mesh, order, f, ret);
}

void sym(const ElementMatrix < double > & A, ElementMatrix < double > & B){
    B.copyFrom(A, false);

    for (auto &m : *B.pMatX()){
        if (m.rows() == 1){
        } else if (m.rows() == 4){
            // __MS(m)
            m.row(1) = 0.5*m.row(1) + 0.5*m.row(2);
            m.row(2) = m.row(1);
            // __MS(m)
        } else if (m.rows() == 9) {
            // _matX[i](0).setVal(dNdx_[i], 0 * nVerts, 1 * nVerts); //dNx/dx
            // _matX[i](1).setVal(dNdy_[i], 0 * nVerts, 1 * nVerts); //dNx/dy
            // _matX[i](2).setVal(dNdz_[i], 0 * nVerts, 1 * nVerts); //dNx/dz

            // _matX[i](3).setVal(dNdx_[i], 1 * nVerts, 2 * nVerts); //dNy/dx
            // _matX[i](4).setVal(dNdy_[i], 1 * nVerts, 2 * nVerts); //dNy/dy
            // _matX[i][5].setVal(dNdz_[i], 1 * nVerts, 2 * nVerts); //dNy/dz

            // _matX[i][6].setVal(dNdx_[i], 2 * nVerts, 3 * nVerts); //dNz/dx
            // _matX[i][7].setVal(dNdy_[i], 2 * nVerts, 3 * nVerts); //dNz/dy
            // _matX[i][8].setVal(dNdz_[i], 2 * nVerts, 3 * nVerts); //dNz/dz

            m.row(1) = 0.5 * m.row(1) + 0.5 * m.row(3);
            m.row(3) = m.row(1);
            m.row(2) = 0.5 * m.row(2) + 0.5 * m.row(6);
            m.row(6) = m.row(2);
            m.row(5) = 0.5 * m.row(5) + 0.5 * m.row(7);
            m.row(7) = m.row(5);
        } else {
            __MS(A)
            log(Critical, "Don't not how to symetrize A");
        }

    }
    B.setValid(true);
    B.integrated(false);
}


ElementMatrix < double > sym(const ElementMatrix < double > & A){
    ElementMatrix < double > B;
    sym(A, B);
    return B;
}


void trace(const ElementMatrix < double > & A, ElementMatrix < double > & B){
    B.copyFrom(A, false);
    RVector tr(A.rows());

    for (auto &m : *B.pMatX()){
        if (m.rows() == 1){
        } else if (m.rows() == 4){
            tr = m.row(0) + m.row(3);
            m *= 0.0;
            m.row(0) = tr;
            m.row(3) = tr;
        } else if (m.rows() == 9) {
            tr = m.row(0) + m.row(4) + m.row(8);
            m *= 0.0;
            m.row(0) = tr;
            m.row(4) = tr;
            m.row(8) = tr;
        } else {
            __MS(A)
            log(Critical, "Don't not how to trace A");
        }
    }
    B.setValid(true);
    B.integrated(false);
}


ElementMatrix < double > trace(const ElementMatrix < double > & A){
    ElementMatrix < double > B;
    trace(A, B);
    return B;
}


template < class Vec >
void createForceVectorPerCell_(const Mesh & mesh, Index order, RVector & ret,
                        const Vec & a, Index nCoeff, Index dofOffset){
    if (nCoeff > 3){
        __M
        log(Critical, "Number of coefficients need to be lower then 4");
    }
    Index dof = mesh.nodeCount() * nCoeff;
    ret.resize(dof);

    ElementMatrix < double > u;
    for (auto &cell: mesh.cells()){
        u.pot(*cell, order, true, nCoeff, mesh.nodeCount(), dofOffset);

        if (a.size() == 1){
            ret.add(u, a[0]);
        } else if (a.size() == mesh.cellCount()){
            ret.add(u, a[cell->id()]);
        } else {
            __M
            log(Critical, "Number of cell coefficients (",a.size(),") does not"
                "match cell count:",  mesh.cellCount());
        }
    }
}
template < class Vec >
void createForceVectorMult_(const Mesh & mesh, Index order, RVector & ret,
                            const Vec & a, Index nCoeff, Index dofOffset){
    if (nCoeff > 3){
        log(Critical, "Number of coefficients need to be lower then 4");
    }
    Index dof = mesh.nodeCount() * nCoeff;
    ret.resize(dof);

    ElementMatrix < double > u;
    ElementMatrix < double > ua;

    for (auto &cell: mesh.cells()){
        u.pot(*cell, order, true, nCoeff, mesh.nodeCount(), dofOffset);

        if (a.size() == 1 && mesh.cellCount() != 1){
            createForceVectorPerCell_(mesh, order, ret,
                                      a[0], nCoeff, dofOffset);
        } else if (a.size() == mesh.cellCount()){
            mult(u, a[cell->id()], ua);
            ret.add(ua);
        } else {
            __M
            log(Critical, "Number of per cell coefficients (",a.size(),") does not"
                "match cell count:",  mesh.cellCount());
        }
    }
}
template < class Vec >
void createMassMatrixPerCell_(const Mesh & mesh, Index order,
                              RSparseMapMatrix & ret, const Vec & a,
                              Index nCoeff, Index dofOffset){
    if (nCoeff > 3){
        log(Critical, "Number of coefficients need to be lower then 4");
    }
    ElementMatrix < double > u;
    ElementMatrix < double > uu;

    for (auto &cell: mesh.cells()){
        u.pot(*cell, order, true, nCoeff, mesh.nodeCount(), dofOffset);

        if (a.size() == 1){
            dot(u, u, a[0], uu);
        } else if (a.size() == mesh.cellCount()){
            dot(u, u, a[cell->id()], uu);
        } else {
            __M
            log(Critical, "Number of cell coefficients (",a.size(),") does not"
                "match cell count:",  mesh.cellCount());
        }
        ret.add(uu);
    }
}
template < class Vec >
void createMassMatrixMult_(const Mesh & mesh, Index order,
                              RSparseMapMatrix & ret, const Vec & a,
                              Index nCoeff, Index dofOffset){
    if (nCoeff > 3){
        log(Critical, "Number of coefficients need to be lower then 4");
    }
    ElementMatrix < double > u;
    ElementMatrix < double > ua;
    ElementMatrix < double > uau;

    for (auto &cell: mesh.cells()){
        u.pot(*cell, order, true, nCoeff, mesh.nodeCount(), dofOffset);

        if (a.size() == 1 && mesh.cellCount() != 1){
            createMassMatrixPerCell_(mesh, order, ret, a[0], nCoeff, dofOffset);
        } else if (a.size() == mesh.cellCount()){
            mult(u, a[cell->id()], ua);
            dot(ua, u, 1.0, uau);
            ret.add(uau);
        } else {
            __M
            log(Critical, "Number of cell coefficients (",a.size(),") does not"
                "match cell count:",  mesh.cellCount());
        }
    }
}
template < class Vec >
void createStiffnessMatrixPerCell_(const Mesh & mesh, Index order,
                                   RSparseMapMatrix & ret, const Vec & a,
                                   Index nCoeff, Index dofOffset,
                                   bool elastic, bool kelvin){
    if (nCoeff > 3){
        __M;
        log(Critical, "Number of coefficients need to be lower then 4");
    }

    ElementMatrix < double > du;
    ElementMatrix < double > dudu;

    for (auto &cell: mesh.cells()){
        //#bool elastic, bool sum, bool div,
        du.grad(*cell, order,
                elastic, false, false,
                nCoeff, mesh.nodeCount(), dofOffset, kelvin);

        if (a.size() == 1){
            dot(du, du, a[0], dudu);
        } else if (a.size() == mesh.cellCount()){
            dot(du, du, a[cell->id()], dudu);
        } else {
            __M;
            log(Critical, "Number of cell coefficients (",a.size(),") does not "
                "match cell count:",  mesh.cellCount());
        }
        ret.add(dudu);
    }
}

template < class Vec >
void createStiffnessMatrixMult_(const Mesh & mesh, Index order,
                                   RSparseMapMatrix & ret, const Vec & a,
                                   Index nCoeff, Index dofOffset,
                                   bool elastic, bool kelvin){
    if (nCoeff > 3){
        __M;
        log(Critical, "Number of coefficients need to be lower then 4");
    }
    ElementMatrix < double > du;
    ElementMatrix < double > dua;
    ElementMatrix < double > duadu;

    //#bool elastic, bool sum, bool div,

    for (auto &cell: mesh.cells()){
        du.grad(*cell, order,
                elastic, false, false,
                nCoeff, mesh.nodeCount(), dofOffset, kelvin);
        if (a.size() == 1 && mesh.cellCount() != 1){
            createStiffnessMatrixPerCell_(mesh, order, ret, a[0],
                                          nCoeff, dofOffset, elastic, kelvin);
        } else if (a.size() == mesh.cellCount()){
            mult(du, a[cell->id()], dua);
            dot(dua, du, 1.0, duadu);
            ret.add(duadu);
        } else {
            __M;
            log(Critical, "Number of cell coefficients (",a.size(),") does not"
                "match cell count:",  mesh.cellCount());
        }
    }
}

//** IMPL constants
#define DEFINE_CREATE_SCALAR_IMPL(A_TYPE, V_TYPE) \
void createForceVector(const Mesh & mesh, Index order, RVector & ret,\
                       A_TYPE a, Index nCoeff, Index dofOffset){\
    createForceVectorPerCell_(mesh, order, ret, V_TYPE(1,a),\
                              nCoeff, dofOffset);\
}\
void createMassMatrix(const Mesh & mesh, Index order, \
                      RSparseMapMatrix & ret, A_TYPE a, \
                      Index nCoeff, Index dofOffset){\
    createMassMatrixPerCell_(mesh, order, ret, V_TYPE(1,a),\
                             nCoeff, dofOffset);\
}\
void createStiffnessMatrix (const Mesh & mesh, Index order, \
                            RSparseMapMatrix & ret, A_TYPE a, \
                            Index nCoeff, Index dofOffset, \
                            bool elastic, bool kelvin){\
    createStiffnessMatrixPerCell_(mesh, order, ret, V_TYPE(1,a),\
                                  nCoeff, dofOffset, elastic, kelvin);\
}

DEFINE_CREATE_SCALAR_IMPL(double, RVector)
DEFINE_CREATE_SCALAR_IMPL(const Pos &, PosVector)
#undef DEFINE_CREATE_SCALAR_IMPL

//** IMPL constant matrix
void createForceVector(const Mesh & mesh, Index order, RVector & ret,
                       const RSmallMatrix  & a, Index nCoeff, Index dofOffset){
    std::vector < RSmallMatrix  > aM(1);
    aM[0] = a;
    createForceVectorPerCell_(mesh, order, ret, aM, nCoeff, dofOffset);
}
void createMassMatrix(const Mesh & mesh, Index order, RSparseMapMatrix & ret,
                      const RSmallMatrix  & a, Index nCoeff, Index dofOffset){
    std::vector < RSmallMatrix  > aM(1);
    aM[0] = a;
    createMassMatrixPerCell_(mesh, order, ret, aM, nCoeff, dofOffset);
}
void createStiffnessMatrix (const Mesh & mesh, Index order,
                            RSparseMapMatrix & ret, const RSmallMatrix  & a,
                            Index nCoeff, Index dofOffset,
                            bool elastic, bool kelvin){
    std::vector < RSmallMatrix  > aM(1);
    aM[0] = a;
    createStiffnessMatrixPerCell_(mesh, order, ret, aM,
                                  nCoeff, dofOffset, elastic, kelvin);
}

//** IMPL per cell values
#define DEFINE_CREATE_PERCELL_IMPL(A_TYPE) \
void createForceVector(const Mesh & mesh, Index order, RVector & ret, \
                       A_TYPE a, Index nCoeff, Index dofOffset){ \
    createForceVectorPerCell_(mesh, order, ret, a, nCoeff, dofOffset); \
} \
void createMassMatrix(const Mesh & mesh, Index order, \
                      RSparseMapMatrix & ret, A_TYPE a, \
                      Index nCoeff, Index dofOffset){\
    createMassMatrixPerCell_(mesh, order, ret, a, nCoeff, dofOffset);\
}\
void createStiffnessMatrix (const Mesh & mesh, Index order, \
                            RSparseMapMatrix & ret, A_TYPE a, \
                            Index nCoeff, Index dofOffset, \
                            bool elastic, bool kelvin){\
    createStiffnessMatrixPerCell_(mesh, order, ret, a, \
                                  nCoeff, dofOffset, elastic, kelvin);\
}

DEFINE_CREATE_PERCELL_IMPL(const RVector &)
DEFINE_CREATE_PERCELL_IMPL(const PosVector &)
DEFINE_CREATE_PERCELL_IMPL(const std::vector< RSmallMatrix  > &)
#undef DEFINE_CREATE_PERCELL_IMPL


//** IMPL per cell per quadrature values
#define DEFINE_CREATE_FORCE_VECTOR_IMPL(A_TYPE) \
void createForceVector(const Mesh & mesh, Index order, \
                       RVector & ret, A_TYPE a, \
                       Index nCoeff, Index dofOffset){ \
    createForceVectorMult_(mesh, order, ret, a, nCoeff, dofOffset); \
} \
void createMassMatrix(const Mesh & mesh, Index order, \
                      RSparseMapMatrix & ret, A_TYPE a, \
                      Index nCoeff, Index dofOffset){ \
    createMassMatrixMult_(mesh, order, ret, a, nCoeff, dofOffset);\
} \
void createStiffnessMatrix(const Mesh & mesh, Index order, \
                           RSparseMapMatrix & ret, A_TYPE a, \
                           Index nCoeff, Index dofOffset, \
                           bool elastic, bool kelvin){ \
    createStiffnessMatrixMult_(mesh, order, ret, a, \
                               nCoeff, dofOffset, elastic, kelvin);\
}

DEFINE_CREATE_FORCE_VECTOR_IMPL(const std::vector< RVector > &)
DEFINE_CREATE_FORCE_VECTOR_IMPL(const std::vector< PosVector > &)
DEFINE_CREATE_FORCE_VECTOR_IMPL(const std::vector< std::vector < RSmallMatrix  > > &)
#undef DEFINE_CREATE_FORCE_VECTOR_IMPL

//** IMPL fallbacks
#define DEFINE_CREATE_FORCE_VECTOR_IMPL(A_TYPE) \
void createForceVector(const Mesh & mesh, Index order, \
                       RVector & ret, A_TYPE a, \
                       Index nCoeff, Index dofOffset){ \
THROW_TO_IMPL \
} \
void createMassMatrix(const Mesh & mesh, Index order, \
                      RSparseMapMatrix & ret, A_TYPE a, \
                      Index nCoeff, Index dofOffset){ \
THROW_TO_IMPL \
                      } \
void createStiffnessMatrix(const Mesh & mesh, Index order, \
                           RSparseMapMatrix & ret, A_TYPE a, \
                           Index nCoeff, Index dofOffset, \
                           bool elastic, bool kelvin){ \
THROW_TO_IMPL \
}
DEFINE_CREATE_FORCE_VECTOR_IMPL(const FEAFunction &)

#undef DEFINE_CREATE_FORCE_VECTOR_IMPL


void createAdvectionMatrix(const Mesh & mesh, Index order,
                                     RSparseMapMatrix & ret,
                                     const PosVector & vel,
                                     Index dofOffset){
             THROW_TO_IMPL
}

} // namespace GIMLI
