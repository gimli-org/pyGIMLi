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

#ifndef _GIMLI_INTEGRATION__H
#define _GIMLI_INTEGRATION__H

#include "gimli.h"

namespace GIMLI{

/*! */
class DLLEXPORT IntegrationRules : public Singleton< IntegrationRules >{
public:
    friend class Singleton< IntegrationRules >;

    /*! Return Gauss-Legendre quadrature point upto order <10. */
    const R3Vector & gauAbscissa(Index order) const;
    
    /*! Return Gauss-Legendre quadrature weights upto order <10. */
    const RVector & gauWeights(Index order) const;

    /*!
     * Generic quadrature positions for a triangle based on Gauss-Legendre quadrature
     H. T. RATHOD1*, K. V. NAGARAJA2, B. VENKATESUDU3 AND N. L. RAMESH4.
     Gauss Legendre quadrature over a triangle.
     J. Indian Inst. Sci., Sept.-Oct. 2004, 84, 183-188
    */
    const R3Vector & triGLAbscissa(Index order) const;

    /*!
     * Generic quadrature weights for a triangle based on Gauss-Legendre quadrature
     H.T. RATHOD, K. V. NAGARAJA, B. VENKATESUDU AND N. L. RAMESH.
     Gauss Legendre quadrature over a triangle.
     J. Indian Inst. Sci., Sept.-Oct. 2004, 84, 183-188
    */
    const RVector & triGLWeights(Index order) const;

    const R3Vector & edgAbscissa(Index order) const;
    const RVector & edgWeights(Index order) const;

    const R3Vector & triAbscissa(Index order) const;
    const RVector & triWeights(Index order) const;

    const R3Vector & tetAbscissa(Index order) const;
    const RVector & tetWeights(Index order) const;

    const R3Vector & quaAbscissa(Index order) const;
    const RVector & quaWeights(Index order) const;

    const R3Vector & hexAbscissa(Index order) const;
    const RVector & hexWeights(Index order) const;

    const R3Vector & priAbscissa(Index order) const;
    const RVector & priWeights(Index order) const;
    /*!
     * Return Gauss-Legendre quadrature positions for a given shape of the \ref MeshEntity upto order 10
     */
    const R3Vector & abscissa(const Shape & shape, uint order) const;
    /*!
     * Return Gauss-Legendre quadrature weights for a given shape of the \ref MeshEntity upto order 10
     */
    const RVector & weights(const Shape & shape, uint order) const;

    /*! Set whether triangle integration use Gauss Legendre polynomials (up to order 9) or
     * native triangle coordinates (up to order 5). Default is GL == true. */
    inline void setTriUseGaussLegendre(bool use){ triUseGaussLegendre_ = use;}
    inline bool triUseGaussLegendre() const { return triUseGaussLegendre_; }

private:
    /*! Default constructor */
    IntegrationRules();

    /*! Default destructor */
    ~IntegrationRules();

protected:

    void initGau_();
    void initTriGL_();

    void initEdg_();
    void initTri_();
    void initTet_();
    void initQua_();
    void initHex_();
    void initPri_();

    bool triUseGaussLegendre_;

    std::vector < R3Vector > gauAbscissa_;
    std::vector < RVector  > gauWeights_;

    /*!Quadrature positions for a triangle based on Gauss-Legendre quadrature */
    std::vector < R3Vector > triGLAbscissa_;
    /*!Quadrature weight for a triangle based on Gauss-Legendre quadrature */
    std::vector < RVector  > triGLWeights_;

    std::vector < R3Vector > edgAbscissa_;
    std::vector < RVector  > edgWeights_;

    std::vector < R3Vector > triAbscissa_;
    std::vector < RVector  > triWeights_;

    std::vector < R3Vector > tetAbscissa_;
    std::vector < RVector  > tetWeights_;

    std::vector < R3Vector > quaAbscissa_;
    std::vector < RVector  > quaWeights_;

    std::vector < R3Vector > hexAbscissa_;
    std::vector < RVector  > hexWeights_;

    std::vector < R3Vector > priAbscissa_;
    std::vector < RVector  > priWeights_;

};

} // namespace GIMLI{

#endif // _GIMLI_INTEGRATION__H
