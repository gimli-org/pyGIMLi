/******************************************************************************
 *   Copyright (C) 2009-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_DC1DMODELLING__H
#define _GIMLI_DC1DMODELLING__H

#include "gimli.h"
#include "mesh.h"
#include "meshgenerators.h"
#include "modellingbase.h"
#include "vectortemplates.h"

namespace GIMLI{

//! DC (direct current) 1D modelling
/*! Classical DC 1D forward operator for given resistivities and thicknesses
    DC1dModelling(nlayers, ab2, mn2, verbose)
    DC1dModelling(nlayers, am, an, bm, bn, verbose)
    DC1dModelling(nlayers, dataContainer, verbose)
*/
class DLLEXPORT DC1dModelling : public ModellingBase {
public:
    /*! General constructor using AM, AN, BM, BN, distances
     * (as stored internally). */
    DC1dModelling(size_t nlayers, const RVector & am, const RVector & an,
                  const RVector & bm, const RVector & bn, bool verbose=false);

    /*! constructor for classical Schlumberger sounding */
    DC1dModelling(size_t nlayers, const RVector & ab2, const RVector & mn2,
                  bool verbose=false);

    /*! constructor using a data container */
    DC1dModelling(size_t nlayers, DataContainer & data, bool verbose=false);

    virtual ~DC1dModelling() { }

    /*! Returns an RVector of the 1dc response for
     * model = [thickness_ 0, ..., thickness_(n-1), rho_0 .. rho_n].
     * For n = nlayers. */
    RVector response(const RVector & model);

    RVector rhoa(const RVector & rho, const RVector & thk);

    RVector kern1d(const RVector & lam, const RVector & rho, const RVector & h);

    RVector pot1d(const RVector & R, const RVector & rho, const RVector & thk);

    inline RVector getK() { return k_; }

    inline RVector geometricFactor() { return k_; }

    template < class Vec > Vec rhoaT(const Vec & rho, const RVector & thk){
        Vec tmp;
        tmp  = pot1dT<Vec>(am_, rho, thk);
        tmp -= pot1dT<Vec>(an_, rho, thk);
        tmp -= pot1dT<Vec>(bm_, rho, thk);
        tmp += pot1dT<Vec>(bn_, rho, thk);
        return tmp * k_ + rho[ 0 ];
    }
    template < class Vec > Vec kern1dT(const RVector & lam, const Vec & rho,
                                       const RVector & h){
        size_t nr = rho.size();
        size_t nl = lam.size();
        Vec z(nl, rho[ nr - 1 ]);
        Vec p(nl);
        RVector th(nl);
        for (int i = nr - 2; i >= 0; i--) {
            p = (z - rho[ i ]) / (z + rho[ i ]);
            th = tanh(lam * h[ i ]);
            z = (z + toComplex(th) * rho[ i ]) / (z * th + rho[ i ]) * rho[ i ];
        }

        Vec ehl(p * RVector(exp(-2.0 * lam * h[0])));
        return ehl / (1.0 - ehl) * rho[0] / 2.0 / PI ;
    }
    template < class Vec > Vec pot1dT(const RVector & R, const Vec & rho,
                                      const RVector & thk){
        Vec z0(R.size());
        //double rabs;
        RVector rabs(abs(R));
        for (size_t i = 0; i < R.size(); i++) {
            //rabs = std::fabs(R[ i ]);
            z0[ i ] = sum(myw_ * kern1dT<Vec>(myx_ / rabs[i],
                                              rho, thk) * 2.0) / rabs[i];
        }
        return z0;
    }

    RVector createDefaultStartModel();

protected:

    /*! init myw and myx */
    void init_();

    void postprocess_();

    size_t nlayers_;
    double meanrhoa_;
    RVector am_;
    RVector an_;
    RVector bm_;
    RVector bn_;
    RVector k_;
    RVector tmp_;

    RVector myx_;
    RVector myw_;
};

/*! DC (direct current) 1D modelling for complex resistivity */
class DLLEXPORT DC1dModellingC : public DC1dModelling {
public:
    /*! General constructor using AM, AN, BM, BN, distances
     * (as stored internally). */
    DC1dModellingC(size_t nlayers,
                   const RVector & am, const RVector & an,
                   const RVector & bm, const RVector & bn, bool verbose=false);

    /*! Constructor for classical Schlumberger sounding */
    DC1dModellingC(size_t nlayers, const RVector & ab2, const RVector & mn2,
                   bool verbose=false);

    virtual ~DC1dModellingC() { }

    /*! Return [|rhoa|, +phi(rad)] for [thicks, res, phi(rad)]*/
    RVector response(const RVector & model);
};

/*! DC1dRhoModelling - Variant of DC 1D modelling with fixed parameterization
 * (occam inversion) */
/*! DC1dRhoModelling(mesh, dataContainer, thicknesses, verbose) */
class DLLEXPORT DC1dRhoModelling : public DC1dModelling {
public:
    DC1dRhoModelling(const RVector & thk, const RVector & am,
                     const RVector & an, const RVector & bm,
                     const RVector & bn, bool verbose=false)
            : DC1dModelling(thk.size(), am, an, bm, bn, verbose), thk_(thk) {
                setMesh(createMesh1D(thk.size() + 1, 1));
            }

    DC1dRhoModelling(const RVector & thk, const RVector & ab2,
                     const RVector & mn2, bool verbose = false)
            : DC1dModelling(thk.size(), ab2, mn2, verbose), thk_(thk) {
                setMesh(createMesh1D(thk.size() + 1, 1));
            }

    DC1dRhoModelling(const RVector & thk, DataContainer & data,
                     bool verbose=false)
            : DC1dModelling(thk.size(), data, verbose), thk_(thk) {
                setMesh(createMesh1D(thk.size() + 1, 1));
            }

    virtual ~DC1dRhoModelling() { }

    RVector response(const RVector & rho) {  return rhoa(rho, thk_); }

    RVector createDefaultStartModel() {
        return RVector(thk_.size() + 1, meanrhoa_);
    }

protected:
    RVector thk_;
};

} // namespace GIMLI{

#endif // _GIMLI_DC1DMODELLING__H
