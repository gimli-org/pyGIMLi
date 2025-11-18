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

#ifndef _GIMLI_EM1DMODELLING__H
#define _GIMLI_EM1DMODELLING__H

#include "gimli.h"
#include "mesh.h"
#include "meshgenerators.h"
#include "modellingbase.h"

namespace GIMLI{
/*! this file holds different electromagnetic forward operators for 1D discretizations */

//! Magnetotellurics (MT) 1D modelling
/*! Classical MT 1D block model yielding rhoa&phase for given resistivities&thicknesses */
class DLLEXPORT MT1dModelling : public ModellingBase {
public:
    MT1dModelling(const RVector & periods, size_t nlay, bool verbose = false)
        : ModellingBase(verbose), periods_(periods), nlay_(nlay) {
        setMesh(createMesh1DBlock(nlay));
    }

//    MT1dModelling(Mesh & mesh, const RVector & periods, size_t nlay, bool verbose = false)
//        : ModellingBase(mesh, verbose), periods_(periods), nlay_(nlay) {
//    }

    virtual ~MT1dModelling() { }

    /*! different sub-forward operators for alternate use */
    virtual RVector rhoaphi(const RVector & rho, const RVector & thk); //! app. res. and phase

    virtual RVector rhoa(const RVector & rho, const RVector & thk);    //! only app. res.

    virtual RVector rhoa(const RVector & model);
    //! app. res. for thk/res vector


    /*! the actual (full) forward operator returning app.res.+phase for thickness+resistivity */
    virtual RVector response(const RVector & model);

protected:
    RVector periods_;
    size_t nlay_;
};

/*! Magnetotellurics (MT) 1D modelling */
/*! variant with fixed thickness vector (for smooth (Occam) inversion) */
/*! MT1dRhoModelling(RVector periods, thk [, verbose]) */
class DLLEXPORT MT1dRhoModelling : public MT1dModelling {
public:
//    MT1dRhoModelling(Mesh & mesh, const RVector & periods, const RVector & thk, bool verbose = false)
//        : MT1dModelling(mesh, periods, thk.size(), verbose), thk_(thk) {
//    }

    MT1dRhoModelling(RVector & periods, RVector & thk, bool verbose = false)
        : MT1dModelling(periods, thk.size(), verbose), thk_(thk) {
            setMesh(createMesh1D(thk.size() + 1, 1));
        }

    virtual ~MT1dRhoModelling() { }

    virtual RVector response(const RVector & rho) { return rhoaphi(rho, thk_); }

    virtual RVector rhoa(const RVector & rho) { return MT1dModelling::rhoa(rho, thk_); }

protected:
    RVector thk_;
};

/*! Frequency domain electromagnetic (FDEM) sounding using a block discretization */
/*! FDEM1dModelling(nlay, RVector freq, coilspacing [, elevation, verbose ]) */
class DLLEXPORT FDEM1dModelling : public ModellingBase {
public:
    //! default constructor creating a block model
    FDEM1dModelling(size_t nlay, const RVector & freqs,
                    const RVector & coilspacing,
                    double z=0.0, bool verbose=false);

    FDEM1dModelling(size_t nlay, const RVector & freqs,
                    double coilspacing,
                    double z=0.0, bool verbose=false);

    virtual ~FDEM1dModelling() { }

    void init();

    virtual RVector response(const RVector & model);

    /*! Return reference to the freeAirSolution.
     * This will be automatic filled by and init() call during
     * instance construction. */
    const RVector & freeAirSolution() const { return freeAirSolution_; }

    RVector calc(const RVector & rho, const RVector & thk);


protected:
    size_t nlay_;
    RVector freqs_;
    RVector coilspacing_;
    double zs_, ze_; // transmitter&receiver heights (minus)
    size_t nfr_;
    RVector freeAirSolution_;
};

//class MaxMinModelling:FDEMModelling
//110,220,440,880,1760,3520,7040,14080,28160,56320

/*! Frequency Domain EM 1d modelling with predefined thickness vector */
/*! FDEM1dRhoModelling(RVector thk, freq, coilspacing[, double elevation, verbose]) */
/*! FDEM1dRhoModelling(RVector thk, freq, double coilspacing[, elevation, verbose]) */
class DLLEXPORT FDEM1dRhoModelling : public FDEM1dModelling {
public:
    //! default constructor creating a block model
    FDEM1dRhoModelling(const RVector & thk, const RVector & freq, const RVector & coilspacing, double z=0.0, bool verbose=false)
        : FDEM1dModelling(thk.size(), freq, coilspacing, z, verbose), thk_(thk) {
            setMesh(createMesh1D(thk.size() + 1, 1));
        }

    FDEM1dRhoModelling(const RVector & thk, const RVector & freq, double coilspacing, double z=0.0, bool verbose=false)
        : FDEM1dModelling(thk.size(), freq, coilspacing, z, verbose), thk_(thk) {
            setMesh(createMesh1D(thk.size() + 1, 1));
        }

    virtual ~FDEM1dRhoModelling() { }

    RVector response(const RVector & model){ return calc(model, thk_); }

protected:
    RVector thk_;
};

/*! Magnetic Resonance Sounding (MRS) modelling */
/*! classical variant using a fixed parameterization */
/*! MRSModelling([mesh,] RMatrix KR, KI [, verbose]) */
class DLLEXPORT MRSModelling : public ModellingBase {
public:
    //! constructor using a predefined mesh and real/imag matrix
    MRSModelling(Mesh & mesh, RMatrix & KR, RMatrix & KI, bool verbose = false) :
            ModellingBase(mesh, verbose), KR_(&KR), KI_(&KI) { }
    //! constructor with real/imag matrix only
    MRSModelling(RMatrix & KR, RMatrix & KI, bool verbose = false) :
            ModellingBase(verbose), KR_(&KR), KI_(&KI) { }
    //! destructor
    virtual ~MRSModelling() { }

    /*! return response voltage for a given water content vector */
    RVector response(const RVector & model);

    /*! create jacobian matrix (AmplitudeJacobian) */
    void createJacobian(const RVector & model);

protected:
    RMatrix *KR_, *KI_;
};

/*! Magnetic Resonance Sounding (MRS) modelling
    non-classical variant using a block model using mapping to classical varian
    MRS1dBlockModelling(int nlay, RMatrix KR, KI, RVector zvec [, verbose]) */
class DLLEXPORT MRS1dBlockModelling : public MRSModelling{
public:
    /*! constructor using mesh, real/imag matrix and zvector */
    MRS1dBlockModelling(Mesh & mesh, RMatrix & KR, RMatrix & KI, RVector & zvec, bool verbose = false)
        : MRSModelling(mesh, KR, KI, verbose), nlay_((mesh.cellCount() + 1) / 2),
                        nvec_(KR.cols()),zvec_(zvec) {
    }
    /*! constructor using real/imag matrix and zvector */
    MRS1dBlockModelling(int nlay, RMatrix & KR, RMatrix & KI, RVector & zvec, bool verbose = false)
        : MRSModelling(KR, KI, verbose), nlay_(nlay), nvec_(KR.cols()), zvec_(zvec) {
            setMesh(createMesh1DBlock(nlay));
        }

    virtual ~MRS1dBlockModelling() { }
    /*! return voltage for a given block model vector */
    RVector response(const RVector & model);

    /*! use default (brute-force) jacobian generator (override the smooth one) */
    void createJacobian(const RVector & model){
        ModellingBase::createJacobian(model);
    }

protected:
    size_t nlay_;
    size_t nvec_;
    RVector zvec_;
};



} // namespace GIMLI{

#endif // _GIMLI_EM1D__H
