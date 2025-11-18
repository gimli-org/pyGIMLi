/******************************************************************************
 *   Copyright (C) 2006-2024 by the resistivity.net development team          *
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

#ifndef _BERT_BERT_MISC__H
#define _BERT_BERT_MISC__H

#include "bert.h"

#include <vector.h>
#include <pos.h>

namespace GIMLI{

DLLEXPORT double exactDCSolution(const RVector3 & pot, const RVector3 & src,
                                 double k, double surfaceZ, double fallback);

DLLEXPORT double exactDCSolution(const RVector3 & pot, const RVector3 & src);

DLLEXPORT RVector exactDCSolution(const Mesh & mesh, const RVector3 & src,
                                   double k, double surfaceZ=0.0);
DLLEXPORT RVector exactDCSolution(const Mesh & mesh,
                                  const Node * nA, const Node * nB,
                                   double k, double surfaceZ=0.0);

/*! Calculate the analytical solution as well as an aproxmiate value
 * for the singular position.*/
DLLEXPORT RVector exactDCSolution(const Mesh & mesh,
                                  const ElectrodeShape * elec,
                                  double k, double surfaceZ=0.0,
                                  bool setSingValue=true);

DLLEXPORT RVector exactDCSolution(const Mesh & mesh, int aID, int bID, double k,
                                  double surfaceZ=0.0);
DLLEXPORT RVector exactDCSolution(const Mesh & mesh, int aID, double k=0.0,
                                  double surfaceZ=0.0);


/*! Helper function to calculate configuration factors for a
 * given \ref DataContainerERT */
DLLEXPORT RVector geometricFactors(const DataContainerERT & data, int dim=3,
                                   bool forceFlatEarth=false);

/*! DEPRECATED due to wrong typo. */
inline RVector geometricFactor(const DataContainerERT & data, int dim=3,
                                  bool forceFlatEarth=false){
    __MS("Deprecated please use 'geometricFactors'")
    return geometricFactors(data, dim, forceFlatEarth);
}

DLLEXPORT void initKWaveList(double rMin, double rMax,
                             int nGauLegendre, int nGauLaguerre,
                             RVector & kValues, RVector & weights);

DLLEXPORT void initKWaveList(const Mesh & mesh,
                             RVector & kValues, RVector & weights,
                             const R3Vector & sources,
                             bool verbose=false);

DLLEXPORT void initKWaveList(const Mesh & mesh,
                             RVector & kValues, RVector & weights,
                             bool verbose=false);

DLLEXPORT int countKWave(const Mesh & mesh);

DLLEXPORT void DCErrorEstimation(DataContainerERT & data,
                                 double errPerc=3.0, double errVolt=100e-6,
                                 double defaultCurrent=100e-3,
                                 bool verbose=false);

DLLEXPORT double DCParaDepth(const DataContainerERT & data);

/*! Set the boundary marker at the outer boundarys.
MARKER_BOUND_HOMOGEN_NEUMANN at the surface (boundary.center == zmax) and otherwise MARKER_BOUND_MIXED. */
DLLEXPORT void setDefaultBERTBoundaryConditions(Mesh & mesh);

/*! Set the boundary marker to MARKER_BOUND_HOMOGEN_NEUMANN at the outer boundarys. */
DLLEXPORT void setAllNeumannBoundaryConditions(Mesh & mesh);

/*! log10 scale of potential data for visualisation. */
DLLEXPORT RVector prepExportPotentialData(const RVector & data, double logdrop=1e-6);

} // namespace BERT

#endif // _BERT_BERT_MISC__H


