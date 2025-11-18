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

#ifndef _BERT_BERT_JACOBIAN__H
#define _BERT_BERT_JACOBIAN__H

#include "bert.h"

#include <vector.h>

namespace GIMLI{

/* Recommended. Need more mem than createSensitivityRow, but is extremly faster.*/
// DLLEXPORT void createSensitivityCol(RMatrix & S, const Mesh & mesh,
//                                      const Container & data,
// 				     const RMatrix & pots, bool verbose);
//
// DLLEXPORT void createSensitivityCol(RMatrix & S, const Mesh & mesh,
//                                      const DataContainer & data,
// 				     const RMatrix & pots, int nThreads, bool verbose);
//
// DLLEXPORT void createSensitivityCol(RMatrix & S, const Mesh & mesh,
//                                      const DataContainer & data,
//                                      const RMatrix & pots,
//                                      const std::map< long, uint > & currPatternIdx,
// 				     int nThreads, bool verbose);

DLLEXPORT void createSensitivityCol(RMatrix & S,
                                    const Mesh & mesh,
                                    const DataContainerERT & data,
                                    const RMatrix & pots,
                                    const RVector & weights,
                                    const RVector & k,
                                    std::vector < std::pair < Index, Index > > & matrixClusterIds,
                                    uint nThreads, bool verbose);

DLLEXPORT void createSensitivityCol(CMatrix & S,
                                    const Mesh & mesh,
                                    const DataContainerERT & data,
                                    const CMatrix & pots,
                                    const RVector & weights,
                                    const RVector & k,
                                    std::vector < std::pair < Index, Index > > & matrixClusterIds,
                                    uint nThreads, bool verbose);

DLLEXPORT void sensitivityDCFEMSingle(const std::vector < Cell * > & para,
                                      const RVector & p1, const RVector & p2,
                                      RVector & sens, bool verbose);

/*! log10 scale of sensitivity data for visualisation.
 mesh need to be a parameter mesh, data is sensitivity matrix row of
 length nModel. Returning vector have the length of mesh.cellSize()
 and can be viewed directly. */
DLLEXPORT RVector prepExportSensitivityData(const Mesh & mesh,
                                            const RVector & data,
                                            double logdrop=1e-3);

DLLEXPORT void exportSensitivityVTK(const std::string & fileName,
                                    const Mesh & mesh,
                                    const RVector & data,
                                    double logdrop=1e-3);

// DLLEXPORT void exportSensMatrixDC(const std::string & filename,
//                                   const Mesh & mesh,
//                                   const RMatrix & S);
DLLEXPORT void exportSensMatrixDC(const std::string & filename,
                                  const Mesh & mesh,
                                  const RMatrix & S,
                                  const IVector & idx,
                                  double logdrop=1e-3);

DLLEXPORT RVector coverageDC(const RMatrix & S);

DLLEXPORT RVector coverageDCtrans(const MatrixBase & S,
                                  const RVector & dd,
                                  const RVector & mm);

DLLEXPORT RVector createCoverage(const MatrixBase & S, const Mesh & mesh);

DLLEXPORT RVector createCoverage(const MatrixBase & S, const Mesh & mesh,
                                 const RVector & response, const RVector & model);

} // namespace BERT

#endif //_BERT_BERT_JACOBIAN__H
