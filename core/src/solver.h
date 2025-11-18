/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *   Thomas Günther thomas@resistivity.net                                    *
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

#ifndef _GIMLI_SOLVER__H
#define _GIMLI_SOLVER__H

#include "gimli.h"
#include "matrix.h"
#include "vectortemplates.h"

namespace GIMLI{

typedef RVector Vec;

DLLEXPORT int solveCGLSCDWWhtrans(const MatrixBase & S, const MatrixBase & C,
                        const Vec & dWeight, const Vec & b, Vec & x,
                        const Vec & wc, const Vec & wm,
                        const Vec & tm, const Vec & td,
                        double lambda, const Vec & roughness,
                        int maxIter=200, double tol=-1.0,
                        bool verbose=false);

DLLEXPORT int solveCGLSCDWWtrans(const MatrixBase & S, const MatrixBase & C,
                       const Vec & dWeight,  const Vec & b, Vec & x,
                       const Vec & wc, const Vec & mc, const Vec & tm,
                       const Vec & td, double lambda, const Vec & deltaX,
                       int maxIter=200, bool verbose=false);

} //namespace _GIMLI_SOLVER__H

#endif
