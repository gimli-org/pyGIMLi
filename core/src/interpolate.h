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

#ifndef _GIMLI_INTERPOLATE__H
#define _GIMLI_INTERPOLATE__H

#include "gimli.h"
#include "matrix.h"
#include <vector>

namespace GIMLI{

/*! Utility function for interpolation. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const RVector & inVec,
                           const R3Vector & destPos, RVector & outVec,
                           bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation. */
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const R3Vector & destPos,
                              bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const RVector & inVec,
                           const Mesh & destMesh, RVector & outVec,
                           bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation. Read in data from fileName and add the
 interpolated data into the destination mesh. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const std::string & fileName,
                           Mesh & destMesh, bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const RVector & x, const RVector & y,
                              const RVector & z, bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const RVector & x, const RVector & y,
                              bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation.*/
DLLEXPORT RVector interpolate(const Mesh & srcMesh, const RVector & inVec,
                              const RVector & x,
                              bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation.
 * Interpolate all export data from srcMesh to the destination mesh.
 * Point and Cell based. */
DLLEXPORT void interpolate(const Mesh & srcMesh, Mesh & destMesh,
                           bool verbose=false, double fillValue=0.0);

/*!Interpolate a given input data regarding the mesh srcMesh to a set of positions
 * and write the interpolated data to outMat. outMat will resized if necessary.
 * Each data vector in inMat have to correspond to mesh.nodeCount().
 * If data length is mesh.cellCount() \ref cellDataToPointData will performed.
 * The interpolation rule depend on the shape functions of mesh cells.
 * Several utility or shortcut functions are defined. */
DLLEXPORT void interpolate(const Mesh & srcMesh, const RMatrix & inMat,
                           const R3Vector & destPos, RMatrix & outMat,
                           bool verbose=false, double fillValue=0.0);

/*! Utility function for interpolation.
 * Interpolate the z-coordinate from mesh to the z-coordinate of the
 * destination mesh.*/
DLLEXPORT void interpolateSurface(const Mesh & srcMesh, Mesh & destMesh,
                                  bool verbose=false, double fillValue=0);

/*! Utility function. Convert cell data to point data with the
 * corresponding the cell interpolation function */
DLLEXPORT RVector cellDataToPointData(const Mesh & mesh,
                                      const RVector & cellData);

DLLEXPORT void triangleMesh_(const Mesh & mesh, Mesh & tmpMesh);

/*! Temporary until implemented using Interpolator matrix.*/
DLLEXPORT void interpolateGradients(const Mesh & mesh,
                          const stdVectorR3Vector & qp,
                          const RVector & u,
                          stdVectorR3Vector & grad,
                          Index dim=3);
/*! Temporary until implemented using Interpolator matrix.*/
DLLEXPORT void interpolateGradients(const Mesh & mesh,
                          const stdVectorR3Vector & qp,
                          const R3Vector & u,
                          stdVectorMatrixVector & grad,
                          Index dim=3);

} // namespace GIMLI


#endif // _GIMLI_INTERPOLATE__H
