/******************************************************************************
 *   Copyright (C) 2005-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_TRIANGLEWRAPPER__H
#define _GIMLI_TRIANGLEWRAPPER__H

#include "gimli.h"

#include <vector>
#include <string>

struct triangulateio;

namespace GIMLI{

class DLLEXPORT TriangleWrapper{
public:
    /*! Construct a PLC from nodes and edges in inMesh and initialize
     * the trianglewrapper. Add region and hole marker later.*/
    TriangleWrapper(const Mesh & inMesh);

    /*! Constructor initialize with input mesh
     * (nodes and edges where used for input and also region marker,
     * if added to the mesh) and generate the 2d mesh with corresponding
     * triSwitches */
    TriangleWrapper(const Mesh & inMesh, Mesh & outMesh,
                    const std::string & triSwitches);

    /*! Default destructur */
    virtual ~TriangleWrapper();

    /*! Set the triangle commandline switches */
    void setSwitches(const std::string & s);

    /*! Return the triangle switches. */
    inline const std::string & switches() const { return switches_; }

    /*! Generate the mesh and store in mesh. */
    void generate(Mesh & mesh);

    /*! Generate and return the new mesh. */
    Mesh generate();

protected:
    /*! For internal use only. */
    void init_();

    /*! For internal use only. */
    void transformTriangleToMesh_(const triangulateio & trimesh, Mesh & mesh);

    /*! For internal use only. Only Edges and nodes(x,y) from mesh are used. */
    void transformMeshToTriangle_(const Mesh & mesh, triangulateio & trimesh);

    /*! For internal use only. */
    void allocateOutMemory_();

    /*! For internal use only. */
    void freeMemory_();

    struct triangulateio * mesh_input_;
    struct triangulateio * mesh_output_;
    struct triangulateio * mesh_voronoi_output_;

    std::string switches_;

    const Mesh * inMesh_;

};

}  // namespace GIMLI

#endif // _GIMLI_TRIANGLEWRAPPER__H
