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

#include "interpolate.h"

#include "meshentities.h"
#include "mesh.h"
#include "node.h"
#include "shape.h"

namespace GIMLI{

void interpolate(const Mesh & mesh, const RMatrix & vData,
                 const R3Vector & ipos, RMatrix & iData,
                 bool verbose, double fillValue){ ALLOW_PYTHON_THREADS

    R3Vector pos(ipos);

    if (mesh.dim() == 2){
        if ((zVari(pos) || max(abs(z(pos))) > 0.) &&
            (!yVari(pos) && max(abs(y(pos))) < 1e-8)) {
            if (verbose)
                std::cout << "Warning! swap YZ coordinates for query "
                            "positions to meet mesh dimensions." << std::endl;
            swapYZ(pos);
        }
    }

    if (iData.rows() != vData.rows()){
        iData.resize(vData.rows(), pos.size());
    }

    std::vector < Cell * > cells(pos.size());
    size_t count = 0;

    Cell * c = 0;
    for (uint i = 0; i < pos.size(); i ++) {
//         __MS(pos[i])
        c = mesh.findCell(pos[i], count, false);
//         if (!c) {
//             __MS(pos[i])
// //             c = mesh.findCell(pos[i], count, true);
// //             if (!c) exit(0);
//         }
        cells[i] = c;
//         if (verbose) std::cout << "\r" << i + 1 << " \t/ " << pos.size()
//                                << "\t searched: " << count << std::endl;
    }
    if (verbose) std::cout << std::endl;

    for (uint i = 0; i < vData.rows(); i ++) {
        if (verbose) std::cout << "\r" << i + 1 << " \t/ " << vData.rows();

        RVector data;

        if (vData[i].size() != 0){

            if (vData[i].size() == mesh.nodeCount()){
                data = vData[i];
            } else if (vData[i].size() == mesh.cellCount()){
                data = cellDataToPointData(mesh, vData[i]);
            } else {
                throwLengthError(
                                 WHERE_AM_I +
                                 " data.size not nodeCount and cellCount " +
                                 str(vData.size()) + " x " +
                                 str(vData[i].size()) + " != " +
                                 str(mesh.nodeCount()) + " != " +
                                 str(mesh.cellCount()));
            }

            for (uint j = 0; j < pos.size(); j ++) {
                if (cells[j]){
                    iData[i][j] = cells[j]->pot(pos[j], data);

//              this check is obsolete if the findCell (from above) is working properly
//                     if (cells[j]->shape().isInside(pos[j])){
//                         iData[i][j] = cells[j]->pot(pos[j], data);
//                     } else {
//                         std::cout << WHERE_AM_I << " here is somethng wrong" << std::endl;
//                     }

//                    std::cout << j << " " << iData[i][j] << std::endl;
                    //** return cell data
//                    iData[i][j] = vData[i][cells[j]->id()];
                } else {
                    iData[i][j] = fillValue;
                    //std::cout << "Cant find cell for " << pos[j]<< std::endl;
//                     for (uint i = 0; i < mesh.cellCount(); i ++){
//                     	if (mesh.cell(i).shape().isInside(pos[j], true)){
//                         	std::cout << mesh.cell(i) << std::endl;
//                         }
//                     }
//                     exit(0);
                }
            }
        } // if vData.size() != 0
    }  // for each in vData
    if (verbose) std::cout << std::endl;
}

void interpolate(const Mesh & mesh, const RVector & data,
                 const Mesh & pos, RVector & iData, bool verbose, double fillValue){

    RMatrix vData; vData.push_back(data);
    RMatrix viData;
    interpolate(mesh, vData, pos.positions(), viData, verbose, fillValue);
    iData = viData[0];
}

RVector interpolate(const Mesh & mesh, const RVector & data,
                    const R3Vector & pos, bool verbose, double fillValue){

    RMatrix vData; vData.push_back(data);
    RMatrix viData;
    interpolate(mesh, vData, pos, viData, verbose, fillValue);
    return viData[0];
}

void interpolate(const Mesh & mesh, const std::string & dataName, Mesh & pos,
                 bool verbose, double fillValue){
    RMatrix vData; vData.push_back(mesh.data(dataName));
    RMatrix viData;
    interpolate(mesh, vData, pos.positions(), viData, verbose, fillValue);
    pos.addData(dataName, viData[0]);
}

void interpolate(const Mesh & mesh, const RVector & data,
                 const R3Vector & pos,
                 RVector & iData, bool verbose, double fillValue){
    RMatrix vData; vData.push_back(data);
    RMatrix viData;
    interpolate(mesh, vData, pos, viData, verbose, fillValue);
    iData = viData[0];
}

RVector interpolate(const Mesh & mesh, const RVector & data,
                    const RVector & x, const RVector & y,
                    const RVector & z, bool verbose, double fillValue){

    if (z.size() == 0) return interpolate(mesh, data, x, y, verbose, fillValue);
    if (z.size() == 0 && y.size() == 0) return interpolate(mesh, data, x, verbose, fillValue);

    if (x.size() != y.size() || x.size() != z.size()) {
        throwLengthError(" x.size invalid y.size invalid z.size() "
                + str(x.size()) + " != " + str(y.size()) + " != " + str(z.size()));
    }

    std::vector < RVector3 > pos(x.size());
    for (uint i = 0; i < x.size(); i ++) pos[i] = RVector3(x[i], y[i], z[i]);

    RVector iData;
    interpolate(mesh, data, pos, iData, verbose, fillValue);
    return iData;
}

RVector interpolate(const Mesh & mesh, const RVector & data,
                    const RVector & x, const RVector & y,
                    bool verbose, double fillValue){
    if (y.size() == 0) return interpolate(mesh, data, x, verbose, fillValue);
    return interpolate(mesh, data, x, y, RVector(x.size(), 0.0), verbose, fillValue);
}

RVector interpolate(const Mesh & mesh, const RVector & data,
                    const RVector & x, bool verbose, double fillValue){
    return interpolate(mesh, data, x,
                       RVector(x.size(), 0.0), RVector(x.size(), 0.0), fillValue);
}

void interpolate(const Mesh & mesh, Mesh & qmesh, bool verbose, double fillValue){
    RMatrix cellData;
    RMatrix nodeData;
    std::vector< std::string > cellDataNames;
    std::vector< std::string > nodeDataNames;

    for (std::map< std::string, RVector >::const_iterator it = mesh.dataMap().begin();
          it != mesh.dataMap().end(); it ++){

        if (it->second.size() == mesh.nodeCount()){
            if (verbose) std::cout << " interpolate node data: " << it->first << std::endl;
            nodeData.push_back(it->second);
            nodeDataNames.push_back(it->first);
        } else if (it->second.size() == mesh.cellCount()){
            if (verbose) std::cout << " interpolate cell data: " << it->first << std::endl;
            cellData.push_back(it->second);
            cellDataNames.push_back(it->first);
        } else {
            if (verbose) std::cout << " omit unknown data: " << it->first << " " <<
                it->second.size()<< std::endl;
        }
    }

    if (cellData.rows() > 0){
        RMatrix qCellData;
        interpolate(mesh, cellData, qmesh.cellCenter(), qCellData, verbose, fillValue);
        for (uint i= 0; i < cellData.rows(); i ++){
            qmesh.addData(cellDataNames[i], qCellData[i]);
        }
    }
    if (nodeData.rows() > 0){
        RMatrix qNodeData;
        interpolate(mesh, nodeData, qmesh.positions(), qNodeData, verbose) ;
        for (uint i= 0; i < nodeData.rows(); i ++){
            qmesh.addData(nodeDataNames[i], qNodeData[i]);
        }
    }
}

void interpolateSurface(const Mesh & mesh, Mesh & qmesh, bool verbose, double fillValue){
    RVector z(mesh.nodeCount());
    for (uint i = 0; i < z.size(); i ++) z[i] = mesh.node(i).pos()[2];
    RVector qz(qmesh.nodeCount());
    interpolate(mesh, z, qmesh, qz, verbose, fillValue);
    for (uint i = 0; i < qz.size(); i ++) qmesh.node(i).pos()[2] = qz[i];
}

void triangleMesh_(const Mesh & mesh, Mesh & tmpMesh){

    for (uint i = 0; i < mesh.nodeCount(); i ++) tmpMesh.createNode(mesh.node(i));
    for (uint i = 0; i < mesh.cellCount(); i ++) {
        switch(mesh.cell(i).rtti()){
            case MESH_TRIANGLE_RTTI:
                tmpMesh.createCell(mesh.cell(i));
            	tmpMesh.cell(tmpMesh.cellCount() -1).setId(mesh.cell(i).id());
            	break;
            case MESH_QUADRANGLE_RTTI:
                tmpMesh.createTriangle(tmpMesh.node(mesh.cell(i).node(0).id()),
                                        tmpMesh.node(mesh.cell(i).node(1).id()),
                                        tmpMesh.node(mesh.cell(i).node(2).id()),
                                        mesh.cell(i).marker());
                tmpMesh.createTriangle(tmpMesh.node(mesh.cell(i).node(0).id()),
                                        tmpMesh.node(mesh.cell(i).node(2).id()),
                                        tmpMesh.node(mesh.cell(i).node(3).id()),
                                        mesh.cell(i).marker());
                break;
          default:
            std::cerr << WHERE_AM_I << " nothing known to cell.rtti " << mesh.cell(i).rtti() << std::endl;
        }
    }
    tmpMesh.createNeighborInfos(true);
}


RVector cellDataToPointData(const Mesh & mesh, const RVector & cellData){
    if (cellData.size() != mesh.cellCount()){
        throwLengthError(" vector size invalid mesh.cellCount "
                        + str(mesh.cellCount()) + " != " + str(cellData.size()));
    }

    RVector ret(mesh.nodeCount());

    std::set < Cell * > cset;
    for (uint i = 0; i < mesh.nodeCount(); i ++){
        cset = mesh.node(i).cellSet();
        for (std::set < Cell * >::iterator it = cset.begin(); it != cset.end(); it ++){
            ret[i] += cellData[(*it)->id()];
        }
        ret[i] /= cset.size();
    }

    return ret;


}

// double interpolate(const RVector3 & queryPos, const MeshEntity & entity, const RVector & sol){
//   return entity->interpolate(queryPos, sol);
// }

} // namespace GIMLI
