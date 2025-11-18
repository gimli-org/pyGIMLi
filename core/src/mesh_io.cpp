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

#include "mesh.h"
#include "node.h"
#include "matrix.h"
#include "pos.h"
#include "vectortemplates.h"

#include <map>
#include <fstream>

namespace GIMLI{

void Mesh::load(const std::string & fbody, bool createNeighbors, IOFormat format){
    if (fbody.find(".mod") != std::string::npos){
        importMod(fbody);
    } else if (fbody.find(".vtk") != std::string::npos){
        importVTK(fbody);
    } else if (fbody.find(".vtu") != std::string::npos){
        importVTU(fbody);
    } else if (format == Binary || fbody.find(MESHBINSUFFIX) != std::string::npos){
        try {
             loadBinaryV2(fbody);
        } catch(std::exception & e){
            //std::cout << "Failed to loadBinary " << e.what() << std::endl;
            //std::cout << "try load bms.v2" << std::endl;
            loadBinary(fbody);
        }
    } else {
        loadAscii(fbody);
    }
    if (createNeighbors) createNeighborInfos();
}

void Mesh::loadAscii(const std::string & fbody){
    THROW_TO_IMPL
}

int Mesh::save(const std::string & fbody, IOFormat format) const {
    if (format == Binary || fbody.find(MESHBINSUFFIX) != std::string::npos) {
        saveBinaryV2(fbody);
    } else return saveAscii(fbody);
    return 1;
}

int Mesh::saveAscii(const std::string & fbody) const {
  //** save nodes
  std::fstream file; if (!openOutFile(fbody + ".n", & file)) {  return 0; }

  file.precision(14);
  for (uint i = 0; i < nodeCount(); i ++){
    for (uint j = 0; j < dim(); j ++) file << node(i).pos()[j] << "\t";
    file << node(i).marker() << std::endl;
  }
  file.close();

  //** save cells;
  if (!openOutFile(fbody + ".e", & file)) { return 0; }

  for (int i = 0, imax = cellCount(); i < imax; i++){
    for (int j = 0, jmax = cell(i).nodeCount(); j < jmax; j++){
      file << cell(i).node(j).id() << "\t" ;
    }
    file << cell(i).marker() << std::endl;
  }
  file.close();

  //** save boundarys;
  if (!openOutFile(fbody + ".s", & file)) { return 0; }

  for (int i = 0, imax = boundaryCount(); i < imax; i++){
    for (int j = 0, jmax = boundary(i).nodeCount(); j < jmax; j++){
      file << boundary(i).node(j).id() << "\t" ;
    }
    file << "-33\t-33\t" << boundary(i).marker() << std::endl;
  }

  file.close();
  return 1;
}

int Mesh::saveBinary(const std::string & fbody) const {
    std::string fileName(fbody.substr(0, fbody.rfind(MESHBINSUFFIX)) + MESHBINSUFFIX);

    FILE *file;
    file = fopen(fileName.c_str(), "w+b");

    int dimension = dim();
    if (!fwrite(&dimension, sizeof(int), 1, file)){
    }

    //! write vertex dummy-infos;
    int dummy[127]; memset(dummy, 0, 127 * sizeof(int));
    dummy[0] = int(this->isGeometry_);
    if (!fwrite(dummy, sizeof(int), 127, file)) {
    }

    int nVerts = nodeCount();
    if (!fwrite(&nVerts, sizeof(int), 1, file)){
    }

    double * coord = new double[dimension * nVerts];
    for (int i = 0; i < nVerts; i ++){
        for (int j = 0; j < dimension; j ++){
            coord[i * dimension + j] = node(i).pos()[j];
        }
    }

    if (!fwrite(coord, sizeof(double), dimension * nVerts, file)){

    }
    delete [] coord;

    int * marker = new int[nVerts];
    for (int i = 0; i < nVerts; i ++) marker[i] = node(i).marker();

    if (!fwrite(marker, sizeof(int), nVerts, file)){

    }
    delete [] marker;

    //! write cell dummy-infos
    if (!fwrite(dummy, sizeof(int), 127, file)){

    }

    int nCells = cellCount();
    if (!fwrite(&nCells, sizeof(int), 1, file)){

    }

    int * cellVerts = new int[nCells];
    int nCellIdx = 0;

    for (int i = 0; i < nCells; i ++){
        cellVerts[i] = cell(i).nodeCount();
        nCellIdx += cellVerts[i];
    }
    if (!fwrite(cellVerts, sizeof(int), nCells, file)){

    }

    int * idx = new int[nCellIdx];

    int count = 0;
    for (int i = 0; i < nCells; i ++){
        for (int j = 0; j < cellVerts[i]; j ++){
            idx[count] = cell(i).node(j).id();
            count ++;
        }
    }

    if (!fwrite(idx, sizeof(int), nCellIdx, file)){

    }
    delete [] idx;
    delete [] cellVerts;

    double * attribute = new double[nCells];
    for (int i = 0; i < nCells; i ++) attribute[i] = cell(i).marker();

    if (!fwrite(attribute, sizeof(double), nCells, file)){
    }

    delete [] attribute;

  //** write boundary dummy-infos
    if (!fwrite(dummy, sizeof(int), 127, file)){

    }

    int nBounds = boundaryCount();
    if (!fwrite(&nBounds, sizeof(int), 1, file)){

    }

    int * boundVerts = new int[nBounds];
    int nBoundIdx = 0;
    for (int i = 0; i < nBounds; i ++){
        boundVerts[i] = boundary(i).nodeCount();
        nBoundIdx += boundVerts[i];
    }
    if (!fwrite(boundVerts, sizeof(int), nBounds, file)){

    }

    idx = new int[nBoundIdx];
    count = 0;
    for (int i = 0; i < nBounds; i ++){
        for (int j = 0; j < boundVerts[i]; j ++){
            idx[count] = boundary(i).node(j).id();
            count ++;
        }
    }

    if (!fwrite(idx, sizeof(int), nBoundIdx, file)){

    }

  delete [] idx;
  delete [] boundVerts;

  marker = new int[nBounds];
  for (int i = 0; i < nBounds; i ++) marker[i] = boundary(i).marker();
 if (!fwrite(marker, sizeof(int), nBounds, file)){}
  delete [] marker;

  idx = new int[nBounds];
  for (int i = 0; i < nBounds; i ++){
    if (boundary(i).leftCell() != NULL){
      idx[i] = boundary(i).leftCell()->id();
    } else {
      idx[i] = -1;
    }
  }
  if (!fwrite(idx, sizeof(int), nBounds, file)){}

  for (int i = 0; i < nBounds; i ++){
    if (boundary(i).rightCell() != NULL){
      idx[i] = boundary(i).rightCell()->id();
    } else {
      idx[i] = -1;
    }
  }
  if (!fwrite(idx, sizeof(int), nBounds, file)) {}
  delete [] idx;

  fclose(file);
  return 1;
}

void Mesh::loadBinary(const std::string & fbody){
//   sizeof(int) = 4 byte
//   int[1] dimension
//   int[127] dummy vertices information, dummy[0] = isGeometry_
//   int[1] nVerts, number of vertices
//   double[dimension * nVerts]; coordinates,  dimension == 2 (x, y), dimension == 3 (x, y, z)
//   int[nVerts] vertex markers
//   int[127] dummy cell information
//   int[1] nCells, number of cell
//   int[nCells] cellVerts; number of nodes for each cell
//   int[sum(cellVerts)] cellsidx
//   double[nCells] attribute, cell attributes
//   int[127] dummy boundary information
//   int[1] nBounds, number of boundarys
//   int[nBounds] boundVerts; numbers of nodes for each boundary
//   int[sum(boundVerts)] boundIdx
//   int[nBounds] boundary markers
//   int[nBounds] leftNeighbor idx (-1) if no neighbor present or info unavailable
//   int[nBounds] rightNeighbor idx (-1) if no neighbor present or info unavailable
    //std::cout << "load binary " << fbody << std::endl;
    clear();
    std::string fileName(fbody.substr(0, fbody.rfind(MESHBINSUFFIX)) + MESHBINSUFFIX);

    FILE *file; file = fopen(fileName.c_str(), "r+b");
    if (!file) {
        throwError(WHERE_AM_I + " " + fileName + ": " + strerror(errno));
    }

    int dim = 0;
    uint ret = fread(&dim, sizeof(int), 1, file);

    if ((dim !=2 && dim !=3) || (ret == 0)){
        throwError(WHERE_AM_I + " cannot determine dimension " + str(dim));
    }
    this->setDimension(dim);
    //** read vertex dummy-infos
    int dummy[127]; ret = fread(dummy, sizeof(int), 127, file);
    this->setGeometry(bool(dummy[0]));

    int nVerts; ret = fread(&nVerts, sizeof(int), 1, file);
    double * coords = new double[dimension_ * nVerts];
    ret = fread(coords, sizeof(double), dimension_ * nVerts, file);

    int * nodeMarker = new int[nVerts]; ret = fread(nodeMarker, sizeof(int), nVerts, file);
    //** read cell dummy-infos
    ret = fread(dummy, sizeof(int), 127, file);

    int nCells; ret = fread(&nCells, sizeof(int), 1, file);
    int * cellVerts = new int[nCells]; ret = fread(cellVerts, sizeof(int), nCells, file);
    int nCellIdx = 0; for (int i = 0; i < nCells; i ++) nCellIdx += cellVerts[i];
    int * cellIdx = new int[nCellIdx]; ret = fread(cellIdx, sizeof(int), nCellIdx, file);
    double * attribute = new double[nCells]; ret = fread(attribute, sizeof(double), nCells, file);

    //** read boundary dummy-infos
    ret = fread(dummy, sizeof(int), 127, file);

    int nBounds; ret = fread(&nBounds, sizeof(int), 1, file);
    int * boundVerts = new int[nBounds]; ret = fread(boundVerts, sizeof(int), nBounds, file);
    int nBoundIdx = 0; for (int i = 0; i < nBounds; i ++) nBoundIdx += boundVerts[i];
    int * boundIdx = new int[nBoundIdx]; ret = fread(boundIdx, sizeof(int), nBoundIdx, file);
    int * boundMarker = new int[nBounds]; ret = fread(boundMarker, sizeof(int), nBounds, file);
    int * left = new int[nBounds]; ret = fread(left, sizeof(int), nBounds, file);
    int * right = new int[nBounds]; ret = fread(right, sizeof(int), nBounds, file);

    //** create Nodes;
    nodeVector_.reserve(nVerts);
    for (int i = 0; i < nVerts; i ++){
        RVector3 pos;
        for (uint j = 0; j < dimension_; j ++){
            pos[j] = coords[i * dimension_ + j];
        }
        createNode(pos, nodeMarker[i]);
    }

    //** create Cells;
    int count = 0;
    std::vector < Node * > pNodeVector;

    cellVector_.reserve(nCells);
    for (int i = 0; i < nCells; i ++){
        std::vector < Node * > nodes(cellVerts[i]);
        for (uint j = 0; j < nodes.size(); j ++) {
            nodes[j] = & node(cellIdx[count + j]);
        }
        createCell(nodes);
        count += cellVerts[i];
    }

    for (int i = 0; i < nCells; i ++) {
        cell(i).setMarker((int)rint(attribute[i]));
        cell(i).setAttribute(attribute[i]);
    }

    //** create Boundaries;
    count = 0;
    boundaryVector_.reserve(nBounds);
    for (int i = 0; i < nBounds; i ++){
        std::vector < Node * > nodes(boundVerts[i]);
        for (uint j = 0; j < nodes.size(); j ++) nodes[j] = & node(boundIdx[count + j]);
        createBoundary(nodes);
        count += boundVerts[i];
    }

    for (int i = 0; i < nBounds; i ++){
        boundary(i).setMarker(boundMarker[i]);
        if (left[i] != -1) {
            boundary(i).setLeftCell(&cell(left[i]));
        }
        if (right[i] != -1) boundary(i).setRightCell(&cell(right[i]));
    }

    delete [] coords;
    delete [] nodeMarker;

    delete [] cellIdx;
    delete [] cellVerts;
    delete [] attribute;

    delete [] boundIdx;
    delete [] boundVerts;
    delete [] boundMarker;
    delete [] left;
    delete [] right;

    fclose(file);

}

template < class ValueType > void writeToFile(FILE * file, const ValueType & v, int count=1){
    if (!fwrite(&v, sizeof(ValueType), count, file)){
        __MS(v << " " << count)
        throwError(WHERE_AM_I + strerror(errno) + " " + str(errno));
    }
}

template < class ValueType > void readFromFile(FILE * file, ValueType & v, int count=1){
    uint ret = fread(&v, sizeof(ValueType), count, file);

    if (ret && ferror(file)){
        throwError(WHERE_AM_I + strerror(errno) + " " + str(errno));
    }
}

void Mesh::saveBinaryV2(const std::string & fbody) const {
//     std::cout << sizeof(uint8) << " "  << sizeof(uint16) << " " << sizeof(uint32) << " " << sizeof(uint64) << std::endl;

//   uint8[1] dimension
//   uint8[1] file format version
//   uint8[128]; // from v3 up
//   uint32[1] nVerts, number of vertices, max 2 ^ 32 (4e9)
//   double[3 * nVerts]; coordinates, dimension == 3 (x, y, z)
//   int32[nVerts] vertex markers [-2e9, .. , 2e9]
//   uint32[1] nCells, number of cell
//   uint8[nCells] cell nodeCount;
//   uint32[sum(cell nodeCount)] cellsidx

//   int32[nCells] cellMarkers [-2e9, .. , 2e9]
//   uint32[1] nBounds, number of boundarys
//   uint8[nBounds] bound nodeCount
//   uint32[sum(bound nodeCount)] boundsidx
//   int32[nBounds] boundaryMarkers [-2e9, .. , 2e9]
//   int32[nBounds] leftNeighbor idx (-1) if no neighbor present or info unavailable
//   int32[nBounds] rightNeighbor idx (-1) if no neighbor present or info unavailable

    std::string fileName(fbody.substr(0, fbody.rfind(MESHBINSUFFIX)) + MESHBINSUFFIX);

    FILE *file;
    file = fopen(fileName.c_str(), "w+b");
    if (!file) {
        throwError(WHERE_AM_I + " " + fileName + ": " + strerror(errno));
    }

    //** write preample
    writeToFile(file, uint8(this->dimension()));
    //** version

    int version = 3;
    writeToFile(file, uint8(version));

    uint8 * dummy = new uint8[128];
    std::memset(dummy, '\0', 128*sizeof(uint8));
    dummy[0] = (uint8)this->isGeometry();
    writeToFile(file, dummy[0], 128);

    //** write nodes
    double * coord = new double[3 * this->nodeCount()];
    for (uint i = 0; i < this->nodeCount(); i ++){
        for (uint j = 0; j < 3; j ++){
            coord[i * 3 + j] = node(i).pos()[j];
        }
    }

    int32 * marker = new int32[this->nodeCount()];
    for (uint i = 0; i < this->nodeCount(); i ++) marker[i] = node(i).marker();

    writeToFile(file, uint32(this->nodeCount()));

    if (this->nodeCount()){
        writeToFile(file, coord[0], 3 * this->nodeCount());
        writeToFile(file, marker[0], this->nodeCount());
    }

    //** write Cells
    uint32 nCells = this->cellCount();
    uint8 * cellVerts = new uint8[nCells];
    uint32 nCellIdx = 0;
    for (uint i = 0; i < nCells; i ++){
        cellVerts[i] = (uint8)cell(i).nodeCount();
        nCellIdx += cellVerts[i];
    }
    uint32 * cellIdx = new uint32[nCellIdx];
    uint count = 0;
    for (uint i = 0; i < nCells; i ++){
        for (uint j = 0; j < cellVerts[i]; j ++){
            cellIdx[count] = cell(i).node(j).id();
            count ++;
        }
    }
    int32 * cellMarker = new int32[nCells];
    for (uint i = 0; i < nCells; i ++) cellMarker[i] = cell(i).marker();

    writeToFile(file, nCells);

    if (nCells > 0){
        writeToFile(file, cellVerts[0], nCells);
        writeToFile(file, cellIdx[0], nCellIdx);
        writeToFile(file, cellMarker[0], nCells);
    }

    //** write boundarys
    uint32 nBound = this->boundaryCount();
    uint8 * boundVerts = new uint8[nBound];
    std::vector < uint32 > boundIdx;
    int32 * boundMarker = new int32[nBound];
    int32 * leftCells = new int32[nBound];
    int32 * rightCells = new int32[nBound];

    count = 0;
    for (uint i = 0; i < nBound;  i ++){
        Boundary * bound = &this->boundary(i);
        boundVerts[i] = (uint8)bound->nodeCount();

        for (uint j = 0; j < boundVerts[i]; j ++){
            boundIdx.push_back((uint32)bound->node(j).id());
        }

        boundMarker[i] = bound->marker();

        if (bound->leftCell()) {
            leftCells[i] = bound->leftCell()->id();
        } else leftCells[i] = -1;

        if (bound->rightCell()) {
            rightCells[i] = bound->rightCell()->id();
        } else rightCells[i] = -1;
    }

    writeToFile(file, nBound);
    if (nBound > 0){
        writeToFile(file, boundVerts[0], nBound);
        writeToFile(file, boundIdx[0], boundIdx.size());
        writeToFile(file, boundMarker[0], nBound);
        writeToFile(file, leftCells[0], nBound);
        writeToFile(file, rightCells[0], nBound);
    }

    writeToFile(file, dataMap_.size());
    if (dataMap_.size() > 0){
        for (auto & x: dataMap_){
            if (x.first.length() > 0 && x.second.size() > 0){
                writeToFile(file, x.first.length());
                writeToFile(file, x.first[0], x.first.length());
                writeToFile(file, x.second.size());
                writeToFile(file, x.second[0], x.second.size());
            } else {
                log(Warning, "Export data map invalid: " + x.first);
            }
        }
    }

    fclose(file);
    delete [] coord;
    delete [] marker;
    delete [] cellVerts;
    delete [] cellIdx;
    delete [] cellMarker;
    delete [] boundVerts;
    delete [] boundMarker;
    delete [] leftCells;
    delete [] rightCells;
}

void Mesh::loadBinaryV2(const std::string & fbody) {
    this->clear();
    std::string fileName(fbody.substr(0, fbody.rfind(MESHBINSUFFIX)) + MESHBINSUFFIX);

    FILE *file;
    file = fopen(fileName.c_str(), "r+b");

    if (!file) {
        throwError(WHERE_AM_I + " " + fileName + ": " + strerror(errno));
    }
    uint8 dim; readFromFile(file, dim);
    if (dim !=2 && dim !=3){
        throwError(WHERE_AM_I + " cannot determine dimension " + str(dim));
    }
    this->setDimension(dim);
    uint8 version; readFromFile(file, version);

    if (version == 3){
        uint8 *dummy = new uint8[128]; readFromFile(file, dummy[0], 128);
        this->setGeometry(bool(dummy[0]));
    } else if (version != 2){
        throwError(WHERE_AM_I + " wrong version " + str(version));
    }

    //** read nodes
    uint32 nVerts; readFromFile(file, nVerts);

    if (nVerts > 1e9){
        throwError(WHERE_AM_I + " probably something wrong: nVerts > 1e9 " + str(nVerts));
    }

    if (nVerts > 0){
        double * coord = new double[3 * nVerts]; readFromFile(file, coord[0], 3 * nVerts);
        int32 * marker = new int32[nVerts]; readFromFile(file, marker[0], nVerts);
        nodeVector_.reserve(nVerts);
        for (uint i = 0; i < nVerts; i ++) {
            this->createNode(coord[i * 3], coord[i * 3 + 1], coord[i * 3 + 2], marker[i]);
        }
        delete [] coord;
        delete [] marker;
    }

    //** read cells
    uint32 nCells; readFromFile(file, nCells);
    uint count = 0;

    if (nCells > 0){
        uint8 * cellVerts = new uint8[nCells]; readFromFile(file, cellVerts[0], nCells);
        uint nCellIdx = 0; for (uint i = 0; i < nCells; i ++) nCellIdx += cellVerts[i];
        uint32 * cellIdx = new uint32[nCellIdx];   readFromFile(file, cellIdx[0], nCellIdx);
        int32 * cellMarker = new int32[nCells]; readFromFile(file, cellMarker[0], nCells);

        //** create cells
        cellVector_.reserve(nCells);
        for (uint i = 0; i < nCells; i ++){
            std::vector < Node * > nodes(cellVerts[i]);
            for (uint j = 0; j < nodes.size(); j ++) nodes[j] = & node(cellIdx[count + j]);
            this->createCell(nodes, cellMarker[i]);
            count += cellVerts[i];
        }
        delete [] cellVerts;
        delete [] cellIdx;
        delete [] cellMarker;
    }

    //** read bounds
    uint32 nBound; readFromFile(file, nBound);
    if (nBound > 0){
        uint8 * boundVerts = new uint8[nBound];       readFromFile(file, boundVerts[0], nBound);
        uint nBoundIdx = 0; for (uint i = 0; i < nBound; i ++) nBoundIdx += boundVerts[i];
        uint32 * boundIdx = new uint32[nBoundIdx];    readFromFile(file, boundIdx[0], nBoundIdx);
        int32 * boundMarker = new int32[nBound];      readFromFile(file, boundMarker[0], nBound);
        int32 * leftCells = new int32[nBound];        readFromFile(file, leftCells[0], nBound);
        int32 * rightCells = new int32[nBound];       readFromFile(file, rightCells[0], nBound);

        //** create cells
        count = 0;
        boundaryVector_.reserve(nBound);
        for (uint i = 0; i < nBound; i ++){
            std::vector < Node * > nodes(boundVerts[i]);
            for (uint j = 0; j < nodes.size(); j ++) nodes[j] = & node(boundIdx[count + j]);

            Boundary * bound = this->createBoundary(nodes, boundMarker[i]);
            count += boundVerts[i];

            if (leftCells[i] > -1) bound->setLeftCell(&this->cell(leftCells[i]));
            if (rightCells[i] > -1) bound->setRightCell(&this->cell(rightCells[i]));
        }
        delete [] boundVerts;
        delete [] boundIdx;
        delete [] boundMarker;
        delete [] leftCells;
        delete [] rightCells;
    }

    size_t nData; readFromFile(file, nData);

    for (uint i = 0; i < nData; i ++){
        size_t strLen; readFromFile(file, strLen);
        std::string str; str.resize(strLen); readFromFile(file, str[0], strLen);
        size_t datLen; readFromFile(file, datLen);

        RVector dat(datLen); readFromFile(file, dat[0], datLen);
        this->addData(str, dat);
        //delete [] str;
    }

    fclose(file);
}

int Mesh::exportSimple(const std::string & fbody, const RVector & data) const {
  //output x y x y x y rhoa file
  std::fstream file; if (!openOutFile(fbody , & file)){ throwError("can't open file"); }

  for (uint i = 0; i < cellCount(); i ++){
    for (uint j = 0; j < 3; j ++){
      file << cell(i).node(j).x() << "\t" << cell(i).node(j).y() << "\t";
    }
    file << data[i] << std::endl;
  }
  file.close();
  return 1;
}

void Mesh::exportVTK(const std::string & fbody, bool writeCells) const {
    return exportVTK(fbody, dataMap_, PosVector(0), writeCells);
}

void Mesh::exportVTK(const std::string & fbody,
                     const PosVector & vec,
                     bool writeCells) const {
    return exportVTK(fbody, dataMap_, vec, writeCells);
}

void Mesh::exportVTK(const std::string & fbody,
                     const std::map< std::string,
                     RVector > & dataMap, bool writeCells) const {
    return exportVTK(fbody, dataMap, PosVector(0), writeCells);
}

void Mesh::exportVTK(const std::string & fbody, const RVector & arr) const {
    std::map< std::string, RVector > m(dataMap_);
    m.insert(std::pair< std::string, RVector >("arr", arr));
    return exportVTK(fbody, m, PosVector(0), true);
}

void Mesh::exportVTK(const std::string & fbody,
                     const std::map< std::string, RVector > & dataMap,
                     const PosVector & vec, bool cells) const {

    log(Debug, "Writing vtk: ", fbody + ".vtk");
    
    std::fstream file;
    if (!openOutFile(fbody.substr(0, fbody.rfind(".vtk")) + ".vtk", & file)) {
        return;
    }

    std::map< std::string, RVector > nData(dataMap);
    std::map< std::string, RVector > cData(dataMap);
    std::map< std::string, RVector > bData(dataMap);
    //std::map< std::string, R3Vector > vData(dataMap);

    for (auto & r: dataMap){
        if (r.second.size() == this->cellCount()){
            if (haveInfNaN(r.second)){
                log(Warning, "data: ", r.first, " contains inf or nan values .. skipping");    
            } else {
                cData[r.first] = r.second;
            }
        } else if (r.second.size() == this->boundaryCount()){
            if (haveInfNaN(r.second)){
                log(Warning, "data: ", r.first, " contains inf or nan values .. skipping");    
            } else {
                bData[r.first] = r.second;
            }
        } else if (r.second.size() == this->nodeCount()){
            if (haveInfNaN(r.second)){
                log(Warning, "data: ", r.first, " contains inf or nan values .. skipping");    
            } else {
                nData[r.first] = r.second;
            }
        } else {
            log(Warning, "data: ", r.first, " not written to vtk. size: ", r.second.size());
        }
    }

    if (cellCount() > 0){
        if (!cData.count("Marker")) cData["Marker"] = this->cellMarkers();
    }

    bool binary = false;
    file.precision(14);
    file << "# vtk DataFile Version 3.0" << std::endl;
    if (commentString_.size() > 0) {
        file << "d-" << dimension_ << "__ " << commentString_ << std::endl;
    } else {
        file << "d-" << dimension_ << "__ created by "
             << versionStr() << std::endl;
    }
    if (binary){
        file << "BINARY" << std::endl;
    } else {
        file << "ASCII" << std::endl;
    }
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    //** write nodes
    log(Debug, "writing points:", this->nodeCount());
    
    file << "POINTS " << nodeCount() << " double" << std::endl;

    for (Index i = 0; i < nodeCount(); i ++){
        if (binary){
            file.write((char*)&node(i).pos()[0], sizeof(double));
            file.write((char*)&node(i).pos()[1], sizeof(double));
            file.write((char*)&node(i).pos()[2], sizeof(double));
        } else {
            file << node(i).pos().x() << "\t"
                 << node(i).pos().y() << "\t"
                 << node(i).pos().z() << std::endl;
        }
    }
    if (binary){ file << std::endl; }

    //** write cells
    if (cells && cellCount() > 0){
        Index idxCount = 0;
        for (Index i = 0; i < cellCount(); i ++) {
            idxCount += cell(i).nodeCount() + 1;
        }
        log(Debug, "writing cell connections", cellCount());

        file << "CELLS " << cellCount() << " " << idxCount << std::endl;

        long iDummy;
        for (Index i = 0, imax = cellCount(); i < imax; i ++){
            if (binary){
                iDummy = cell(i).nodeCount();
                file.write((char*)&iDummy, sizeof(iDummy));
            } else {
                file << cell(i).nodeCount() << "\t";
            }
            if (cell(i).rtti() == MESH_TETRAHEDRON10_RTTI && oldTet10NumberingStyle_){
                file << cell(i).node(0).id() << "\t" << cell(i).node(1).id() << "\t"
                    << cell(i).node(2).id() << "\t" << cell(i).node(3).id() << "\t"
                    << cell(i).node(4).id() << "\t" << cell(i).node(7).id() << "\t"
                    << cell(i).node(5).id() << "\t" << cell(i).node(6).id() << "\t"
                    << cell(i).node(9).id() << "\t" << cell(i).node(8).id();
            } else {

                for (uint j = 0, jmax = cell(i).nodeCount(); j < jmax; j ++){
                    if (binary){
                        iDummy = cell(i).node(j).id();
                        file.write((char*)&iDummy, sizeof(iDummy));
                    } else {
                        file << cell(i).node(j).id() << "\t";
                    }
                }
            }
            if (!binary) file << std::endl;
        }
        if (binary) file << std::endl;
        log(Debug, "writing cell types");
        file << "CELL_TYPES " << cellCount() << std::endl;
        iDummy = 10;
        for (uint i = 0, imax = cellCount(); i < imax; i ++) {
            if (binary){
                file.write((char*)&iDummy, sizeof(iDummy));
            } else {
                switch (cell(i).rtti()){
                    case MESH_EDGE_CELL_RTTI:
                    case MESH_EDGE_RTTI: file          << "3 "; break;
                    case MESH_EDGE3_RTTI:
                    case MESH_EDGE3_CELL_RTTI: file    << "21 "; break;
                    case MESH_TRIANGLE_RTTI: file      << "5 "; break;
                    case MESH_TRIANGLE6_RTTI: file     << "22 "; break;
                    case MESH_QUADRANGLE_RTTI: file    << "9 "; break;
                    case MESH_QUADRANGLE8_RTTI: file   << "23 "; break;
                    case MESH_TETRAHEDRON_RTTI: file   << "10 "; break;
                    case MESH_TETRAHEDRON10_RTTI: file << "24 "; break;
                    case MESH_TRIPRISM_RTTI: file      << "13 "; break;
                    case MESH_TRIPRISM15_RTTI: file    << "13 "; break;
                    case MESH_PYRAMID_RTTI: file       << "14 "; break;
                    case MESH_PYRAMID13_RTTI: file     << "14 "; break;
                    case MESH_HEXAHEDRON_RTTI: file    << "12 "; break;
                    case MESH_HEXAHEDRON20_RTTI: file  << "25 "; break;
                    default:
                        std::cerr << WHERE_AM_I << " nothing known about."
                        << cell(i).rtti() << std::endl;
                }
            }
        }
        file << std::endl;

        //** write cell data
        if (cData.size() > 0){
            file << "CELL_DATA " << cellCount() << std::endl;
            for (auto & cd: cData){

                log(Debug, "writing cell data: " + cd.first + " " + str(cd.second.size()));

                if (cd.second.size() == (uint)cellCount()){
                    file << "SCALARS " << replace(cd.first, ' ', '_')
                            << " double 1" << std::endl;
                    file << "LOOKUP_TABLE default" << std::endl;

                    for (Index i = 0, imax = cd.second.size(); i < imax; i ++) {
                        if (binary){
                            //file.write((char*)&scaledValues[i], sizeof(double));
                        } else {
                            file << cd.second[i] << " ";
                        }
                    }
                    file << std::endl;
                }
            }
        }
    } else {  //   if !(cells && cellCount() > 0){
        //** write boundaries
        if (boundaryCount() > 0){
            Index idxCount = 0;
            for (Index i = 0; i < boundaryCount(); i ++) {
                idxCount += boundary(i).nodeCount() + 1;
            }

            log(Debug, "writing boundary connections");
            file << "CELLS " << boundaryCount() << " " << idxCount << std::endl;

            long iDummy;
            for (Index i = 0, imax = boundaryCount(); i < imax; i ++){
                if (binary){
            //        iDummy = cell(i).nodeCount();
            // 	      file.write((char*)&iDummy, sizeof(iDummy));
                } else {
                    file << boundary(i).nodeCount() << "\t";
                }

                for (Index j = 0, jmax = boundary(i).nodeCount();
                     j < jmax; j ++){

                    if (binary){
// 	          iDummy = cell(i).node(j).id();
// 	          file.write((char*)&iDummy, sizeof(iDummy));
                    } else {
                        file << boundary(i).node(j).id() << "\t";
                    }
                }
                if (!binary) file << std::endl;
            }
            if (binary) file << std::endl;

            log(Debug, "writing boundary types");
            file << "CELL_TYPES " << boundaryCount() << std::endl;
            iDummy = 10;
            for (Index i = 0, imax = boundaryCount(); i < imax; i ++) {
                if (binary){
                    file.write((char*)&iDummy, sizeof(iDummy));
                } else {
                    switch (boundary(i).rtti()){
                        case MESH_BOUNDARY_NODE_RTTI: file     << "1 "; break;
                        case MESH_EDGE_RTTI: file              << "3 "; break;
                        case MESH_EDGE3_RTTI: file             << "21 "; break;
                        case MESH_TRIANGLEFACE_RTTI: file      << "5 "; break;
                        case MESH_TRIANGLEFACE6_RTTI: file     << "22 "; break;
                        case MESH_QUADRANGLEFACE_RTTI: file    << "9 "; break;
                        case MESH_QUADRANGLEFACE8_RTTI: file   << "23 "; break;
                        case MESH_POLYGON_FACE_RTTI: file      << "7 "; break;
                    default:
                        std::cerr << WHERE_AM_I
                                  << " nothing known about." << boundary(i).rtti()
                                  << std::endl;
                    }
                }
            }
            file << std::endl;

            if (!bData.count("Marker")) bData["Marker"] = this->boundaryMarkers();
            
            if (bData.size() > 0){
            //** write boundary data
                file << "CELL_DATA " << boundaryCount() << std::endl;

                for (auto & bd: bData){

                    log(Debug, "writing boundary data: " + bd.first + " " + str(bd.second.size()));

                    if (bd.second.size() == (uint)boundaryCount()){

                        file << "SCALARS " << replace(bd.first, ' ', '_')
                            << " double 1" << std::endl;
                        file << "LOOKUP_TABLE default" << std::endl;

                        for (Index i = 0, imax = bd.second.size(); i < imax; i ++){
                            if (binary){
                            //file.write((char*)&scaledValues[i], sizeof(double));
                            } else {
                                file << bd.second[i] << " ";
                            }
                        }
                        file << std::endl;
                    }
                }
            }
            //** write boundary vector data
            if (vec.size() == boundaryCount()){
                log(Debug, "write vector field data");
                file << "VECTORS vec double" << std::endl;

                for (Index i = 0; i < vec.size(); i ++){
                    file << vec[i][0] << " "
                         << vec[i][1] << " "
                         << vec[i][2] << " " << std::endl;
                }
            }
        } // if (boundaryCount() > 0)

    } // else write boundaries

    //** write point data
    if (nData.size() > 0){
        log(Debug, "writing point data");
        file << "POINT_DATA " << nodeCount() << std::endl;

        for (auto & nd: nData){

            if (nd.second.size() == (uint)nodeCount()){
                log(Debug, "writing point data: " + nd.first + " " + str(nd.second.size()));

                file << "SCALARS " << replace(nd.first, ' ', '_')
                    << " double 1" << std::endl;
                file << "LOOKUP_TABLE default" << std::endl;

                for (Index i = 0, imax = nodeCount(); i < imax; i ++) {
                    file << nd.second[i] << " ";
                }
                file << std::endl;
            }
        }
    }

    //** write point vector data
    if (vec.size() == nodeCount()){
        log(Debug, "write vector field data: vec");
        file << "VECTORS vec double" << std::endl;

        for (Index i = 0; i < vec.size(); i ++){
            file << vec[i][0] << " "
                 << vec[i][1] << " "
                 << vec[i][2] << " " << std::endl;
        }
    } else {
        if (vec.size() > 0){
            std::cerr << "Vector data size does not match node size: "
                      << vec.size() << " " << nodeCount() << std::endl;
        }
    }

    file.close();
}

void Mesh::importVTK(const std::string & fbody) {
//     __MS(dimension_)
    this->clear();
    std::fstream file; openInFile(fbody.substr(0, fbody.rfind(".vtk")) + ".vtk",
                                  &file);

    std::vector < std::string > row;
    getline(file, commentString_); //** vtk version line
    getline(file, commentString_); //** comment line


    if (commentString_.find("d-2__") != std::string::npos){
        dimension_ = 2;
    } else if (commentString_.find("d-3__") != std::string::npos){
        dimension_ = 3;
    }

    while (!file.eof()){
        row = getRowSubstrings(file);
      //  std::cout << row.size() << std::endl;
        if (row.size()){
            if (row[0] == "ASCII"){
                break;
            } else if (row[0] == "BINARY"){
                THROW_TO_IMPL
                break;
            }
        }
    }
    row = getRowSubstrings(file);
    //std::cout << row.size() << std::endl;
    if (row.back() == "UNSTRUCTURED_GRID"){

        //** End reading header
        while (!file.eof()){
            row = getRowSubstrings(file);
            if (row.size()){
                if (row[0] == "POINTS") readVTKPoints_(file, row);
                else if (row[0] == "CELLS") readVTKCells_(file, row);
                else if (row[0] == "SCALARS") readVTKScalars_(file, row);
                else if (row[0] == "FIELD") {
                    if (row.size() == 3){
                        std::string dataName(row[1]);
                        Index numArray(toInt(row[2]));
                        for (Index i = 0; i < numArray; i ++ ){
                            std::vector < std::string > r(getNonEmptyRow(file));

                            if (r.size() == 4){
                                std::string arrayName(r[0]);
                                Index numComponents(toInt(r[1]));
                                Index numTuples(toInt(r[2]));
                                std::string dataType(r[3]);

                                for (Index j = 0; j < numComponents; j ++ ){
                                    //u 1 2601 double
                                    RVector field(numTuples);
                                    for (Index k = 0; k < numTuples; k ++) {
                                        file >> field[k];
                                    }

                                    addData(dataName + "_" + str(i) + "_" +
                                            arrayName + "_" + str(j), field);
                                }
                            }
                        }
                    }
                }
                else if (row[0] == "VECTORS") {
                    log(Warning, " .. clean up VTK vector data read");
                    std::string name(row[1]);
                    std::vector < std::string > r(getRowSubstrings(file));
                    if (r.size()){
                        if (r[0] == "LOOKUP_TABLE"){
                            r = getRowSubstrings(file);
                        }
                    }
                    RVector data_x(nodeCount());
                    RVector data_y(nodeCount());
                    RVector data_z(nodeCount());
                    for (uint i = 0; i < data_x.size(); i ++) {
                        // __MS(r.size())
                        if (r.size() == 3){
                            // __MS(r[0] << " " << r[1] << " " << r[2])
                            data_x[i] = toDouble(r[0]);
                            data_y[i] = toDouble(r[1]);
                            data_z[i] = toDouble(r[2]);
                        }
                        r = getRowSubstrings(file);
                    }
                    addData(name + "_x", data_x);
                    addData(name + "_y", data_y);
                    addData(name + "_z", data_z);
                }
            }
        }

    } else if (row.back() == "POLYDATA"){
        while (!file.eof()){
            row = getRowSubstrings(file);
            if (row.size()){
                if (row[0] == "POINTS") readVTKPoints_(file, row);
                else if (row[0] == "POLYGONS") readVTKPolygons_(file, row);
            }
        }
    } else if (row.back() == "STRUCTURED_GRID"){
        Index nx = 0;
        Index ny = 0;
        Index nz = 0;

        while (!file.eof()){
            row = getRowSubstrings(file);

            if (row.size()){
                if (row[0] == "DIMENSIONS"){
                    if (row.size() == 4){
                        nx = toInt(row[1]);
                        ny = toInt(row[2]);
                        nz = toInt(row[3]);
                    } else {
                        __MS(row)
                        THROW_TO_IMPL
                    }
                } else if (row[0] == "POINTS"){
                    //POINTS 1331 double
                    Index nVerts = toInt(row[1]);
                    RVector vx(nVerts);
                    RVector vy(nVerts);
                    RVector vz(nVerts);
                    for (Index i = 0; i < nVerts; i ++) {
                        file >> vx[i] >> vy[i] >> vz[i];
                    }

                    RVector gvx(nx+1), gvy(ny+1), gvz(nz+1);
                    for (Index i = 0; i < nx+1; i ++){
                        gvx[i] = min(vx) + i * (max(vx) - min(vx))/nx;
                    }
                    for (Index i = 0; i < ny+1; i ++){
                        gvy[i] = min(vy) + i * (max(vy) - min(vy))/ny;
                    }
                    for (Index i = 0; i < nz+1; i ++){
                        gvz[i] = min(vz) + i * (max(vz) - min(vz))/nz;
                    }

                    this->create3DGrid(gvx, gvy, gvz);

                } else if (row[0] == "SCALARS") {
                    readVTKScalars_(file, row);
                }
            }
        }
    } else {
        __MS(row)
        THROW_TO_IMPL
    }
//     __MS(dimension_)
    //this->showInfos();
    file.close();

}

void Mesh::readVTKPoints_(std::fstream & file,
                          const std::vector < std::string > & row){
    uint nVerts = toInt(row[1]);
    //std::cout << "nVerts: " << nVerts << std::endl;
    double x = 0.0, y = 0.0, z = 0.0;
    for (uint i = 0; i < nVerts; i ++) {
        file >> x >> y >> z;
        this->createNode(x, y, z);
    }

//     std::cout << GIMLI::z(positions()) << std::endl;
//     std::cout << (min(abs(GIMLI::z(positions()))) << std::endl;
//     __MS(nonZero(GIMLI::x(positions())))
//     __MS(nonZero(GIMLI::y(positions())))
//     __MS(nonZero(GIMLI::z(positions())))
    if (!nonZero(GIMLI::y(positions())) && nonZero(GIMLI::z(positions()))){
        dimension_ = 2;
        // swap y and z
        for (Index i = 0; i < nodeCount(); i ++ ){
            nodeVector_[i]->pos()[1] = nodeVector_[i]->pos()[2];
            nodeVector_[i]->pos()[2] = 0.0;
        }

    } else if (!nonZero(GIMLI::z(positions()))) dimension_ = 2;
}

void Mesh::readVTKCells_(std::fstream & file,
                         const std::vector < std::string > & row){
    uint nCells = toInt(row[1]);
    //std::cout << "nCells: " << nCells << std::endl;
    uint nNodes = 0;
    uint id;
    std::vector < Node * > nodes;
    for (uint i = 0; i < nCells; i ++) {
        file >> nNodes;
        nodes.resize(nNodes);
        for (uint j = 0; j < nNodes; j ++) {
            file >> id;
            nodes[j] = &this->node(id);
        }
        this->createCell(nodes);
    }
}

void Mesh::readVTKPolygons_(std::fstream & file, const std::vector < std::string > & row){
    uint nPoly = toInt(row[1]);
    uint nNodes = 0;
    uint id;
    std::vector < Node * > nodes;
    for (uint i = 0; i < nPoly; i ++) {
        file >> nNodes;
        nodes.resize(nNodes);
        for (uint j = 0; j < nNodes; j ++) {
            file >> id;
            nodes[j] = &this->node(id);
        }
        this->createBoundary(nodes);
    }
}

// template <class Arg, class Result> class bind : std::unary_function< Arg, Result >{
// public:
//     Result operator() (Arg & o) { return toDouble(o); }
// };
//   struct unary_function{
//     typedef Arg argument_type;
//     typedef Result result_type;
//   };
// template < class ret > class bind{
// public:
//     bind() {
//         funct_ = &toDouble;
//     }
//     template < class type > ret operator() (type & o) { return toDouble(o); }
// protected:
//     //void (Region::*)
//     size_t &funct_;
//     //(toDouble*) funct_;
// };

void Mesh::readVTKScalars_(std::fstream & file, const std::vector < std::string > & row){
    std::string name(row[1]);
    // __MS(str(row[0]) + " " + str(row[1]))
    std::vector < std::string > r(getRowSubstrings(file));
    // __MS(r)
    if (r.size()){
        if (r[0] == "LOOKUP_TABLE"){
            r = getRowSubstrings(file);
        }
    }
    // __MS(r)
    // __MS(r.size())
    // __MS((cellCount() > 1 || boundaryCount() > 1 || nodeCount() > 1))

    std::vector < std::string > ri;
    
    if (r.size() == 1 && (cellCount() > 1 || boundaryCount() > 1 || nodeCount() > 1)){
        // __M
        bool go=true;
        while (go){
            ri = getNonEmptyRow(file);
            if (ri.size() == 0){
                go=false;
            } else if (ri.size() > 1){
                file.unget();
                go=false;
            } else{
                r.push_back(ri[0]);
            }
        }
        // std::vector < std::string > r(getRowSubstrings(file));
        // THROW_TO_IMPL
    }

    // __MS(r)
    // __MS(r.size())

    RVector data(r.size());
//     std::cout.precision(14);
//     std::cout << *r.begin() << std::endl;
//     std::cout << toDouble(*r.begin()) << std::endl;
//     std::cout << std::bind1st(&toDouble)(*r.begin()) << std::endl;
    //std::cout << bind< double, toDouble >(&toDouble)(*r.begin()) << std::endl;

    //std::copy(r.begin(), r.end(), data.begin(), bind< double >(toDouble);
    for (uint i = 0; i < data.size(); i ++) data[i] = toDouble(r[i]);
    
    addData(name, data);
    if (name == "Marker"){
        if (data.size() == this->cellCount()) this->setCellMarkers(data);
        if (data.size() == this->boundaryCount()) this->setBoundaryMarkers(data);
    }
}


void Mesh::importVTU(const std::string & fbody) {
    this->clear();
    THROW_TO_IMPL
}

void Mesh::exportVTU(const std::string & fbody, bool binary) const {
    std::string filename(fbody);
    if (filename.rfind(".vtu") == std::string::npos){
        filename = fbody.substr(0, filename.rfind(".vtk")) + ".vtu";
    }
    std::fstream file; if (!openOutFile(filename, & file)) { return ; }
    file.precision(14);
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    file << "<UnstructuredGrid>" << std::endl;

    std::map< std::string, RVector > data(dataMap_);
    if (cellCount() > 0){
        if (!data.count("_Marker")) {
            data.insert(std::make_pair("_Marker",  this->cellMarkers()));
        }
        if (!data.count("_Attribute")) data.insert(std::make_pair("_Attribute",  cellAttributes()));
    }
    addVTUPiece_(file, *this, data);

    file << "</UnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;
    file.close();
}

void Mesh::exportBoundaryVTU(const std::string & fbody, bool binary) const {
    std::string filename(fbody);
    if (filename.rfind(".vtu") == std::string::npos){
        filename = fbody.substr(0, filename.rfind(".vtk")) + ".vtu";
    }
    std::fstream file; if (!openOutFile(filename, & file)) { return ; }

    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    file << "<UnstructuredGrid>" << std::endl;

    std::vector < Boundary * > bs;
    for (uint i = 0; i < boundaryCount(); i ++) {
        if (boundary(i).marker() != 0.0) {
            bs.push_back(&boundary(i));
        }
//         if (boundary(i).marker() == -2) {
//             std::cout << boundary(i)<< " " << boundary(i).center() << " " << boundary(i).marker() << std::endl;
//         }
    }

    Mesh boundMesh;
    boundMesh.createMeshByBoundaries(*this, bs);
//    boundMesh.showInfos();
    //for (uint i =0; i < boundMesh.nodeCount(); i ++) std::cout << boundMesh.node(i)<< std::endl;
    std::map< std::string, RVector > boundData;

    if (!boundData.count("_BoundaryMarker")) {
        boundData.insert(std::make_pair("_BoundaryMarker", 
                         boundMesh.boundaryMarkers()));
    }

    //boundMesh.exportVTK(fbody, boundData);
    addVTUPiece_(file, boundMesh, boundData);

    file << "</UnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;
    file.close();
}

void Mesh::addVTUPiece_(std::fstream & file, const Mesh & mesh,
                  const std::map < std::string, RVector > & data) const{

    bool binary = false;
    bool cellsAreBoundaries = false;

    std::vector < MeshEntity * > cells;

    if (mesh.cellCount() == 0 && mesh.boundaryCount() > 0){
        cellsAreBoundaries = true;
        cells.reserve(mesh.boundaryCount());
        for (uint i = 0; i < mesh.boundaryCount(); i ++) cells.push_back(& mesh.boundary(i));
    } else {
        cells.reserve(mesh.cellCount());
        for (uint i = 0; i < mesh.cellCount(); i ++) cells.push_back(& mesh.cell(i));
    }

    uint nNodes = mesh.nodeCount();
    uint nCells = cells.size();

    file << "<Piece NumberOfPoints=\"" << nNodes << "\" NumberOfCells=\"" << nCells << "\">" << std::endl;

    file << "<Points>" << std::endl;
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for (uint i = 0; i < mesh.nodeCount(); i ++) {
        file << mesh.node(i).x() << " " << mesh.node(i).y() << " " << mesh.node(i).z() << " ";
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "</Points>" << std::endl;

    file << "<Cells>" << std::endl;
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    for (uint i = 0; i < nCells; i ++) {
        MeshEntity * cell = cells[i];
        if (cell->rtti() == MESH_TETRAHEDRON10_RTTI){

            file << cell->node(0).id() << " " << cell->node(1).id() << " "
                << cell->node(2).id() << " " << cell->node(3).id() << " "
                << cell->node(4).id() << " " << cell->node(7).id() << " "
                << cell->node(5).id() << " " << cell->node(6).id() << " "
                << cell->node(9).id() << " " << cell->node(8).id() << " ";
        } else for (uint j = 0; j < cell->nodeCount(); j ++) file << cell->node(j).id() << " ";
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    uint count = 0;

    for (uint i = 0; i < nCells; i ++) {
        count += cells[i]->nodeCount();
        file << count << " ";
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;

    for (uint i = 0; i < nCells; i ++) {
        switch (cells[i]->rtti()){
            case MESH_BOUNDARY_NODE_RTTI: file     << "1 "; break;
            case MESH_EDGE_CELL_RTTI:
            case MESH_EDGE_RTTI: file          << "3 "; break;
            case MESH_EDGE3_CELL_RTTI:
            case MESH_EDGE3_RTTI: file         << "21 "; break;
            case MESH_TRIANGLEFACE_RTTI:
            case MESH_TRIANGLE_RTTI: file      << "5 "; break;
            case MESH_TRIANGLEFACE6_RTTI:
            case MESH_TRIANGLE6_RTTI: file     << "22 "; break;
            case MESH_QUADRANGLEFACE_RTTI:
            case MESH_QUADRANGLE_RTTI: file    << "9 "; break;
            case MESH_QUADRANGLEFACE8_RTTI:
            case MESH_QUADRANGLE8_RTTI: file   << "23 "; break;
            case MESH_TETRAHEDRON_RTTI: file   << "10 "; break;
            case MESH_TETRAHEDRON10_RTTI: file << "24 "; break;
            case MESH_HEXAHEDRON_RTTI: file    << "12 "; break;
            case MESH_HEXAHEDRON20_RTTI: file  << "25 "; break;
            case MESH_POLYGON_FACE_RTTI: file  << "7 "; break; // VTK_POLYGON
            default: std::cerr << WHERE_AM_I << " nothing know about." << cells[i]->rtti() << std::endl;
        }
    }
    file << std::endl << "</DataArray>" << std::endl;
    file << "</Cells>" << std::endl;

    uint nodeData = 0, cellData = 0;

    for (std::map < std::string, RVector >::const_iterator it = data.begin();
         it != data.end(); it ++){


        if (it->second.size() == nNodes && !cellsAreBoundaries) {
            //NodeCount == Cellcount for cellsAreBoundaries(2d)
            nodeData++;
        } else if (it->second.size() == nCells) {
            cellData++;
        } else {
            std::cerr << WHERE_AM_I << " dont know how to handle data array: " << it->first
                        << " with size " << it->second.size() << " nodesize = " << nNodes
                        << " cellsize = " << nCells << std::endl;
        }
    }

    if (nodeData > 0){
        file << "<PointData>" << std::endl;
        for (std::map < std::string, RVector >::const_iterator it = data.begin();
              it != data.end(); it ++){
            if (it->second.size() == nNodes && !cellsAreBoundaries ) {
                 file << "<DataArray type=\"Float32\" Name=\""
                      << it->first << "\" format=\"ascii\">" << std::endl;
                 for (uint i = 0; i < it->second.size(); i ++) file << it->second[i] << " ";
                 file << std::endl << "</DataArray>" << std::endl;
            }
        }
        file << "</PointData>" << std::endl;
    }

     if (cellData > 0){
        file << "<CellData>" << std::endl;
        for (std::map < std::string, RVector >::const_iterator it = data.begin();
             it != data.end(); it ++){

            if (it->second.size() == nCells) {
                 file << "<DataArray type=\"Float64\" Name=\"" << it->first;
                 if (!binary){
                    file << "\" format=\"ascii\">" << std::endl;
                    for (uint i = 0; i < it->second.size(); i ++) file << it->second[i] << " ";
                 } else {
                    file << "\" format=\"binary\">";
                    file.write((char*)&it->second[0], it->second.size() * sizeof(double));
                 }


                 file << std::endl << "</DataArray>" << std::endl;
            }
        }
        file << "</CellData>" << std::endl;
    }

    file << "</Piece>" << std::endl;
}

void Mesh::importMod(const std::string & filename){
    RMatrix mat;
    std::vector < std::string > comments;
    loadMatrixCol(mat, filename, comments);
//     std::cout << comments.size() << std::endl;
//     std::cout << mat.size() << std::endl;
//     std::cout << mat[0].size() << std::endl;

    RVector x(unique(sort(cat(mat[0], mat[1]))));
    RVector y(unique(sort(cat(mat[2], mat[3]))) * -1.0);
    create2DGrid(x, y);
    if (comments.size() > 4 && mat.cols() > 4) addData(comments[4], mat[4]);
    if (comments.size() > 5 && mat.cols() > 5) addData(comments[5], mat[5]);
}

void Mesh::importSTL(const std::string & fileName, bool isBinary, double snap){
    double tolerance = snap;

    std::vector< std::vector < RVector3 > > allVerts;

    if (!isBinary){ // try import ascii
        std::fstream file; openInFile(fileName, & file);

        std::vector < std::string > row;

        row = getNonEmptyRow(file);
        if (row[0] != "solid") {
            file.close();
            importSTL(fileName, true);
        }
        allVerts.push_back(std::vector < RVector3 >());

        while (!file.eof()){
            row = getNonEmptyRow(file);
            if (row.size() > 0){
                if (row[0] == "solid"){
                    allVerts.push_back(std::vector < RVector3 >());

                } else if (row[0] == "facet" && row[1] == "normal"){ //** facet normal  0.0  0.0  0.0
                    row = getNonEmptyRow(file); //** outer loop
                    row = getNonEmptyRow(file); //** vertex x y z;
                    allVerts.back().push_back(RVector3(toDouble(row[1]), toDouble(row[2]), toDouble(row[3])));
                    row = getNonEmptyRow(file); //** vertex x y z
                    allVerts.back().push_back(RVector3(toDouble(row[1]), toDouble(row[2]), toDouble(row[3])));
                    row = getNonEmptyRow(file); //** vertex x y z
                    allVerts.back().push_back(RVector3(toDouble(row[1]), toDouble(row[2]), toDouble(row[3])));
                    row = getNonEmptyRow(file);  //** endloop
                    row = getNonEmptyRow(file);  //** endfacet;
                } else if (row[0] == "endsolid"){
                }
            }
        }
        file.close();
    } else { // import Binary Format
        FILE * file; file = fopen(fileName.c_str(), "r+b");
        char header[80];
        Index ret = 0;
        ret = fread(&header, 1, 80, file);
        if (ret == 0) throwError(WHERE_AM_I + " Oops");
        int nFaces = 0;
        ret = fread(&nFaces, 4, 1, file);
        if (ret == 0) throwError(WHERE_AM_I + " Oops");

        float rd[48];
        char padding[2];
        for (int i = 0; i < nFaces; i ++){
            std::vector < RVector3 > face;
            ret = fread(&rd, 4, 12, file);
            if (ret == 0) throwError(WHERE_AM_I + " Oops");

            face.push_back(RVector3(rd[3], rd[4],  rd[5]));
            face.push_back(RVector3(rd[6], rd[7],  rd[8]));
            face.push_back(RVector3(rd[9], rd[10], rd[11]));
            allVerts.push_back(face);

            ret = fread(&padding, 1, 2, file);
            if (ret == 0) throwError(WHERE_AM_I + " Oops");
        }
        fclose(file);
    } // end import binary STL format

    Node *n1, *n2, *n3;
    // __MS(allVerts.size())
    for (Index j = 0; j < allVerts.size(); j ++ ){
        if (allVerts[j].size() % 3 == 0 && allVerts[j].size() > 0){
            for (uint i = 0; i < allVerts[j].size() / 3; i ++){
                n1 = createNodeWithCheck(allVerts[j][i * 3], tolerance);
                n2 = createNodeWithCheck(allVerts[j][i * 3 + 1], tolerance);
                n3 = createNodeWithCheck(allVerts[j][i * 3 + 2], tolerance);
                this->createTriangleFace(*n1, *n2, *n3, j);
            }
        } else {
            throwError(WHERE_AM_I + " there is something wrong in ascii-stl-format "
                    + str(allVerts.size()) + " " + str(allVerts.size() % 3));
        }
    }
}

int Mesh::exportMidCellValue(const std::string & fileName,
                              const RVector & data1, const RVector & data2) const {

    RMatrix mat(dimension_, cellCount());
    for (uint j = 0; j < cellCount(); j ++) {
        for (uint i = 0; i < dimension_; i ++) {
            mat[i][j] = cell(j).center()[i];
        }
    }

    if (data1.size() == cellCount()) mat.push_back(data1);
    if (data2.size() == cellCount()) mat.push_back(data2);

    return saveMatrixCol(mat, fileName);
}

void Mesh::exportAsTetgenPolyFile(const std::string & filename){
    std::fstream file; openOutFile(filename.substr(0, filename.rfind(".poly")) + ".poly", & file);

    uint nverts = nodeCount();
    uint nfacets = boundaryCount();

    //  file << "# nverties dimension nattrib boolbndrymarker" << endl;
    file << nverts << "\t3\t0\t1" << std::endl;

    //file << "#    pointNr. x y z marker" << endl;
    file.setf(std::ios::scientific, std::ios::floatfield);
    file.precision(12);
    for (uint i = 0; i < nverts; i++){
        file << i << "\t" << node(i).x()
                << "\t" << node(i).y()
                << "\t" << node(i).z()
                << "\t" << node(i).marker() << std::endl;
    }

    //file << "# nfacets boolbndrymarker" << endl;
    size_t nPolyFacet = 0;
    size_t nNodesPoly = 0;
    file << nfacets << "\t1" << std::endl;
    for (size_t i = 0, imax = nfacets;i < imax; i ++){

        nPolyFacet = 1;
        file << nPolyFacet << "\t0\t" << boundary(i).marker() << std::endl;
        //file << nPolyFacet << "\t" << facet(i).holeCount() << "\t" << facet(i).boundaryMarker() << endl;

        nNodesPoly = boundary(i).nodeCount();
        file << nNodesPoly << " \t";

        for (size_t k = 0, kmax = nNodesPoly; k < kmax; k++){
            file << boundary(i).node(k).id() << "\t";
        }
        file << std::endl;
    }

    //file << "# nholes" << endl;
    file << 0 << std::endl;
    //file << "#    hole x y z" << endl;
//     for (int i = 0, imax = vecpHoles_.size(); i < imax ; i ++){
//         file << i << "\t"
//         << vecpHoles_[i]->x() << "\t"
//         << vecpHoles_[i]->y() << "\t"
//         << vecpHoles_[i]->z() << endl;
//     }

  //file << "# nregions" << endl;
  file << 0 << std::endl;
  //file << "#    region x y z attribute maxarea" << endl;
//   for (int i = 0, imax = regions_.size(); i < imax ; i ++){
//     file << i << "\t"
//      << regions_[i]->x() << "\t"
//      << regions_[i]->y() << "\t"
//          << regions_[i]->z() << "\t"
//      << regions_[i]->attribute() << "\t"
//      << regions_[i]->dx() << endl;
//   }
  file.close();

}

} // namespace GIMLI{
