/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@gimli.org                                         *
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

#ifndef _GIMLI_MESH__H
#define _GIMLI_MESH__H

#include "gimli.h"
#include "meshentities.h"
#include "node.h"
#include "pos.h"
#include "utils.h"
#include "shape.h"

#include <vector>
#include <list>
#include <set>
#include <map>
#include <fstream>

namespace GIMLI{

class KDTreeWrapper;

//! A BoundingBox
/*! A BoundingBox which contains a min and max Vector3< double >*/
class DLLEXPORT BoundingBox{
public:
    /*! Default constructor, with BoundingBox[0.0, 0.0, 0.0; 1.0, 1.0, 1.0] */
    BoundingBox(const Pos & min=Pos(0, 0, 0),
                const Pos & max=Pos(1.0, 1.0, 1.0))
        : min_(min), max_(max){
    }

    /*! Construct BBox from position vector */
    BoundingBox(const PosVector & vPos){
        min_ = Pos((double)MAX_DOUBLE, (double)MAX_DOUBLE, (double)MAX_DOUBLE);
        max_ = min_ * -1.0;
        for (uint i = 0; i < vPos.size(); i ++){
            min_[0] = std::min(vPos[i][0], min_[0]);
            min_[1] = std::min(vPos[i][1], min_[1]);
            min_[2] = std::min(vPos[i][2], min_[2]);
            max_[0] = std::max(vPos[i][0], max_[0]);
            max_[1] = std::max(vPos[i][1], max_[1]);
            max_[2] = std::max(vPos[i][2], max_[2]);
        }
    }

    /*! Copy constructor. */
    BoundingBox(const BoundingBox & bbox){ this->copy_(bbox); }

    /*! Assignment operator. */
    BoundingBox & operator = (const BoundingBox & bbox){
    	if (this != & bbox){
            this->copy_(bbox);
        } return * this;
    }

    /*! Check if a point lie inside (with boundary). */
    bool isInside(const Pos & p){
        return ((p[0] <= max_[0] && p[0] >= min_[0]) &&
                (p[1] <= max_[1] && p[1] >= min_[1]) &&
                (p[2] <= max_[2] && p[2] >= min_[2]));
    }

    /*! Set minimum Position. */
    void setMin(const Pos & min) { min_ = min; }
    /*! Return a copy of the minimum position. */
    const Pos & min() const { return min_; }

    /*! Set maximum Position. */
    void setMax(const Pos & max) { max_ = max; }
    /*! Return a copy of the maximum position. */
    const Pos & max() const { return max_; }

    /*! Return minimal x coordinate.*/
    inline double xMin() const { return min_[0];}
    /*! Return minimal y coordinate.*/
    inline double yMin() const { return min_[1];}
    /*! Return minimal z coordinate.*/
    inline double zMin() const { return min_[2];}

    /*! Return maximal x coordinate.*/
    inline double xMax() const { return max_[0];}
    /*! Return maximal y coordinate.*/
    inline double yMax() const { return max_[1];}
    /*! Return maximal z coordinate.*/
    inline double zMax() const { return max_[2];}

    /*! Return maximal x length.*/
    inline double xSize() const { return abs(this->xMax()-this->xMin());}
    /*! Return maximal y length.*/
    inline double ySize() const { return abs(this->yMax()-this->yMin());}
    /*! Return maximal z length.*/
    inline double zSize() const { return abs(this->zMax()-this->zMin());}

protected:

    void copy_(const BoundingBox & bbox){
    	min_ = bbox.min();
        max_ = bbox.max();
    }

    Pos min_;
    Pos max_;
};

inline std::ostream & operator << (std::ostream & str, const BoundingBox & bb){
    str << "BoundingBox [" << bb.min() << ", " << bb.max() << "]" ;
    return str;
}

DLLEXPORT std::ostream & operator << (std::ostream & str, const Mesh & mesh);

class DLLEXPORT Mesh {

public:
    typedef std::vector< RegionMarker > RegionMarkerList;
    typedef RVector3 HoleMarker;
    typedef PosVector HoleMarkerList;

    /*! Default constructor, create empty mesh with dimension dim
    If this mesh is supposed to be a geometry definition, all
    created nodes will be checked for duplicates.*/
    Mesh(Index dim=2, bool isGeometry=false);

    /*! Constructor, read mesh from filename */
    Mesh(const std::string & filename, bool createNeighborInfos=true);

    /*! Copy constructor. */
    Mesh(const Mesh & mesh);

    /*! Copy assignment operator. */
    Mesh & operator = (const Mesh & mesh);

    /*! Default destructor. */
    ~Mesh();

    /*! Clear all data, inclusive all caches.*/
    void clear();

    /*!If the mesh is static in geometry and shape some useful information are cached.
     * (cell sizes, boundary sizes, ...)
     For dynamic meshes, i.e., node positions can be moved, you have to set staticGeometry to false to avoid any caching.*/
    void setStaticGeometry(bool stat);

    /*! Return true if this mesh have static geometry. [Default=True]*/
    inline bool staticGeometry() const { return staticGeometry_; }

    /*! Mesh is marked as geometry definition or PLC
    so createNode will allways with check. */
    void setGeometry(bool b);

    /*! Return if the mesh is a geometry definition.*/
    bool isGeometry() const { return isGeometry_; }

    /*! Some parts of the geometry changed so the mesh is supposed to be dynamic.*/
    void geometryChanged();

    /*! Set the dimension of the mesh. [Default = 2] */
    void setDimension(uint dim){ dimension_ = dim;}

    /*! Return the dimension of this mesh.*/
    uint dimension() const { return dimension_; }

    /*! Shortcut for \ref dimension.*/
    uint dim() const { return dimension_; }

    //** start creation stuff
    Node * createNode(double x, double y, double z, int marker=0);
    Node * createNode(const Node & node);
    Node * createNode(const RVector3 & pos, int marker=0);

    /*! Create a secondary node, which is stored in an additional list for additional use.
    If tolerance tol set to a value > 0, then it will be checked if there is already a node
    at this position and return a ptr to the existing node instead of creating a new. */
    Node * createSecondaryNode(const RVector3 & pos, double tol=-1);

    /*! Create new Node with duplication checks. Returns the already existing node when its within a tolerance distance to pos.
    If edgeCheck is set, any 2d (p1) boundary edges will be checked for any intersection with pos and splitted if necessary.*/
    Node * createNodeWithCheck(const RVector3 & pos, double tol=1e-6,
                               bool warn=false, bool edgeCheck=false);

    Boundary * createBoundary(std::vector < Node * > & nodes, int marker=0, bool check=true);
    /*! Create a boundary from the given node indices */
    Boundary * createBoundary(const IndexArray & nodes, int marker=0, bool check=true);
    Boundary * createBoundary(const Boundary & bound, bool check=true);
    Boundary * createBoundary(const Cell & cell, bool check=true);
    Boundary * createNodeBoundary(Node & n1, int marker=0, bool check=true);
    Boundary * createEdge(Node & n1, Node & n2, int marker=0, bool check=true);
    Boundary * createEdge3(Node & n1, Node & n2, Node & n3, int marker=0, bool check=true);
    Boundary * createTriangleFace(Node & n1, Node & n2, Node & n3, int marker=0, bool check=true);
    Boundary * createQuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4, int marker=0, bool check=true);
    Boundary * createPolygonFace(std::vector < Node * > & nodes, int marker, bool check=true);

    /*! Create empty cell without a node or a shape. */
    Cell * createCell(int marker=0);
    Cell * createCell(std::vector < Node * > & nodes, int marker=0);
    /*! Create a cell from the given node indices */
    Cell * createCell(const IndexArray & nodes, int marker=0);
    Cell * createCell(const Cell & cell);
    Cell * createTriangle(Node & n1, Node & n2, Node & n3, int marker=0);
    Cell * createQuadrangle(Node & n1, Node & n2, Node & n3, Node & n4, int marker=0);
    Cell * createTetrahedron(Node & n1, Node & n2, Node & n3, Node & n4, int marker=0);

    /*! Create a cell as a copy of a cell from an alternative mesh.
     * Each node of the new cell will be created with duplication check and
     * reused if there is already a node withing the tolerance distance tol.
     * tol=-1 disables this duplication check. */
    Cell * copyCell(const Cell & cell, double tol=1e-6);
    /*! Create a boundary as a copy of a boundary from an alternative mesh.
     * Each node of the new cell will be created with duplication check and
     * reused if there is already a node withing the tolerance distance tol.
     * tol=-1 disables this duplication check. */
    Boundary * copyBoundary(const Boundary & bound, double tol=1e-6, bool check=true);

    void create1DGrid(const RVector & x);

    /*! Default boundary marker are -x[1], +x[2], +z[3], -z[4].
     If worldBoundaryMarker is set it becomes +z[-1], else[-2]. */
    void create2DGrid(const RVector & x, const RVector & y, int markerType=0,
                      bool worldBoundaryMarker=false);

    /*! Default boundary marker are -x[1], +x[2], +z[3], -z[4], -y[5], +z[6].
     If worldBoundaryMarker is set it becomes +z[-1], else[-2]. You can \ref exportBoundaryVTU
     and take a look with Paraview. */
    void create3DGrid(const RVector & x, const RVector & y, const RVector & z,
                      int markerType=0, bool worldBoundaryMarker=false);

    /*! Create one dimensional grid.
     * Boundary on the domain border will get
     * marker: 1,2 for: left, right.*/
    void createGrid(const RVector & x) { create1DGrid(x); }

    /*! Create two dimensional grid. Boundaries on the domain border will
    get the markers: 1,2,3,4 for: left, right, bottom, top.*/
    void createGrid(const RVector & x, const RVector & y,
                    int markerType=0, bool worldBoundaryMarker=false) {
        create2DGrid(x, y, markerType, worldBoundaryMarker);
    }
    /*! Create three dimensional grid.
    Boundaries on the domain border will get the markers:
    1,2,3,4,5,6 for: left, right, bottom, top, front, back*/
    void createGrid(const RVector & x, const RVector & y, const RVector & z,
                    int markerType=0, bool worldBoundaryMarker=false){
        create3DGrid(x, y, z, markerType, worldBoundaryMarker);
    }

    void createHull_(const Mesh & mesh);

    /*! Create 3D mesh with 3D boundary elements from this 2D mesh cells.
    Increase mesh dimension. Mesh should contain 2D cells. */
    Mesh createHull() const;

    /*! Create and copy global H2 mesh of this mesh.*/
    Mesh createH2() const;

    /*! Create and copy global P2 mesh of this mesh.*/
    Mesh createP2() const;

    /*! Create a partly mesh without cells from mesh, based on a vector of ptrs to boundaries */
    void createMeshByCells(const Mesh & mesh, const std::vector < Cell * > & cells);

    /*! Create a partly mesh without cells from mesh, based on a vector of ptrs to boundaries */
    void createMeshByBoundaries(const Mesh & mesh, const std::vector < Boundary * > & bounds);

    /*! Create a new mesh that is a part from this mesh, based on cell-ids */
    Mesh createMeshByCellIdx(const IndexArray & idxList);

    /*! Create a partly mesh from mesh, based on cell-ids */
    void createMeshByCellIdx(const Mesh & mesh, const IndexArray & idxList);

    /*! Create a partly mesh from mesh, based on meshs attributes.
    For a single attribute set to to 0, for unlimited set to to -1 */
    void createMeshByMarker(const Mesh & mesh, int from, int to=-1);

    /*! Syntactic sugar to extract a part of the mesh based on cells.*/
    Mesh createSubMesh(const std::vector< Cell * > & cells) const;

    /*! Syntactic sugar to extract a part of the mesh based on boundaries.*/
    Mesh createSubMesh(const std::vector< Boundary * > & bounds) const;

    /*! Syntactic sugar to extract a part of the mesh based on
    nodes with associated cells and boundaries.*/
    Mesh createSubMesh(const std::vector< Node * > & nodes) const;

    //** end creation stuff

    //! Show some infos
    void showInfos(){ std::cout << *this << std::endl; }

    const std::vector< Node * > & nodes() const { return nodeVector_; }

    const std::vector< Cell * > & cells() const { return cellVector_; }

    /*! Return const reference to all boundaries*/
    const std::vector< Boundary * > & boundaries() const { return boundaryVector_; }

    /*! Return vector of node from index list */
    std::vector< Node * > nodes(const IndexArray & ids) const;

    /*! Return a vector of nodes ptrs matching BVector b.*/
    std::vector< Node * > nodes(const BVector & b) const;

    /*! Return vector of cells from index list */
    std::vector< Cell * > cells(const IndexArray & ids) const;

    /*! Return a vector of cells ptrs matching BVector b.*/
    std::vector< Cell * > cells(const BVector & b) const;

    /*! Return vector of boundaries from index list */
    std::vector< Boundary * > boundaries(const IndexArray & ids) const;

    /*! Return a vector of boundary ptrs matching BVector b.*/
    std::vector< Boundary * > boundaries(const BVector & b) const;

    Index nodeCount(bool withSecNodes=false) const;
    Node & node(Index i) const;
    Node & node(Index i);

    inline Index secondaryNodeCount() const { return secNodeVector_.size(); }
    Node & secondaryNode(Index id) const;
    Node & secondaryNode(Index id);

    /*!Return ids for all nodes. Optionally including for secondary ndoes.*/
    IndexArray nodeIDs(bool withSecNodes=false) const;
    /*!Set all node ids.
    Size of IndexArray indicated if it should be set the secondary nodes too.
    */
    void setNodeIDs(IndexArray & ids);

    Index cellCount() const { return cellVector_.size(); }
    Cell & cell(Index i) const;
    Cell & cell(Index i);

    Index boundaryCount() const { return boundaryVector_.size(); }
    Boundary & boundary(Index i) const;
    Boundary & boundary(Index i);

    /*! Return a vector of all node positions */
    PosVector positions(bool withSecNodes=false) const;

    /*! Return a vector of node positions for an index vector */
    PosVector positions(const IndexArray & idx) const;

    /*! Return a vector of all cell center positions*/
    PosVector cellCenters() const;
    PosVector cellCenter() const { return cellCenters(); }

    /*! Return a vector of all center positions for all boundaries */
    PosVector boundaryCenters() const;

    /*! Return the reference to a RVector of all cell sizes. Cached for static geometry.*/
    const RVector & cellSizes() const;

    /*! Return the reference to a RVector of all boundary sizes. Cached for static geometry. */
    const RVector & boundarySizes() const;

    /*! Return the reference to the vector of scaled normal directions for each boundary.
     * Cached for static geometry and will be build on first call. Not thread safe, perhaps not python GC safe.
     Return \f$ \{A_i \vec{n}_i\} \forall i = [0..N_B]\f$.
     Where \f$ A_i\f$ is the size and \f$ \vec{n}_i\f$ the normal direction for the i-th boundary.
     If you want to use this, i.e. for the calculation of inside or outside flow through the boundary, you need to recognize the orientation of this boundary to the cell the flow goes into or comes from.
     For the left cell neighbor the normal direction should be always the outer normal.*/
    const PosVector & boundarySizedNormals() const;

    /*! Set the marker to all boundaries in index array. */
    void setBoundaryMarkers(const IndexArray & ids, int marker);

    /*! Set all cell marker the values in attribute. */
    void setBoundaryMarkers(const IVector & marker);

    /*! Return a vector of all boundary marker */
    IVector boundaryMarkers() const;

    /*! Return a vector of all node marker */
    IVector nodeMarkers() const;

    /*! Return an index vector of all nodes that match the marker */
    IndexArray findNodesIdxByMarker(int marker) const;

//     /*! Return an index list of all nodes that match the marker */
//     std::list < Index > findListNodesIdxByMarker(int marker) const;

    /*! Return a vector of boundary ptrs with the boundary marker equal marker.*/
    std::vector < Boundary * > findBoundaryByMarker(int marker) const;

    /*! Return a vector of boundary ptrs with the boundary marker between [from and to). \n
        for to equal open end set to = MAX_INT */
    std::vector < Boundary * > findBoundaryByMarker(int from, int to) const;

    /*! Return ptr to the cell that match position pos, counter holds amount of touch tests.
        Searching is done first by nearest-neighbor-kd-tree search,
        followed by slope-search if extensive is set. Return NULL if no cell can be found. */
    Cell * findCell(const RVector3 & pos, size_t & counter, bool extensive) const ;

    /*! Shortcut for \ref findCell(const RVector3 & pos, size_t & counter, bool extensive) */
    Cell * findCell(const RVector3 & pos, bool extensive=true) const {
        size_t counter; return findCell(pos, counter, extensive); }

    /*! Return the index to the node of this mesh with the smallest distance to pos. */
    Index findNearestNode(const RVector3 & pos);

    /*! Return vector of cell ptrs with marker match the range [from .. to). \n
        For single marker match to is set to 0, for open end set to = -1 */
    std::vector < Cell * > findCellByMarker(int from, int to=0) const;

    /*! Return vector of cell ptrs with attribute match the range [from .. to). \n
        For single attribute match to is set to 0.0, for open end set to = -1.0 */
    std::vector < Cell * > findCellByAttribute(double from, double to=0.0) const;

    /*! Return vector of cells that are intersected with a given ray from start
     * to end. Intersecting positions, i.e., the travel path are stored in pos.
     Note this will not yet check if the ray lies completely along a boundary.
     This will probably fail and need to be implemented.
     */
    std::vector < Cell * > findCellsAlongRay(const RVector3 & start,
                                             const RVector3 & dir,
                                             PosVector & pos) const;
    //** end get infos stuff

    //** start mesh modification stuff

    /*! Prolongate the empty (lower than \ref TOLERANCE.) cell values in vals
     * from its neighboring cells.
     * This function is called recursively until all zero-attribute-values in
     * vals are filled with an attribute greater than Zero.
     * RVector vals need to be of size \ref cellCount().
     * If Background is unequal -9e99 all empty values will be set to background.
     */
    void prolongateEmptyCellsValues(RVector & vals, double background=-9e99) const;

    void recountNodes();

    void sortNodes(const IndexArray & perm);

    /*! Return true if createNeighborInfos is called once */
    inline bool neighborsKnown() const { return neighborsKnown_; }

    /*! Remove from each boundary the ptr to the corresponding left and right cell*/
    void cleanNeighborInfos();

    /*! Search and set to each boundary the corresponding left and right cell.*/
    void createNeighborInfos(bool force=false);

    /*! Create and store boundaries and neighboring information for this cell.*/
    void createNeighborInfosCell_(Cell *c);

    void relax();

    /*! Smooth the mesh via moving all free nodes into the average of all neighboring nodes. Repeat this smoothIteration times. There is currently only this smoothFunction. EdgeSwapping is deactivated.*/
    void smooth(bool nodeMoving=true, bool edgeSwapping=true, uint smoothFunction=1, uint smoothIteration=10);

    /*! Scales the mesh with \ref RVector3 s.
    Returns a reference to the mesh (no copy).*/
    Mesh & scale(const RVector3 & s);

    /*! Scales the mesh with s. Shortcut for scale(RVector3(s,s,s))
    Returns a reference to the mesh (no copy).*/
    Mesh & scale(const double & s){ return scale(RVector3(s, s, s));}

    /*! Translates the mesh with \ref RVector3 t.
    Returns a reference to the mesh (no copy).*/
    Mesh & translate(const RVector3 & t);

    /*! Rotates the mesh the with \ref RVector3 r, r in radian.
    If you want to rotate in degree, use \ref degToRad(const RVector3 & deg).
    Returns a reference to the mesh (no copy).*/
    Mesh & rotate(const RVector3 & r);

    /*! Apply a 4x4 transformation matrix to the whole mesh.
    Returns a reference to the mesh (no copy).*/
    Mesh & transform(const RMatrix & mat);

    /*! Apply deformation epsilon to all nodes. Optional magnify the deformation.
    Returns a reference to the mesh (no copy).*/
    Mesh & deform(const R3Vector & eps, double magnify=1.0);

    /*! Apply deformation epsilon (with squeezed array) to all nodes.
    Optional magnify the deformation.
    Returns a reference to the mesh (no copy).*/
    Mesh & deform(const RVector & eps, double magnify=1.0);

    /*! Swap coordinate i with j for i and j lower then dimension of the mesh.
    Returns a reference to the mesh (no copy).*/
    void swapCoordinates(Index i, Index j);

    //** end mesh modification stuff
    //** start I/O stuff

    /*! Save mesh to file.
     * The file format is determined by the file suffix.
     * Supported formats are: \ref Binary, \ref VTK, \ref VTU.
     * Default is \ref Binary, the file becomes the suffix .bms.
     * If the file name ends with .vtk or .vtu, the mesh will be saved in
     * Visualization Toolkit format.
     * The filename is returned. */
    std::string save(const std::string & fileName) const;

    /*! Load Mesh from file and try to import fileformat regarding file suffix.
    * If createNeighborInfos is set, the mesh is checked for consistency and
    * missing boundaries will be created. */
    void load(const std::string & fileName, bool createNeighbors=true);

    /*! Serialize mesh into a string. Same structure like save.*/
    //std::string serialize() const;
    ByteBuffer serialize() const;
    /*! Deserialize mesh from a string. Same like loading a mesh.*/
    void deserialize(const ByteBuffer & m);

    void writeToStream(std::ostream & out) const;

    void readFromStream(std::istream & in);

    void importMod(const std::string & fileName);

    void importVTK(const std::string & fbody);

    void importVTU(const std::string & fbody);

    /*! Import Ascii STL as 3D mesh and save triangles as \ref Boundary Faces.
    Node positions can be snap to a tolerance.*/
    void importSTL(const std::string & fileName, bool isBinary=false,
                   double snap=1e-3);

    int exportSimple(const std::string & fbody, const RVector & data) const ;

    /*! Very simple export filter. Write to file fileName:
        x_0 y_0 [z_0] data1_0 [data2]_0 \n
        x_n y_n [z_n] data1_n [data2]_n \n
        with n = number of cells, and x y [z] == cell center
    */
    int exportMidCellValue(const std::string & fileName, const RVector & data1, const RVector & data2) const ;
    // no default arg here .. pygimli@win64 linker bug
    int exportMidCellValue(const std::string & fileName, const RVector & data1) const{
        RVector data2(0);
        return exportMidCellValue(fileName, data1, data2);
    }

    void exportVTK(const std::string & fbody,
                   const std::map< std::string, RVector > & data,
                   const PosVector & vec,
                   bool writeCells=true) const;

    void exportVTK(const std::string & fbody,
                   const std::map< std::string, RVector > & data,
                   bool writeCells=true) const;

    /*! Export mesh and whole data map */
    void exportVTK(const std::string & fbody, bool writeCells=true) const;

    /*! Export mesh and whole data map and vector data in vec*/
    void exportVTK(const std::string & fbody,
                   const PosVector & vec,
                   bool writeCells=true) const;

    /*! Export mesh with one additional array that will called 'arr' */
    void exportVTK(const std::string & fbody, const RVector & arr) const;

    void readVTKPoints_(std::fstream & file, const std::vector < std::string > & row);
    void readVTKCells_(std::fstream & file, const std::vector < std::string > & row);
    void readVTKScalars_(std::fstream & file, const std::vector < std::string > & row);
    void readVTKPolygons_(std::fstream & file, const std::vector < std::string > & row);

    /*! Export the mesh in filename using vtu format:
    Visualization Toolkit Unstructured Points Data (http://www.vtk.org)
    Set binary to true writes the data content in binary format.
    The file suffix .vtu will be added or substituted if .vtu or .vtk is found.
    \ref data, cell.markers and cell.attribute will be exported as data. */
    void exportVTU(const std::string & filename, bool binary = false) const ;

    /*! Export the boundary of this mesh in vtu format: Visualization Toolkit Unstructured Points Data (http://www.vtk.org) Set Binary to true writes the datacontent in binary format. The file suffix .vtu will be added or substituted if .vtu or .vtk is found. */
    void exportBoundaryVTU(const std::string & fbody, bool binary = false) const ;

    /*! Internal function for exporting VTU */
    void addVTUPiece_(std::fstream & file, const Mesh & mesh,
                        const std::map < std::string, RVector > & data) const;

    void exportAsTetgenPolyFile(const std::string & filename);
    //** end I/O stuff

    /*! All outer boundaries, i.e., all boundaries with only one cell (the left) need to be sorted that the norm vector shows outside the mesh. */
    void fixBoundaryDirections();

    void addData(const std::string & name, const CVector & data){
        this->addData(name+"_Re", real(data));
        this->addData(name+"_Im", imag(data));
    }

    /*! Add data to the mesh that will be saved with by using the binary mesh
     * format v.2. or exported with the appropriate name in VTK format,
     * if the size of data equals the amount of nodes, cells or boundaries.
     */
    void addData(const std::string & name, const RVector & data);

    /*! Add data to the mesh that will be saved with by using the binary mesh
     * format v.2. or exported with the appropriate name in VTK format,
     * if the size of data equals the amount of nodes, cells or boundaries.
     */
    void addData(const std::string & name, const PosVector & data){
        this->addData(name+"_x", x(data));
        this->addData(name+"_y", y(data));
        this->addData(name+"_z", z(data));
    }

    /*! Return the data with a given name.
     * If there is no such data an exception is thrown.*/
    RVector data(const std::string & name) const;

    /*! Return True if date with such a name exists.*/
    bool haveData(const std::string & name) const {
        return dataMap_.count(name) > 0;
    }

    /*! Return the full data map read only. */
    const std::map< std::string, RVector > & dataMap() const {
        return this->dataMap_;
    }
    /*! Replace the datamap by m */
    void setDataMap(const std::map< std::string, RVector > m) {
        this->dataMap_ = m;
    }
    /*! Print data map info.*/
    void dataInfo() const;

    /*! Empty the data map.*/
    void clearData();

    /*! Set the comment for VTK Ascii export headline.*/
    void setCommentString(const std::string & commentString) {commentString_ = commentString;}

    /*! Return comment for VTK Ascii export headline.*/
    const std::string & commentString() const {return commentString_;}

    void mapCellAttributes(const std::map < float, float > & aMap);

    //void mapParameterToAttribute(const std::vector< int > & cellMapIndex);
    /*! Change all boundary marker that match bMap.first to bMap.second. */
    void mapBoundaryMarker(const std::map < int, int > & aMap);

    /*! Set all cell attributes to the values in vector attribute.*/
    void setCellAttributes(const RVector & attribute);

    /*! Set all cell attributes to the scalar value: attribute */
    void setCellAttributes(double attribute);

    /*! Return a RVector of all cell attributes */
    RVector cellAttributes() const;

    /*! Set the cell marker of all indices in ids to marker. */
    void setCellMarkers(const IndexArray & ids, int marker);

    /*! Set all cell marker the values in attribute. */
    void setCellMarkers(const IVector & marker);

    /*! Set all cell marker the values in attribute (casting to int)*/
    void setCellMarkers(const RVector & attribute);

    /*! Return a vector of all cell marker */
    IVector cellMarkers() const;

    //** probably deprecated 20191120
    // double xmin() const { findRange_(); return minRange_[0]; }
    // double ymin() const { findRange_(); return minRange_[1]; }
    // double zmin() const { findRange_(); return minRange_[2]; }
    // double xmax() const { findRange_(); return maxRange_[0]; }
    // double ymax() const { findRange_(); return maxRange_[1]; }
    // double zmax() const { findRange_(); return maxRange_[2]; }

    // better use your bounding box
    double xMin() const { findRange_(); return minRange_[0]; }
    double yMin() const { findRange_(); return minRange_[1]; }
    double zMin() const { findRange_(); return minRange_[2]; }
    double xMax() const { findRange_(); return maxRange_[0]; }
    double yMax() const { findRange_(); return maxRange_[1]; }
    double zMax() const { findRange_(); return maxRange_[2]; }


    const BoundingBox boundingBox() const { findRange_(); return BoundingBox(minRange_, maxRange_);}

    /*! Return the interpolation matrix I from all node positions
     * to the query points q. I is a (len(q) x nodeCount()) SparseMapMatrix.
     * To perform the interpolation just calculate the matrix vector product.
     * uInterpolated = I.mult(uPerNode) or uInterpolated = I * uPerNode */
    RSparseMapMatrix interpolationMatrix(const PosVector & q);
    std::vector < RSparseMapMatrix > interpolationMatrix(
                                        const std::vector< PosVector > & q);

    /*! Inplace version of \ref interpolationMatrix(const PosVector & q) */
    void interpolationMatrix(const PosVector & q, RSparseMapMatrix & I);

    /*! Inplace version of \ref interpolationMatrix(const PosVector & q) */
    void interpolationMatrix(const std::vector < PosVector > & q,
                             std::vector < RSparseMapMatrix > & I);


    /*! Return the reference to the matrix for cell value to boundary value interpolation matrix. */
    RSparseMapMatrix & cellToBoundaryInterpolation() const;

    /*!Return the divergence for each cell of a given vector field for each
     * boundary.
     * The divergence is calculated by simple 1 point boundary integration
     * over each cell.
     * \f$ d(cell) = \sum_boundaries V(boundary center) \cdot n(boundary)\f$
     * Higher order integration needs to be implemented.
     * Contact the author if you need this.*/
    RVector divergence(const PosVector & V) const;

    /*! Interpolate boundary based values to cell based gradients. */
    PosVector boundaryDataToCellGradient(const RVector & boundaryData) const;

    /*! Interpolate cell based values to boundary based gradients. */
    PosVector cellDataToBoundaryGradient(const RVector & cellData) const;

    /*! Interpolate cell based values to boundary based gradients with a given cell Gradient.*/
    PosVector cellDataToBoundaryGradient(const RVector & cellData,
                                        const PosVector & cellGradient) const;

    /*! Add a region marker for tetgen or triangle creation if the mesh
     *is a PLC, if area is < 0 a hole is added. */
    void addRegionMarker(const RVector3 & pos, int marker, double area=0);
    void addRegionMarker(const RegionMarker & reg);

    const RegionMarkerList & regionMarkers() const { return regionMarker_; }

    /*! Return the pointer to region marker with the marker is i or throws
    an exception of there is no such marker.*/
    RegionMarker * regionMarker(SIndex i);

    /*! Add a hole marker for tetgen or triangle creation if the mesh
     * is a PLC */
    void addHoleMarker(const RVector3 & pos);

    /*!Return read only reference for all defined hole regions. */
    const HoleMarkerList & holeMarker() const { return holeMarker_; }

    Index hash() const;

protected:
    void copy_(const Mesh & mesh);

    void findRange_() const ;

    /*!Ensure is geometry check*/
    Node * createNodeGC_(const RVector3 & pos, int marker);

    Node * createNode_(const RVector3 & pos, int marker);

    Node * createSecondaryNode_(const RVector3 & pos);

    template < class B > Boundary * createBoundary_(
        std::vector < Node * > & nodes, int marker, int id){

        if (id == -1) id = boundaryCount();
        boundaryVector_.push_back(new B(nodes));
        boundaryVector_.back()->setMarker(marker);
        boundaryVector_.back()->setId(id);
        ShapeFunctionCache::instance().shapeFunctions(*boundaryVector_.back());
        return boundaryVector_.back();
    }

    template < class B > Boundary * createBoundaryChecked_(
        std::vector < Node * > & nodes, int marker, bool check=true){

        if (!check) return createBoundary_< B >(nodes, marker, boundaryCount());

        Boundary * b = findBoundary(nodes);
        if (!b) {
            b = createBoundary_< B >(nodes, marker, boundaryCount());
        } else {
            if (marker != 0) b->setMarker(marker);
        }
        return b;
    }

    template < class C > Cell * createCell_(
        std::vector < Node * > & nodes, int marker, int id){

        if (id == -1) id = cellCount();
        cellVector_.push_back(new C(nodes));
        cellVector_.back()->setMarker(marker);
        cellVector_.back()->setId(id);
        ShapeFunctionCache::instance().shapeFunctions(*cellVector_.back());
        return cellVector_.back();
    }

//    Node * createRefinementNode_(Node * n0, Node * n1, SparseMapMatrix < Node *, Index > & nodeMatrix);
    Node * createRefinementNode_(Node * n0, Node * n1, std::map< std::pair < Index, Index >, Node * > & nodeMatrix);

    void createRefined_(const Mesh & mesh, bool p2, bool r2);

    Cell * findCellBySlopeSearch_(const RVector3 & pos, Cell * start, size_t & count, bool tagging) const;

    void fillKDTree_() const;

    std::vector< Node * >     nodeVector_;
    std::vector< Node * >     secNodeVector_;
    std::vector< Boundary * > boundaryVector_;
    std::vector< Cell * >     cellVector_;

    uint dimension_;

    mutable RVector3 minRange_;
    mutable RVector3 maxRange_;
    mutable bool rangesKnown_;

    bool neighborsKnown_;

    mutable KDTreeWrapper * tree_;

    /*! A static geometry mesh caches geometry informations. */
    bool staticGeometry_;
    bool isGeometry_; // mesh is marked as PLC
    mutable RVector cellSizesCache_;
    mutable RVector boundarySizesCache_;
    mutable PosVector boundarySizedNormCache_;

    mutable RSparseMapMatrix * cellToBoundaryInterpolationCache_;

    bool oldTet10NumberingStyle_;

    std::map< std::string, RVector > dataMap_;
    std::string commentString_;

    // for PLC creation
    RegionMarkerList regionMarker_;
    HoleMarkerList holeMarker_;

}; // class Mesh

} // namespace GIMLI;

#endif // _GIMLI_MESH__H
