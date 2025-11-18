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

#ifndef _GIMLI_MESHENTITIES__H
#define _GIMLI_MESHENTITIES__H

#include "gimli.h"

#include "vector.h"
#include "matrix.h"
#include "elementmatrix.h"
#include "baseentity.h"

#include "polynomial.h"

#include <vector>
#include <set>

namespace GIMLI{

template < class Typname > bool lesserId(const Typname * a, const Typname * b){
    return (a->id() < b->id());
}

DLLEXPORT Boundary * findBoundary(const Node & n1);
DLLEXPORT Boundary * findBoundary(const Node & n1, const Node & n2);
DLLEXPORT Boundary * findBoundary(const Node & n1, const Node & n2, const Node & n3);

DLLEXPORT Boundary * findBoundary(const std::vector < Node * > & n);

/*! Find all boundaries that have these nodes in common. */
DLLEXPORT std::set < Boundary * >  
findBoundaries(const std::vector < Node * > & n);

/*! */
DLLEXPORT Boundary * findCommonBoundary(const Cell & c1, const Cell & c2);

DLLEXPORT Cell * findCommonCell(const std::vector < Node * > & n, bool warn = true);

/*! Functor to collect all nodes of a MeshEntity into a given std::set. */
class CollectNodeFunctor{
public:
    CollectNodeFunctor(std::set< Node * > & c):c_(& c){}

    template < class E > void operator()(E * e){
        for (Index i = 0; i < e->nodeCount(); i++) c_->insert(& e->node(i));
    }
    std::set< Node * > * c_;
};

/*! Collect all common \ref Node's for a container of MeshEntities, e.g.,
 * commonNodes(std::list< Cell > & cellList), returns all common nodes that are
 * associated to each cell at the cellList. */
template < class ContainerOfMeshEntities >
std::set< Node * > commonNodes(const ContainerOfMeshEntities & c){
    std::set< Node * > commons;

    CollectNodeFunctor gc(commons);
    for_each(c.begin(), c.end(), gc);
    return commons;
}

class DLLEXPORT RegionMarker : public RVector3{
public:
    RegionMarker(const RVector3 & pos, int marker, double area=0.0,
                 bool hole=false)
    : RVector3(pos), marker_(marker), area_(area), isHole_(hole){}

    ~RegionMarker(){}

    inline void setMarker(SIndex marker) {marker_ = marker;}
    inline int marker() const {return marker_;}

    inline void setArea(double area) {area_ = area;}
    inline double area() const {return area_;}

    inline void setPos(const Pos & pos) {copy_(pos);}

    bool isHole() const { return isHole_; }
    inline void setHole(bool h) { isHole_ = h; }

protected:
    int marker_;
    double area_;
    bool isHole_;
};

class DLLEXPORT MeshEntity : public BaseEntity {
public:

    /*! Default constructor.*/
    MeshEntity();

    /*! Default destructor.*/
    virtual ~MeshEntity();

    /*! Return the dimension for this MeshEntity. */
    virtual uint dim() const { return 0; }

    /*! Return the runtime identification for this MeshEntity. */
    virtual uint rtti() const { return MESH_MESHENTITY_RTTI; }

    /*! To separate between major MeshEntity families e.g. Cell and Boundary. */
    virtual uint parentType() const { return MESH_MESHENTITY_RTTI; }

    virtual void setNodes(const std::vector < Node * > & nodes);

    const std::vector< Node * > & nodes() const { return nodeVector_; }

    inline Node & node(uint i) {
        ASSERT_RANGE(i, 0, nodeCount()); return *nodeVector_[i];
    }
    
    inline Node & node(uint i) const {
        ASSERT_RANGE(i, 0, nodeCount()); return *nodeVector_[i];
    }

    inline uint nodeCount() const { return nodeVector_.size(); }

    inline Shape & shape() { return *shape_; }

    inline Shape & shape() const { return *shape_; }

    inline Shape * pShape() { return shape_; }

    /*! Return rst-coordinates for the i-th node. See Shape::rst. */
    RVector3 rst(uint i) const;

    /*! Return the center coordinates of this MeshEntity. */
    RVector3 center() const;

    /*! Return the size (i.e., length, area, volume) of this MeshEntity. */
    double size() const;

    virtual double attribute() const { return -1.0; }

    /*! Return IndexArray of all node ids. */
    IndexArray ids() const;

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    /*! Return a \ref RVector(n) for
     *\f$ N(L) = {N_i} i = [0, \mathrm{nodeCount()}] \f$ shape functions
     *\f$ N_n(L_1,L_2,L_3)\f$ in local coordinates \f$ L(L1, L2, L3) = L(r, s, t)\f$ */
    virtual RVector N(const RVector3 & rst) const;

    /*! Inplace variant of \ref N */
    virtual void N(const RVector3 & rst, RVector & n) const;

    /*! Return a \ref RVector of the derivation for the \f$ n=[0,\mathrm{nodeCount()}] \f$ shape functions \f$ N_n(L_1,L_2,L_3)\f$ for the local coordinate \f$ (L_1,L_2,L_3)\f$
     *   regarding to the local coordinates \f$ L_i \f$ \n
     * \f$ \frac{\partial N_n(L_1,L_2,L_3)}{\partial L_i} \f$ with may be \f$ i = 0,1,2 \f$*/
    virtual RVector dNdL(const RVector3 & rst, uint i) const;

    /*! Return the complete (n,3) matrix for all \f$ n=[0,\mathrm{nodeCount()}] \f$
     * shape functions of
     * all nodes of the current MeshEntity
     * \f$ \frac{\partial N_n(L_1,L_2,L_3)}{\partial L_i} \f$
     * with \f$ i = 0,1,2 \f$ */
    virtual RMatrix dNdL(const RVector3 & rst) const;

    /*! Interpolate a scalar field at position p for the scalar field u regarding to the shape functions of the entity.
     * \param p Cartesian coordinates (x,y,z) need to be inside, or on the boundary, of the entity.
     * \param u The field vector u need to be of size mesh.nodeCount() for the corresponding mesh.  */
    double pot(const RVector3 & p, const RVector & u) const;

    /*! Interpolate a vector field at position p for the vector field v regarding to the shape functions of the entity.
     * \param p Cartesian coordinates (x,y,z) need to be inside, or on the boundary, of the entity.
     * \param v The vector field vector v need to be of size mesh.nodeCount() for the corresponding mesh. */
    RVector3 vec(const RVector3 & p, const R3Vector & v) const;

    /*! Return gradient at position pos for field u regarding to the shape functions of the entity.
    * The field vector u need to be of size mesh.nodeCount() for the corresponding mesh.
    * The position pos, in Cartesian coordinates (x,y,z), need to be inside, or on the boundary, of the entity.  */
    RVector3 grad(const RVector3 & p, const RVector & u) const;

    friend std::ostream & operator << (std::ostream & str, const MeshEntity & c);

//    void setUxCache(const ElementMatrix < double >  & mat) const { uxCache_ = mat; }

 //   const ElementMatrix < double > & uxCache() const { return uxCache_; }

    void setUxCache(const RMatrix & mat) const { uxCache_ = mat; }

    const RMatrix & uxCache() const { return uxCache_; }

    ElementMatrix < double > & uCache(){ return uCache_; }

    ElementMatrix < double > & gradUCache(){ return gradUCache_; }

    /*! Geometry has been changed. Deletes cache.*/
    void changed();

    void addSecondaryNode(Node * n);

    void delSecondaryNode(Node * n);

    const std::vector < Node * > & secondaryNodes() const;

    /*! Return primary and secondary nodes */
    const std::vector < Node * > allNodes() const;

    Index allNodeCount() const;

    /*! Reverse node sequence order to enforce positive Jacobian determinant.
     * Please use with care! Return True if the order has been changed.*/
    virtual bool enforcePositiveDirection();

protected:
    void fillShape_();

    virtual void registerNodes_();
    virtual void deRegisterNodes_();

    virtual void registerSecNode_(Node *n);
    virtual void deRegisterSecNode_(Node *n);
    

    Shape * shape_;

    std::vector < Node * > nodeVector_;
    std::vector < Node * > secondaryNodes_;

    /*! Cache for derivation matrixes */
    mutable ElementMatrix < double > uCache_;
    mutable ElementMatrix < double > gradUCache_;

    mutable RMatrix uxCache_;

protected:
    /*! do not copy a mesh entity at all */
    MeshEntity(const MeshEntity & ent){
        THROW_TO_IMPL
    }

    /*! do not assign a mesh entity at all */
    MeshEntity & operator = (const MeshEntity & ent){
         if (this != & ent){
             THROW_TO_IMPL
         } return *this;
    }
};

//! A abstract cell
/*! Interface class for all cells. */
class DLLEXPORT Cell : public MeshEntity {
public:
    DLLEXPORT friend std::ostream & operator << (std::ostream & str, const Cell & c);

    /*! Default constructor. */
    Cell();

    /*! Construct cell from vector of nodes. */
    Cell(const std::vector < Node * > & nodes);

    /*! Default destructor. */
    virtual ~Cell();

    /*! For pygimli bindings to allow simple check*/
    bool operator==(const Cell & cell){
        return &cell == this;
    }

    virtual uint rtti() const { return MESH_CELL_RTTI; }
    virtual uint parentType() const { return MESH_CELL_RTTI; }
    virtual uint neighborCellCount() const { return 0; }
    inline uint boundaryCount() const { return neighborCellCount(); }

    void cleanNeighborInfos();

    Cell * neighborCell(const RVector & sf);

    /*! Return the direct neighbor cell corresponding to local node i.
     * The cell will be searched and stored by the virtual method
     * \ref findNeighborCell.
     * All neighboring relationships have to be initialized ones by calling
     * \ref Mesh::createNeighborInfos().
     * If no cell can be found NULL is returned. */
    inline Cell * neighborCell(uint i){ return neighborCells_[i]; }

    /*! Find neighbor cell regarding to the i-th Boundary and store them
     * in neighborCells_. */
    virtual void findNeighborCell(uint i);

    inline double attribute() const { return attribute_; }

    inline void setAttribute(double attr) { attribute_ = attr; }

    /*! DEPRECATED????
      Find the node of this cell which is in opposite position to the given boundary. Returns a pointer to the node. The boundary must be part of the cell otherwise, a NULL pointer returns. Works for triangle/edge and tetrahedron/triangleFace*/
    Node * oppositeTo(const Boundary & bound);

    /*! Find the nearest boundary to be crossed in direction to the point
     * responsible for the shape function. */
    Boundary * boundaryTo(const RVector & sf);

    /*! Return the i-th boundary. Experimental! */
    Boundary * boundary(Index i);

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i) const{
        CERR_TO_IMPL
        std::cout << rtti() << std::endl;
        std::vector < Node * > n;
        return n;
    }

protected:
    virtual void registerNodes_();
    virtual void deRegisterNodes_();
    virtual void registerSecNode_(Node *n);
    virtual void deRegisterSecNode_(Node *n);
    
    std::vector < Cell * > neighborCells_;
    double attribute_;

protected:
    /*! Don't call this class directly */
    Cell(const Cell & cell){
        std::cerr << "cell(const cell & cell)" << std::endl;
        THROW_TO_IMPL
    }

    /*! Don't call this class directly */
    Cell & operator = (const Cell & cell){
        if (this != &cell) {
            THROW_TO_IMPL
            std::cerr << "cell=cell" << std::endl;
        }
        return *this;
    }

};

class DLLEXPORT Boundary : public MeshEntity{
public:
    Boundary();
    Boundary(const std::vector < Node * > & nodes);
    virtual ~Boundary();

    virtual uint rtti() const { return MESH_BOUNDARY_RTTI; }
    virtual uint parentType() const { return MESH_BOUNDARY_RTTI; }

    /*! Return these coordinates manual until boundary coordinate transformation is done. */
    virtual RVector3 rst(uint i) const;

    /*! Normal vector of the boundary shows outside for left cell.
     * Every boundary needs a left cell for a valid mesh. */
    inline const Cell & leftCell() const { return *leftCell_; }
    inline Cell * leftCell() { return leftCell_; }

    inline const Cell & rightCell() const { return *rightCell_; }
    inline Cell * rightCell() { return rightCell_; }

    inline void setLeftCell(Cell * cell) { leftCell_  = cell; }
    inline void setRightCell(Cell * cell) { rightCell_ = cell; }

    friend std::ostream & operator << (std::ostream & str, const Boundary & e);

    /*!Return normal vector for this boundary.*/
    virtual RVector3 norm() const;

    /*!Return outer normal vector for this boundary regarding the given cell.
     *The boundary should part of this cell.*/
    virtual RVector3 norm(const Cell & cell) const;

    /*! Return true if the normal vector of this boundary shown from the cell away (outside-direction) */
    bool normShowsOutside(const Cell & cell) const;

    /*! Reverse node order to swap normal direction. */
    void swapNorm(bool withNeighbours=true);

    /*!Is the boundary is on the outside of the mesh.*/
    bool outside() const { return (leftCell_ != 0) && (rightCell_ == 0); }


protected:
    virtual void registerNodes_();
    virtual void deRegisterNodes_();
    virtual void registerSecNode_(Node *n);
    virtual void deRegisterSecNode_(Node *n);
    

    Cell *leftCell_;
    Cell *rightCell_;

protected:
    /*! Don't call this class directly */
    Boundary(const Boundary & bound){
        std::cerr << "Boundary(const Boundary & bound)" << std::endl;
        THROW_TO_IMPL
    }

    /*! Don't call this class directly */
    Boundary & operator = (const Boundary & boundary){
        if (this != &boundary) {
            std::cerr << "Assignment for boundaries not yet supported." << std::endl;
            THROW_TO_IMPL
        }
        return *this;
    }

};

class DLLEXPORT NodeBoundary : public Boundary{
public:
    NodeBoundary(Node & n1);

    NodeBoundary(const std::vector < Node * > & nodes);

    virtual ~NodeBoundary();

    virtual uint dim() const { return 1; }

    virtual uint rtti() const { return MESH_BOUNDARY_NODE_RTTI; }

    void setNodes(Node & n1);

    virtual double size() const { return 1.0; }

    friend std::ostream & operator << (std::ostream & str, const NodeBoundary & e);

    /*! Returns the normal vector for this boundary that shows outside along the
     * tangent of the left neighbouring cell which is an Edge. If there are no neighbour infos. [1.0, 0.0, 0.0] is returned. */
    virtual RVector3 norm() const;

protected:
};

class DLLEXPORT Edge : public Boundary{
public:
    Edge(Node & n1, Node & n2);

    Edge(const std::vector < Node * > & nodes);

    virtual ~Edge();

    virtual uint dim() const { return 1; }

    virtual uint rtti() const { return MESH_EDGE_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    void setNodes(Node & n1, Node & n2);

    /*! Swap edge between two triangular neighbor cells. Only defined if both neighbors are triangles. */
    int swap();

//     /*! See ref MeshEntity::N() */
//     virtual RVector N(const RVector3 & L) const;

//     void shapeFunctionsL(const RVector3 & L, RVector & funct) const;

    friend std::ostream & operator << (std::ostream & str, const Edge & e);

protected:
};

class DLLEXPORT Edge3 : public Edge{
public:

    Edge3(const std::vector < Node * > & nodes);

    virtual ~Edge3();

    virtual uint rtti() const { return MESH_EDGE3_RTTI; }

    /*! return these coordinates manual until boundary coordinate transformation is done. */
    virtual RVector3 rst(uint i) const;

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

//     /*! See MeshEntity::N() */
//     virtual RVector N(const RVector3 & L) const;
    /*
    void shapeFunctionsL(const RVector3 & L, RVector & funct) const;*/

// friend std::ostream & operator << (std::ostream & str, const Edge & e);

protected:
};


class DLLEXPORT TriangleFace : public Boundary{
public:
    TriangleFace(Node & n1, Node & n2, Node & n3);

    TriangleFace(const std::vector < Node * > & nodes);

    virtual ~TriangleFace();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_TRIANGLEFACE_RTTI; }

    void setNodes(Node & n1, Node & n2, Node & n3);

    friend std::ostream & operator << (std::ostream & str, const TriangleFace & e);

    /*! this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:

     // /*! Don't call this class directly */
    // TriangleFace(const TriangleFace & bound){
        // std::cerr << "TriangleFace(const Boundary & bound)" << std::endl;
    // }

    // /*! Don't call this class directly */
    // TriangleFace & operator = (const TriangleFace & boundary){
        // if (this != &boundary) {
            // std::cerr << "TriangleFace boundary=boundary" << std::endl;
        // }
        // return *this;
    // }
};

class DLLEXPORT Triangle6Face : public TriangleFace{
public:
    Triangle6Face(const std::vector < Node * > & nodes);

    ~Triangle6Face();

    virtual uint rtti() const { return MESH_TRIANGLEFACE6_RTTI; }

    /*! this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    /*! this method need refactoring */
    virtual RVector3 rst(uint i) const;

protected:
};

class DLLEXPORT QuadrangleFace : public Boundary{
public:
    QuadrangleFace(Node & n1, Node & n2, Node & n3, Node & n4);

    QuadrangleFace(const std::vector < Node * > & nodes);

    virtual ~QuadrangleFace();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_QUADRANGLEFACE_RTTI; }

    /*! this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    void setNodes(Node & n1, Node & n2, Node & n3, Node & n4);

    friend std::ostream & operator << (std::ostream & str, const TriangleFace & e);

protected:

    QuadrangleFace(const QuadrangleFace & quad){
        std::cerr << "QuadrangleFace(const QuadrangleFace & quad)" << std::endl;
    }
};

class DLLEXPORT Quadrangle8Face : public QuadrangleFace{
public:

    Quadrangle8Face(const std::vector < Node * > & nodes);

    virtual ~Quadrangle8Face();

    virtual uint rtti() const { return MESH_QUADRANGLEFACE8_RTTI; }

    /*! this method need refactoring */
    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    /*! this method need refactoring */
    virtual RVector3 rst(uint i) const;

protected:

private:

};

class DLLEXPORT PolygonFace : public Boundary {
    /*! */
public:
    typedef Pos HoleMarker;
    typedef std::vector < Pos >  HoleMarkerList;

    PolygonFace(const std::vector < Node * > & nodes);

    ~PolygonFace();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_POLYGON_FACE_RTTI; }

    /*! Insert node into the polygon. Node needs to touch the polygon.
    The node will be inserted in the nodeList or as secondary node if its
    not on an edge.*/
    void insertNode(Node * node, double tol=TOLERANCE);

    /*! Insert nodes for a subpolygon.
    All nodes regarding the parent mesh and need to be inside the face.
    If marked as hole the boundary does not become part if nodes boundary set.*/
    void addSubface(const std::vector < Node * > & nodes, bool isHole=false);

    Index subfaceCount() const {return this->subfaces_.size();}

    const std::vector < Node * > & subface(Index i) const;

    /*! Add a hole marker for tetgen or triangle creation if the mesh
     * is a PLC */
    void addHoleMarker(const RVector3 & pos);

    void delHoleMarker(const RVector3 & pos);

    /*!Return read only reference for all defined hole regions. */
    const HoleMarkerList & holeMarkers() const;

    /*!Return reference to all defined hole markers. */
    HoleMarkerList & holeMarkers(){ return holeMarker_;}

protected:
    std::vector < std::vector < Node * > > subfaces_;
    HoleMarkerList holeMarker_;
private:

};

class DLLEXPORT EdgeCell : public Cell {
public:
    EdgeCell(Node & n1, Node & n2);

    EdgeCell(const std::vector < Node * > & nodes);

    virtual ~EdgeCell();

    virtual uint dim() const { return 1; }

    virtual uint rtti() const { return MESH_EDGE_CELL_RTTI; }

    virtual uint neighborCellCount() const { return 2; }

//     virtual void findNeighborCell(uint id);

    void setNodes(Node & n1, Node & n2);

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    /*! Index : NodeBoundary
     * 0 : 1 \n
     * 1 : 0
     */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

    friend std::ostream & operator << (std::ostream & str, const EdgeCell & t);

protected:
};

/*! count: 0-1, 2(0-1) */
class DLLEXPORT Edge3Cell : public EdgeCell {
public:
    Edge3Cell(const std::vector < Node * > & nodes);

    virtual ~Edge3Cell();

    virtual uint rtti() const { return MESH_EDGE3_CELL_RTTI; }

    virtual RVector3 rst(uint i) const;

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

//! Triangle
/*! Triangle

Node direction:
    2
   / \
  /   \
 /     \
0-------1
*/
class DLLEXPORT Triangle : public Cell {
public:
    Triangle(Node & n1, Node & n2, Node & n3);

    Triangle(const std::vector < Node * > & nodes);

    virtual ~Triangle();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_TRIANGLE_RTTI; }

    virtual uint neighborCellCount() const { return 3; }

//     virtual void findNeighborCell(uint i);

    void setNodes(Node & n1, Node & n2, Node & n3);

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    friend std::ostream & operator << (std::ostream & str, const Triangle & t);

    /*! Index : Edge Boundaries
     * 0 : 1-2 \n
     * 1 : 2-0 \n
     * 2 : 0-1 \n
     */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:
};

//! Triangle6
/*! Triangle6

Node direction:
    2
   / \
  5   4
 /     \
0---3---1
*/
class DLLEXPORT Triangle6 : public Triangle {
public:
    Triangle6(const std::vector < Node * > & nodes);

    virtual ~Triangle6();

    virtual uint rtti() const { return MESH_TRIANGLE6_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

//! Quadrangle
/*! Quadrangle

Node direction:
   3-----2
  /     /
 /     /
0-----1

Neighborship relations:

Neighbor Nr, on Boundary a->b
    0           2->3
    1           3->0
    2           0->1
    3           1->2
*/
class DLLEXPORT Quadrangle : public Cell {
public:
    Quadrangle(Node & n1, Node & n2, Node & n3, Node & n4);

    Quadrangle(const std::vector < Node * > & nodes);

    virtual ~Quadrangle();

    virtual uint dim() const { return 2; }

    virtual uint rtti() const { return MESH_QUADRANGLE_RTTI; }

    void setNodes(Node & n1, Node & n2, Node & n3, Node & n4);

    virtual uint neighborCellCount() const { return 4; }

//     virtual void findNeighborCell(uint i);

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    friend std::ostream & operator << (std::ostream & str, const Quadrangle & t);

    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:
};

//! Quadrangle8 for serendipity type
/*! Quadrangle serendipity type with 8 nodes for quadratic shape functions
 *
 * Node direction: \n
 *     3---6---2  \n
 *    /       /   \n
 *   7       5    \n
 *  /       /     \n
 * 0---4---1      \n
 *
 * count: 0-1-2-3, 4(0-1), 5(1-2), 6(2-3), 7(3-0)
*/
class DLLEXPORT Quadrangle8 : public Quadrangle {
public:
    Quadrangle8(const std::vector < Node * > & nodes);

    virtual ~Quadrangle8();

    virtual uint rtti() const { return MESH_QUADRANGLE8_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

//! A Tetrahedron
/*! A Tetrahedron

Node direction:

3           \n
| 2         \n
|/          \n
0-----1     \n

Neighborship relations:
Boundary normal shows outside .. so the boundary left neighbor is this cell

Neighbor Nr., on Boundary a-b-c. Boundary to neighbor cell is opposite to NodeNr.
    0           1-2-3     le -- view from outer
    1           2-0-3     ri -- view from inner
    2           0-1-3     le -- view from outer
    3           0-2-1     ri -- view from inner
*/
static const uint8 TetrahedronFacesID[4][3] = {
    {1, 2, 3},
    {2, 0, 3},
    {0, 1, 3},
    {0, 2, 1}
};

//! A Tetrahedron
/*!
 * A Tetrahedron for linear shape functions
*/
class DLLEXPORT Tetrahedron : public Cell {
public:
    Tetrahedron(Node & n1, Node & n2, Node & n3, Node & n4);

    Tetrahedron(const std::vector < Node * > & nodes);

    virtual ~Tetrahedron();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_TETRAHEDRON_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    void setNodes(Node & n1, Node & n2, Node & n3, Node & n4);

    virtual uint neighborCellCount() const { return 4; }

//     virtual void findNeighborCell(uint i);

    friend std::ostream & operator << (std::ostream & str, const Tetrahedron & t);

    /*! Index : Triangle Boundaries
     * 0:
     *
     *
     */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:

    /*! Don't call this class directly */
    Tetrahedron(const Tetrahedron& cell){ std::cerr << "Tetrahedron cell(const cell & cell)" << std::endl; }

    /*! Don't call this class directly */
    Tetrahedron & operator = (const Tetrahedron & cell){
        if (this != &cell) {
            std::cerr << "Tetrahedron cell=cell" << std::endl;
        }
        return *this;
    }

};

//*! VTK,Flaherty,Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //
static const uint8 Tet10NodeSplit[10][2] = {
    {0,0},{1,1},{2,2},{3,3},
    {0,1},{1,2},{2,0},
    {0,3},{1,3},{2,3}
};

//*! Zienkiewicz count: 1-2-3-4, 5(1-2), 6(1-3), 7(1-4), 8(2-3), 9(3-4), 10(4-2)* //
static const uint8 Tet10NodeSplitZienk[10][2] = {
    {0,0},{1,1},{2,2},{3,3},
    {0,1},{0,2},{0,3},
    {1,2},{2,3},{3,1}
};

class DLLEXPORT Tetrahedron10 : public Tetrahedron {
public:
    Tetrahedron10(const std::vector < Node * > & nodes);

    virtual ~Tetrahedron10();

    virtual uint rtti() const { return MESH_TETRAHEDRON10_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

//! A Hexahedron
/*! A Hexahedron

Node direction:

  7------6  \n
 /|     /|  \n
4------5 |  \n
| 3----|-2  \n
|/     |/   \n
0------1    \n

Neighborship relations:
Boundary normal shows outside .. so the boundary left neighbor is this cell

Neighbor Nr, on Boundary a-b-c-d
    0           1-2-6-5  // le
    1           2-3-7-6  // ri
    2           3-0-4-7  // ri
    3           0-1-5-4  // le
    4           4-5-6-7  // le
    5           0-3-2-1  // ri

T.~Apel and N.~Düvelmeyer, Transformation of Hexahedral Finite Element Meshes into Tetrahedral Meshes According to Quality Criteria,
Computing Volume 71, Number 4 / November, 2003, DOI 10.1007/s00607-003-0031-5, Pages   293-304
5-Tet-split: type 6(2) 1-4-5-6, 3-7-6-4, 1-4-0-3, 1-2-3-6, 1-6-4-3
6-Tet-split: type 1    0-1-2-6, 0-2-3-6, 0-1-6-5, 0-4-5-6, 0-3-7-6, 0-4-6-7

*/

static const uint8 HexahedronFacesID[6][4] = {
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {3, 0, 4, 7},
    {0, 1, 5, 4},
    {4, 5, 6, 7},
    {0, 3, 2, 1}
};

class DLLEXPORT Hexahedron: public Cell {
public:
    Hexahedron(const std::vector < Node * > & nodes);

    virtual ~Hexahedron();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_HEXAHEDRON_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    virtual uint neighborCellCount() const { return 6; }

    friend std::ostream & operator << (std::ostream & str, const Hexahedron & t);

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:
};

static const uint8 Hexahedron20FacesID[6][8] = {
    {0,1,5,4,8,17,12,16},
    {1,2,6,5,9,18,13,17},
    {2,3,7,6,10,19,14,18},
    {3,0,4,7,11,16,15,19},
    {0,3,2,1,11,10,9,8},
    {4,5,6,7,12,13,14,15},
};

static const uint8 Hex20NodeSplit[20][2] = {
    {0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},
    {0,1},{1,2},{2,3},{3,0},
    {4,5},{5,6},{6,7},{7,4},
    {0,4},{1,5},{2,6},{3,7}
};

//! A Hexahedron with 20 nodes
/*! A Hexahedron with 20 nodes for quadratic base functions of serendipity style

Node direction:

          7-----14------6  \n
         /|            /|  \n
    t   / |    s      / |  \n
    | 15 19   /     13 18  \n
    | /   |  /      /   |  \n
    |/    | /      /    |  \n
    4-----12------5     |  \n
    |     3-----10|-----2  \n
    |    /        |    /   \n
   16   /         17  /    \n
    | 11          |  9     \n
    | /           | /      \n
    |/            |/       \n
    0------8-----1-------r \n

*/
class DLLEXPORT Hexahedron20: public Hexahedron {
public:
    Hexahedron20(const std::vector < Node * > & nodes);

    virtual ~Hexahedron20();

    virtual uint rtti() const { return MESH_HEXAHEDRON20_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:
};

//! Triangular prism
/*!
 * A Triangular prism is a three-sided prism. Equivalently, it is a pentahedron of which two faces are parallel.
 * Node direction:

  5    \n
 /|\   \n
3---4  \n
| 2 |  \n
|/ \|  \n
0---1  \n

*/

static const uint8 TriPrismFacesID[5][4] = {
    {1, 2, 5, 4},        // r
    {2, 0, 3, 5},        // r
    {0, 1, 4, 3},        // l
    {3, 4, 5, 255},      // l
    {0, 2, 1, 255},      // r
};


class DLLEXPORT TriPrism : public Cell {
public:
    TriPrism(const std::vector < Node * > & nodes);

    virtual ~TriPrism();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_TRIPRISM_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    virtual uint neighborCellCount() const { return 5; }

    friend std::ostream & operator << (std::ostream & str, const Hexahedron & t);

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:
};

static const uint8 Prism15NodeSplit[15][2] = {
    {0,0},{1,1},{2,2},{3,3},{4,4},{5,5},
    {0,1},{1,2},{2,0},
    {3,4},{4,5},{5,3},
    {0,3},{1,4},{2,5}
};

//! Triangular15 prism
/*!
 * A Triangular prism with 15 Nodes is a three-sided prism for quadratic base functions.

        5      \n
       / \     \n
     11 14 10  \n
     /  |  \   \n
    3---9---4  \n
    |   2   |  \n
    |  / \  |  \n
   12 8   7 13 \n
    |/     \|  \n
    0---6---1  \n

*/

class DLLEXPORT TriPrism15 : public TriPrism {
public:
    TriPrism15(const std::vector < Node * > & nodes);

    virtual ~TriPrism15();

    virtual uint rtti() const { return MESH_TRIPRISM15_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

//! A Pyramid
/*! A Pyramid cell

Node direction:


    t      3------------2  \n
    |     /            /   \n
    |    /            /    \n
    |   /            /     \n
    |  /      5     /      \n
    | /            /       \n
    |/            /        \n
    0------------1-------r \n

*/

static const uint8 PyramidFacesID[5][4] = {
    {1, 2, 5, 255},    // l
    {2, 3, 5, 255},    // l
    {0, 5, 3, 255},    // l
    {0, 1, 5, 255},    // l
    {0, 3, 2, 1},      // r
};

class DLLEXPORT Pyramid : public Cell {
public:
    Pyramid(const std::vector < Node * > & nodes);

    virtual ~Pyramid();

    virtual uint dim() const { return 3; }

    virtual uint rtti() const { return MESH_PYRAMID_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    virtual uint neighborCellCount() const { return 5; }

    /*! Experimental */
    virtual std::vector < Node * > boundaryNodes(Index i) const;

protected:
};

//*! VTK, Flaherty, Gimli count: 1-2-3-4, 5(1-2), 6(2-3), 7(3-1), 8(1-4), 9(2-4), 10(3-4)* //
static const uint8 Pyramid13NodeSplit[13][2] = {
    {0,0},{1,1},{2,2},{3,3},{4,4},
    {0,1},{1,2},{2,3},{3,0},
    {0,4},{1,4},{2,4},{3,4}
};

//! A Pyramid
/*! A Pyramid cell with 13 nodes for quadratic base functions

Node direction:


    t       3-----7------2  \n
    |      /            /    \n
    |     /  12    11  /     \n
    |    /            /      \n
    |   8     4      6       \n
    |  /            /        \n
    | /  9      10 /         \n
    |/            /          \n
    0------5-----1-------r   \n

*/

class DLLEXPORT Pyramid13 : public Pyramid {
public:
    Pyramid13(const std::vector < Node * > & nodes);

    virtual ~Pyramid13();

    virtual uint rtti() const { return MESH_PYRAMID13_RTTI; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

protected:
};

DLLEXPORT std::ostream & operator << (std::ostream & str, const std::set < GIMLI::MeshEntity * > & ents);

DLLEXPORT std::ostream & operator << (std::ostream & str, const std::set < GIMLI::Boundary * > & bounds);

} // namespace GIMLI

#endif // MESHENTITIES__H
