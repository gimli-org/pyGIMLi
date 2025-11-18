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

#ifndef _GIMLI_Shape__H
#define _GIMLI_Shape__H

#include "gimli.h"
#include "polynomial.h"
#include "curvefitting.h"

#ifndef PYGIMLI_CAST // fails because of boost threads and clang problems
    #if USE_BOOST_THREAD
        #include <boost/thread.hpp>
        extern boost::mutex ShapeFunctionWriteCacheMutex__;
        //#error AVOID BOOST
    #else
        #include <mutex>
        extern std::mutex ShapeFunctionWriteCacheMutex__;
    #endif
#endif

namespace GIMLI{

/*! Return Volume of the Tetrahedron given by 4 RVector3 */
DLLEXPORT double tetVolume(const RVector3 & p0, const RVector3 & p1,
                           const RVector3 & p2, const RVector3 & p3);

/*! Return Size of the Triangle given by 3 RVector3 */
DLLEXPORT double triSize(const RVector3 & p0, const RVector3 & p1,
                         const RVector3 & p2);


DLLEXPORT std::vector < PolynomialFunction < double > >
createPolynomialShapeFunctions(const std::vector < RVector3 > & pnts,
                               uint dim, uint nCoeff,
                               bool pascale, bool serendipity,
                               const RVector & startVector);


template < class Ent > std::vector < PolynomialFunction < double > >
    createPolynomialShapeFunctions(const Ent & ent, uint nCoeff,
                                   bool pascale, bool serendipity,
                                   const RVector & startVector){
// __MS(ent)
    std::vector < RVector3 > pnts;
    for (Index i = 0; i < ent.nodeCount(); i ++){
        pnts.push_back(ent.rst(i));
    }

    return createPolynomialShapeFunctions(pnts, ent.dim(), nCoeff, pascale,
                                          serendipity, startVector);
}

// no default arg here .. pygimli@win64 linker bug
template < class Ent > std::vector < PolynomialFunction < double > >
    createPolynomialShapeFunctions(const Ent & ent, uint nCoeff,
                                   bool pascale, bool serendipity){
    RVector start(0);
    return createPolynomialShapeFunctions(ent, nCoeff, pascale, serendipity, start);
}

class DLLEXPORT ShapeFunctionCache : public Singleton< ShapeFunctionCache > {
public:
    friend class Singleton< ShapeFunctionCache >;

    template < class Ent > const std::vector < PolynomialFunction < double > > &
    shapeFunctions(const Ent & e) const {

        std::map < uint8, std::vector < PolynomialFunction < double > > >::const_iterator it = shapeFunctions_.find(e.rtti());

        if (it == shapeFunctions_.end()){
            this->createShapeFunctions_(e);
            it = shapeFunctions_.find(e.rtti());
        }

        return (*it).second;
    }

    template < class Ent > const std::vector < PolynomialFunction < double > > &
    deriveShapeFunctions(const Ent & e , uint dim) const {
        std::map < uint8, std::vector < std::vector < PolynomialFunction < double > > > >::const_iterator it = dShapeFunctions_.find(e.rtti());

        if (it == dShapeFunctions_.end()) {
            this->createShapeFunctions_(e);
            it = dShapeFunctions_.find(e.rtti());
        }
        return (*it).second[dim];
    }

    /*! Clear the cache. */
    void clear() {
        shapeFunctions_.clear();
        dShapeFunctions_.clear();
    }

    std::vector< RMatrix3 > & RMatrix3Cache();
    std::vector< RMatrix > & RMatrixCache(uint rtti);

    RMatrix3 & cachedRMatrix3(uint i);
    RMatrix & cachedRMatrix(uint rtti, uint i);
private:

    /*! probably threading problems .. pls check*/
    template < class Ent > void createShapeFunctions_(const Ent & e) const {
        std::vector < PolynomialFunction < double > > N = e.createShapeFunctions();

        #if USE_BOOST_THREAD
        //#ifdef WIN32_LEAN_AND_MEAN
        //        __MS("pls check missing mutex")
            //boost::mutex::scoped_lock lock(ShapeFunctionWriteCacheMutex__);
        //#else
            #ifndef PYGIMLI_CAST // fails because of boost threads and clang problems
            boost::mutex::scoped_lock lock(ShapeFunctionWriteCacheMutex__);
            #endif
        //#endif

        #else
        // __MS("lock deactiated " << ShapeFunctionWriteCacheMutex__)

            std::unique_lock < std::mutex > lock(ShapeFunctionWriteCacheMutex__);
        #endif


        shapeFunctions_[e.rtti()] = N;
        dShapeFunctions_[e.rtti()] = std::vector < std::vector < PolynomialFunction < double > > >();

        dShapeFunctions_[e.rtti()].push_back(std::vector < PolynomialFunction < double > >());
        dShapeFunctions_[e.rtti()].push_back(std::vector < PolynomialFunction < double > >());
        dShapeFunctions_[e.rtti()].push_back(std::vector < PolynomialFunction < double > >());

        for (uint i = 0; i < N.size(); i ++){
            dShapeFunctions_[e.rtti()][0].push_back(N[i].derive(0));
            dShapeFunctions_[e.rtti()][1].push_back(N[i].derive(1));
            dShapeFunctions_[e.rtti()][2].push_back(N[i].derive(2));
        }
    }

    /*! Private so that it can not be called */
    ShapeFunctionCache(){}
    /*! Private so that it can not be called */
    virtual ~ShapeFunctionCache(){}
    /*! Copy constructor is private, so don't use it */
    ShapeFunctionCache(const ShapeFunctionCache &){};
    /*! Assignment operator is private, so don't use it */
    void operator = (const ShapeFunctionCache &){};

protected:

    /*! Cache for shape functions. */
    mutable std::map < uint8, std::vector < PolynomialFunction < double > > > shapeFunctions_;

    /*! Cache for shape functions derivatives. */
    mutable std::map < uint8, std::vector< std::vector < PolynomialFunction < double > > > > dShapeFunctions_;

    mutable std::vector< RMatrix3 > _rMatrix3Cache;
    mutable std::map< uint, std::vector< RMatrix > > _rMatrixCache;

};

static const double NodeCoordinates[1][3] = {
    {0.0, 0.0, 0.0}
};

//! A Shape defines a geometrical primitive.
/*! A Shape defines a geometrical primitive. And is the geometrical base of a \ref MeshEntity. */
class DLLEXPORT Shape {
public:
    /*! Default constructor. */
    Shape(MeshEntity * ent);

    /*! Default destructor. */
    virtual ~Shape();

    /*! Pure virtual methode for runtime identification. */
    virtual int rtti() const = 0;

    virtual int dim() const = 0;

    /*! Return an identification name for the shape. */
    virtual std::string name() const { return "Shape"; }

    /*! Return the amount of nodes for this shape. */
    Index nodeCount() const { return this->nodeCount_; }

    /*! Return a read only reference to the i-th \ref Node of this shape. */
    Node & node(Index i) const;

    // /*! Return a reference to the i-th \ref Node of this shape. */
    // Node & node(Index i);

    /*! Set the ptr to nodes array from, which the MeshEntity holds. */
    void setNodesPtr(const std::vector< Node * > & n) {
        this->nodeVector_ = &n;
    }

    /*! Return a read only reference to all nodes. */
    const std::vector< Node * > & nodes() const { return * this->nodeVector_; }

    /*! Return a reference to all nodes. */
    // std::vector < Node * > & nodes() { return nodeVector_ ; }

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const ;

    /*! Create and return the Jacobian matrix for this shape at the local coordinate \f$ L[r,s,t]\f$. \n
     * The Jacobian matrix is a 3x3 matrix representing the absolute derivative matrix for this shape regarding the \f$ \mathcal{N} \f$ shape functions:\n
     \f{eqnarray*}{
        J(0,i-1) & = \sum_n^{\mathcal{N}}\frac{dN_n(r,s,t)}{\partial r,s,t}x_n \\
        J(1,i-1) & = \sum_n^{\mathcal{N}}\frac{dN_n(r,s,t)}{\partial r,s,t}y_n \\
        J(2,i-1) & = \sum_n^{\mathcal{N}}\frac{dN_n(r,s,t)}{\partial r,s,t}z_n
     \f} for i = [1,2,3]

     \f{eqnarray*}{
                    & \frac{\partial x}{\partial r}, \frac{\partial x}{\partial s}, \frac{\partial x}{\partial t}\\
        J      =    & \frac{\partial y}{\partial r}, \frac{\partial y}{\partial s}, \frac{\partial y}{\partial t}\\
                    & \frac{\partial z}{\partial r}, \frac{\partial z}{\partial s}, \frac{\partial z}{\partial t}\\
     \f}
     the inverse of the Jacobian results in:
     \f{eqnarray*}{
                               & \frac{\partial r}{\partial x}, \frac{\partial r}{\partial y}, \frac{\partial r}{\partial z}\\
        J^{-1} = \frac{1}{|J|} & \frac{\partial s}{\partial x}, \frac{\partial s}{\partial y}, \frac{\partial s}{\partial z}\\
                               & \frac{\partial t}{\partial x}, \frac{\partial t}{\partial y}, \frac{\partial t}{\partial z}
        \f}
    */
    void createJacobian(RMatrix3 & J) const;

    RMatrix3 createJacobian() const;

    /*! Return the inverse of the Jacobian Matrix. And create and cache it on demand.
     * The matrix is no more valid if the shape was transformed.
     */
    const RMatrix3 & invJacobian() const;

    /*! Return a \ref RVector for the \f$ n=[0,\mathrm{nodeCount()}] \f$ shape functions
     *\f$ N_n(r,s,t)\f$ at the local coordinate \f$ L(r,s,t)\f$. */
    virtual RVector N(const RVector3 & L) const;

    /*! TODO replace this with expressions. \n
     * Fill the allocated \ref RVector n of size \f$ \mathcal{N} \f$ for the \f$ n=[0,\f$ \mathcal{N} \f$) \f$ shape functions \f$ N_n(r,s,t)\f$ for the local coordinate \f$ L(r,s,t)\f$. For performance reasons the size of L will not be checked.
     This method needs to be specialized in the corresponding child classes. */
    virtual void N(const RVector3 & L, RVector & ret) const;

    /*! Return the derivative matrix of size (\f$ 3\times\mathcal{N} \f$) for the shape functions at the local coordinate.
     * \f$ L(r,s,t) \f$. Result is independent of L for linear shape function (TODO Remove on cleanup)
     * \f$ [[\frac{dN_i(r,s,t)}{\partial r}],[\frac{dN_i(r,s,t)}{\partial s}],[\frac{dN_i(r,s,t)}{\partial t}]^{\mathrm{T}}] \f$ for \f$ i = [0,\mathcal{N}\f$ */
    virtual void dNdrst(const RVector3 & rst, RMatrix & MdNdrst) const;

    virtual RMatrix dNdrst(const RVector3 & L) const;

    /*! Perform coordinate transformation from the locale coordinates \f$ (r,s,t)=(r,s,t)=([0..1,0..1,0..1]) \f$ of this shape to Cartesian coordinates \f$ (x,y,z) \f$ regarding to the \f$ \mathcal{N} \f$ shape functions \ref N
     * \f$ N_i \f$ with \f$ i=[0,\mathcal{N})\f$ \n
     * This is the opposite to \ref xyz2rst().
     * \f{eqnarray*}{
     *  x &=& \sum_i^{\mathcal{N}} N_i(r,s,t)x_i \\
     *  y &=& \sum_i^{\mathcal{N}} N_i(r,s,t)y_i \\
     *  z &=& \sum_i^{\mathcal{N}} N_i(r,s,t)z_i
        \f}
     * This is a generic function and may be overwritten by faster methods for shape simplexes.
     */
    virtual void rst2xyz(const RVector3 & rst, RVector3 & xyz) const;

    /*! Return the Cartesian coordinates for the locale coordinates rst. See \ref rst2xyz. */
    virtual RVector3 xyz(const RVector3 & rst) const;

    /*! Convert Cartesian coordinates into locale coordinates regarding the shape functions.
     * This is the opposite to \ref xyz2rst.
     * This is a generic function and may be overwritten by faster methods for shape simplexes.
     * Solve the nonlinear system of equations by newton method \f$ \mathrm{d}(x,y,z)= J*\mathrm{d}(r,s,t) \f$
     */
    virtual void xyz2rst(const RVector3 & xyz, RVector3 & rst) const;

    /*! Return local coordinates for Cartesian coordinates regarding the shape function. */
    virtual RVector3 rst(const RVector3 & xyz) const;

    /*! Return local coordinates for node i. */
    virtual RVector3 rst(Index i) const;

    /*! Return derivative from local coordinates to Cartesian coordinates.
     * These are the elements of the inverse Jacobian matrix.
     */
    inline double drstdxyz(uint rstI, uint xyzJ) const {
//         return invJacobian()[rstI][xyzJ];}
        return invJacobian()[rstI * 3 + xyzJ];}

    /*! Return true if the Cartesian coordinates xyz are inside the shape.
     * On boundary means inside too.
     Works only for shapes dedicated as cells because they need to
     be aligned to the dimension. See also \ref touch. */
    virtual bool isInside(const RVector3 & xyz, bool verbose=false) const;

    /*! Return true if the Cartesian coordinates xyz are inside the shape.
     * On boundary means inside too.
     * sf contains the complete shape function to identify next neighbor.
     Works only for shapes dedicated as cells because they need to
     be aligned to the dimension. See also \ref touch.*/
    virtual bool isInside(const RVector3 & xyz, RVector & sf,
                          bool verbose=false) const;

    /*! Check if the position touches the entity.
    Works only for shapes dedicated as boundaries.
    On edge of the boundary means inside too. See also \ref isInside.*/
    virtual bool touch(const RVector3 & pos, double tol=1e-6, bool verbose=false) const;

    /*!* Return true if the ray intersects the shape.
     * On boundary means inside too. The intersection position is stored in pos.
     * */
    virtual bool intersectRay(const RVector3 & start, const RVector3 & dir,
                              RVector3 & pos){
        __MS(rtti() << " " << name())
        THROW_TO_IMPL
        return false;
    };

    /*! Get the domain size of this shape, i.e., length, area or volume */
    double domainSize() const;

    /*! Returns the middle position of this shape */
    RVector3 center() const;

    /*! Returns the norm vector if possible otherwise returns non valid Vector3 */
    virtual RVector3 norm() const;

    /*! Returns maximum distance between 2 nodes.*/
    double h() const;

    /*! Returns the a plane for this shape if its possible (2D or 3D plane shapes) otherwise returns non valid Plane. */
    virtual Plane plane() const;

    /*! Notify this shape that the inverse Jacobian matrix and the domain size are not longer valid and need recalculation. This method is called if a node has bee transformed. */
    void changed();

    double jacobianDeterminant() const { return det(this->createJacobian()); }

    inline void resizeNodeSize_(Index n) { this->nodeCount_ = n; }


protected:

    /*! Virtual method to calculate the domain size i.e length, area, volume of the shapes */
    virtual double domainSize_() const { return 0.0; }

    Index nodeCount_;

    mutable double domSize_;
    mutable bool hasDomSize_;
    mutable double _h;

    mutable RMatrix3 invJacobian_;

    // const std::vector< Node * > & nodes()
    const std::vector < Node * > * nodeVector_;
};

DLLEXPORT std::ostream & operator << (std::ostream & str, const Shape & c);

class DLLEXPORT NodeShape : public Shape{
public:
    NodeShape(MeshEntity * ent) : Shape(ent) { resizeNodeSize_(1); }

    virtual ~NodeShape(){ }

    virtual int rtti() const { return MESH_SHAPE_NODE_RTTI; }

    virtual int dim() const { return 0; }

    virtual std::string name() const { return "NodeShape"; }

    virtual RVector3 rst(Index i) const;

    virtual double domainSize_() const { return 1.0; }

    virtual RVector3 norm() const;

protected:
    //virtual double jacobianDeterminant_() const { return 0.0; }
};

static const double EdgeCoordinates[2][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0}
};

class DLLEXPORT EdgeShape : public Shape {
public:
    EdgeShape(MeshEntity * ent) : Shape(ent) { resizeNodeSize_(2); }

    virtual ~EdgeShape(){ }

    virtual int rtti() const { return MESH_SHAPE_EDGE_RTTI; }

    virtual int dim() const { return 1; }

    virtual std::string name() const { return "EdgeShape"; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

    /*!* Return true if the ray intersects the shape.
     * On boundary means inside too. The intersection position is stored in pos.
     * */
    virtual bool intersectRay(const RVector3 & start, const RVector3 & dir,
                              RVector3 & pos);

    /*! Check if the position touches the entity.
    Works only for shapes dedicated as boundaries.
    On edge of the boundary means inside too. See also \ref isInside.*/
    virtual bool touch(const RVector3 & pos, double tol=1e-6, bool verbose=false) const;

    double length() const;

    virtual RVector3 norm() const;

protected:

    virtual double domainSize_() const { return length(); }
};

static const double TriCoordinates[3][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0}
};

/*! Triangle shape.

Cartesian \f$ (x,y) \f$ to natural coordinates \f$ (r,s,t) \f$ relation:\n
\f{eqnarray*}{
        x-x_1 &=& x_{21}s + x_{31} t\\
        y-y_1 &=& y_{21}s + y_{31} t\\
        r   &=& 1-s-t
   \f}
whereas \f$ x_{ij} \f$ reads \f$ x_i - x_j \f$ for the \f$x\f$-coordinate of \ref Node \f$i\f$ and \f$j\f$ respectively.
*/
class DLLEXPORT TriangleShape : public Shape {
public:

    TriangleShape(MeshEntity * ent) : Shape(ent) { resizeNodeSize_(3); }

    virtual ~TriangleShape(){ }

    virtual int rtti() const { return MESH_SHAPE_TRIANGLE_RTTI; }

    virtual int dim() const { return 2; }

    virtual std::string name() const { return "TriangleShape"; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

    /*! See Shape::xyz2rst. this is a specialized override for speedup. */
    virtual void xyz2rst(const RVector3 & pos, RVector3 & rst) const;

    double area() const;

    virtual RVector3 norm() const;

    /*!* Return true if the ray intersects the shape.
     * On boundary means inside too. The intersection position is stored in pos.
     * */
    virtual bool intersectRay(const RVector3 & start, const RVector3 & dir,
                              RVector3 & pos);

protected:

    /*! Interface to get the size i.e. area of this \ref TriangleShape */
    virtual double domainSize_() const { return area(); }
};


static const double QuadCoordinates[4][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 0.0}
};

//! Quadrangle shape.
/*!
 * Quadrangle shape. Convex quadrangle. Every inner angle might be arbitrary but lower than 180 degrees.
 *  3-------2
 *  |       |
 *  s       |
 *  |       |
 *  0---r---1
 *
 * r from 0 to 1
 * s from 0 to 1
*/
class DLLEXPORT QuadrangleShape : public Shape {
public:

    QuadrangleShape(MeshEntity * ent) : Shape(ent) { resizeNodeSize_(4); }

    virtual ~QuadrangleShape(){ }

    virtual int rtti() const { return MESH_SHAPE_QUADRANGLE_RTTI; }

    virtual int dim() const { return 2; }

    virtual std::string name() const { return "QuadrangleShape"; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

    double area() const;

//     virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

//     /*! See Shape::N. */
//     virtual void N(const RVector3 & L, RVector & n) const;

//     /*! See Shape::dNdrst. */
//     virtual RMatrix dNdrst(const RVector3 & rst) const;
//
    virtual RVector3 norm() const;

    virtual bool intersectRay(const RVector3 & start, const RVector3 & dir,
                              RVector3 & pos);

protected:
    virtual double domainSize_() const { return area(); }
};


class DLLEXPORT PolygonShape : public Shape {
public:
    PolygonShape(MeshEntity * ent);

    virtual ~PolygonShape();

    virtual int rtti() const { return MESH_SHAPE_POLYGON_FACE_RTTI; }

    virtual int dim() const { return 3; }

    virtual std::string name() const { return "PolygonFaceShape"; }

    virtual RVector3 norm() const;

    virtual RVector3 rst(Index i) const;

    virtual bool isInside(const RVector3 & xyz, bool verbose) const;

protected:
};


static const double TetCoordinates[4][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
};

//! Tetrahedral shape.
/*!
 * Tetrahedral shape.
*/
class DLLEXPORT TetrahedronShape : public Shape {
public:
    TetrahedronShape(MeshEntity * ent) : Shape(ent) { resizeNodeSize_(4); }

    virtual ~TetrahedronShape(){ }

    virtual int rtti() const { return MESH_SHAPE_TETRAHEDRON_RTTI; }

    virtual int dim() const { return 3; }

    virtual std::string name() const { return "TetrahedronShape"; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

    /*! See Shape::xyz2rst. Specialization for speedup */
    void xyz2rst(const RVector3 & pos, RVector3 & rst) const;

//     /*! See Shape::N. */
//     virtual void N(const RVector3 & L, RVector & n) const;
//
//     /*! See Shape::dNdrst. */
//     virtual RMatrix dNdrst(const RVector3 & rst) const;

    double volume() const;

protected:

    virtual double domainSize_() const { return volume(); }
};

static const int HexahedronSplit5TetID[5][4] = {
    {1, 4, 5, 6},
    {3, 6, 7, 4},
    {1, 0, 4, 3},
    {1, 2, 3, 6},
    {1, 4, 6, 3}
};

static const int HexahedronSplit6TetID[6][4] = {
    {0, 1, 2, 6},
    {0, 2, 3, 6},
    {0, 1, 6, 5},
    {0, 4, 5, 6},
    {0, 3, 7, 6},
    {0, 4, 6, 7}
};

static const double HexCoordinates[8][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 0.0, 1.0},
    {1.0, 1.0, 1.0},
    {0.0, 1.0, 1.0}
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

*/

class DLLEXPORT HexahedronShape : public Shape {
public:
    HexahedronShape(MeshEntity * ent) : Shape(ent){ resizeNodeSize_(8); }

    virtual ~HexahedronShape(){ }

    virtual int rtti() const { return MESH_SHAPE_HEXAHEDRON_RTTI; }

    virtual int dim() const { return 3; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

    virtual std::string name() const { return "HexahedronShape"; }

    double volume() const;

    /*! Special version of since simple order reverse isn't enough here.
     * We try to simple swap up and down side. */
    virtual bool enforcePositiveDirection();

//     /*! See Shape::N. */
//     virtual void N(const RVector3 & L, RVector & n) const;
//
//     /*! See Shape::dNdrst. */
//     virtual RMatrix dNdrst(const RVector3 & rst) const;

protected:

    virtual double domainSize_() const { return volume(); }
};

static const uint8 TriPrimSplit3TetID[3][4] = {
    {0, 5, 1, 2},
    {3, 1, 5, 4},
    {3, 5, 1, 0}
};

static const double PrismCoordinates[6][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 0.0, 1.0},
    {0.0, 1.0, 1.0}
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

class DLLEXPORT TriPrismShape : public Shape {
public:
    TriPrismShape(MeshEntity * ent) : Shape(ent){ resizeNodeSize_(6); }

    virtual ~TriPrismShape(){ }

    virtual int rtti() const { return MESH_SHAPE_TRIPRISM_RTTI; }

    virtual int dim() const { return 3; }

    virtual std::string name() const { return "TriagonalPrismShape"; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

    virtual std::vector < PolynomialFunction < double > > createShapeFunctions() const;

    double volume() const;

protected:

    virtual double domainSize_() const { return volume(); }
};

static const double PyramidCoordinates[5][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 0.1}
};

//! Pyramid
/*!
 *         4    \n
 *              \n
 *     3------2 \n
 *    /      /  \n
 *   0------1   \n
 */
class DLLEXPORT PyramidShape : public Shape {
public:
    PyramidShape(MeshEntity * ent) : Shape(ent){ resizeNodeSize_(5); }

    virtual ~PyramidShape(){ }

    virtual int rtti() const { return MESH_SHAPE_PYRAMID_RTTI; }

    virtual int dim() const { return 3; }

    /*! See Shape::rst */
    virtual RVector3 rst(Index i) const;

//     /*! See Shape::N. */
//     virtual void N(const RVector3 & L, RVector & n) const;
//
//     /*! See Shape::dNdrst. */
//     virtual RMatrix dNdrst(const RVector3 & rst) const;

};

} // namespace GIMLI

#endif // _GIMLI_SHAPE__H

