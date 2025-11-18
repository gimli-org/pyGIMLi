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

#ifndef _GIMLI_ELEMENTMATRIX__H
#define _GIMLI_ELEMENTMATRIX__H

#include "gimli.h"
#include "vector.h"
#include "matrix.h"

namespace GIMLI{

class FEAFunction;

template < class ValueType > class DLLEXPORT ElementMatrix {
public:
    /*! If dof != 0 then scalar field approximation is to be supposed.
    For vector field solution give a dof, means be the number of nodes of the current mesh. */
    ElementMatrix(Index dof=0){
        this->_nDof = dof;
        this->_newStyle = false;
        this->_nCoeff = 0;
        this->_dofPerCoeff = 0;
        this->_dofOffset = 0;
    }

    ~ElementMatrix() {}

    inline const Vector< ValueType > & operator[](Index row) const {
        return mat_[row]; }

    void resize(Index rows, Index cols=0) {
        if (cols == 0) cols = rows;
        _idsR.resize(rows);
        _idsC.resize(cols);
        _ids.resize(rows);
        mat_.resize(rows, cols);
    }

    ElementMatrix < ValueType > & operator += (const ElementMatrix < ValueType > & E){
        for (uint i = 0; i < size(); i ++){ mat_[i] += E.row(i); }
        return *this;
    }

    #define DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(OP)                   \
        ElementMatrix < ValueType > & operator OP##= (ValueType val) { \
            if (this->_newStyle){ \
                if (this->_integrated){ \
                    for (Index i = 0; i < size(); i ++) mat_[i] OP##= val; \
                } \
                for (auto & m: _matX){ \
                    m OP##= val; \
                } \
                return *this;\
            } \
            for (Index i = 0; i < size(); i ++) mat_[i] OP##= val; \
            return *this;\
        } \

        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(+)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(-)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(/)
        DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__(*)

    #undef DEFINE_ELEMENTMATRIX_UNARY_MOD_OPERATOR__

    inline Index size() const { return mat_.rows(); }
    inline Index rows() const { return mat_.rows(); }
    inline Index cols() const { return mat_.cols(); }

    inline const ValueType & getVal(Index i, Index j) const {
        return mat_[i][j]; }
    inline void setVal(Index i, Index j, const ValueType & v) {
        mat_[i][j] = v; }
    inline void addVal(Index i, Index j, const ValueType & v) {
        mat_[i][j] += v; }

    /*! Set data matrix. */
    inline void setMat(const Matrix < ValueType > & m) { mat_ = m; }
    /*! Return data matrix. */
    inline Matrix < ValueType > * pMat() { return & mat_; }
    /*! Return data for row i. */
    inline const Matrix < ValueType > & mat() const { return mat_; }
    /*! Return data for row i. */
    inline const Vector < ValueType > & row(Index i) const { return mat_[i]; }
    /*! Return data for col i. */
    inline Vector < ValueType > col(Index i) const {
        Vector < ValueType > ret(this->rows());
        for (Index j = 0; j < ret.size(); j ++){ ret[j] = mat_[j][i];}
        return ret;
    }

    /*! Fill the node ids with a number of coefficents.
    For vector field approximation give field dimension 2 or 3.
    Please note that you need to give the number of nodes to the ElementMatrix constructor.
    */
    void fillIds(const MeshEntity & ent, Index nC=1);

    /*! Set all node indices for row and columns. Can be unsymmetric.*/
    inline void setIds(const IndexArray & idsR, const IndexArray & idsC) {
        _idsR = idsR; _idsC = idsC; _ids = idsR;
    }

    /*! Set all node indices.*/
    inline void setIds(const IndexArray & ids) {
        _idsR = ids; _idsC = ids; _ids = ids;
    }
    /*! Return all node indices.*/
    inline const IndexArray & ids() const { return _ids; }
    /*! Return all row node indices.*/
    inline const IndexArray & rowIDs() const { return _idsR; }
    /*! Return all column node indices.*/
    inline const IndexArray & colIDs() const { return _idsC; }
    /*! Return the node index for node i.*/
    inline const Index idx(Index i) const {
        return _ids[i]; }

    /*! Fill this element matrix with int_boundary C * u */
    ElementMatrix < ValueType > & u(const MeshEntity & ent
                                    // const Matrix< ValueType > & C
                                    );

    /*! Fill this element matrix with int_domain C * u * u */
    ElementMatrix < ValueType > & u2(const MeshEntity & ent
                                    //  const Matrix< ValueType > & C
                                     );

    /*! Get integration weights and points for the entity. */
    void getWeightsAndPoints(const MeshEntity & ent,
                const RVector * &w, const PosVector * &x, int order);


    /*! Return the stress matrix for this entity.*/
    Vector < ValueType > stress(const MeshEntity & ent,
                                const Matrix< ValueType > & C,
                                const RVector & u, bool voigtNotation=false);

    /*! Return gradient base for the last entity, i.e., integrate over gradient u over the cell. */
    const Matrix < ValueType > & gradientBase() const { return this->_grad;}

    /*! Fill Element Gradients Matrix for all integration points. */
    void fillGradientBase(const MeshEntity & ent,
                          const RVector & w,
                          const PosVector & x,
                          Index nC,
                          bool voigtNotation);

    /*! Fill this element matrix with int_domain C * grad u*/
    ElementMatrix < ValueType > & gradU(const Cell & cell,
                                        Index nC,
                                        bool voigtNotation=false);

    /*! Fill this element matrix with int_domain C * grad u */
    ElementMatrix < ValueType > & gradU(const MeshEntity & ent,
                                        const RVector & w,
                                        const PosVector & x,
                                        Index nC,
                                        bool voigtNotation=false);


    /*! Fill this element matrix with int_domain C * grad u * grad u.
    For scalar field approximation define C.size() = (1x1) isotropic or anisotropic
    C.size() = (cell.dim() x cell.dim()) parameter. For vector field approximation
    create the ElementMatrix with appropriate dof and C can be of size = (1x1)
    for isotropic, (cell.dim() x cell.dim()) for anisotropic, or (3x3)
    for 2D elastic and (6x6) 3D elastic parameter.
    Notation for elastic parameters Kelvin notation as default and can be Voigt's notation if needed.
     */
    ElementMatrix < ValueType > & gradU2(const Cell & cell,
                                         const Matrix< ValueType > & C,
                                         bool voigtNotation=false);

    /*! Fill this element matrix with int_domain C * grad u * grad u.*/
    ElementMatrix < ValueType > & gradU2(const Cell & cell, ValueType c){
        Matrix < ValueType > C(1, 1);
        C[0][0] = c;
        return this->gradU2(cell, C, false);
    }

    /*! Fill this element matrix with int_domain C * grad u * grad u. */
    ElementMatrix < ValueType > & gradU2(const MeshEntity & ent,
                                         const Matrix< ValueType > & C,
                                         const RVector & w,
                                         const PosVector & x,
                                         bool voigtNotation=false);

    ElementMatrix < ValueType > & ux2uy2uz2(const Cell & cell,
                                            bool useCache=false);

    ElementMatrix < ValueType > & u(const MeshEntity & ent,
                                    const RVector & w,
                                    const PosVector & x,
                                    bool verbose=false);
    ElementMatrix < ValueType > & u2(const MeshEntity & ent,
                                     const RVector & w,
                                     const PosVector & x,
                                     bool verbose=false);
    ElementMatrix < ValueType > & ux2(const MeshEntity & ent,
                                      const RVector & w,
                                      const PosVector & x,
                                      bool verbose=false);
    ElementMatrix < ValueType > & ux2uy2(const MeshEntity & ent,
                                         const RVector & w,
                                         const PosVector & x,
                                         bool verbose=false);
    ElementMatrix < ValueType > & ux2uy2uz2(const MeshEntity & ent,
                                            const RVector & w,
                                            const PosVector & x,
                                            bool verbose=false);

    ElementMatrix < double > & dudi(const MeshEntity & ent,
                                  const RVector & w,
                                  const PosVector & x,
                                  Index i, bool verbose=false);

    ElementMatrix < double > & ux(const MeshEntity & ent,
                                  const RVector & w,
                                  const PosVector & x,
                                  bool verbose=false){
        return dudi(ent, w, x, 0, verbose);
    }
    ElementMatrix < double > & uy(const MeshEntity & ent,
                                  const RVector & w,
                                  const PosVector & x,
                                  bool verbose=false){
        return dudi(ent, w, x, 1, verbose);
    }
    ElementMatrix < double > & uz(const MeshEntity & ent,
                                  const RVector & w,
                                  const PosVector & x,
                                  bool verbose=false){
        return dudi(ent, w, x, 2, verbose);
    }

    Vector < ValueType > mult(const Vector < ValueType > & v){
        Vector < ValueType > ret(this->size());
        this->mult(v, ret);
        return ret;
    }

    /*! Return S * a */
    void mult(const Vector < ValueType > & a, Vector < ValueType > & ret){
        ASSERT_EQUAL(size(), ret.size())
        for (Index i = 0; i < size(); i ++) {
            for (Index j = 0; j < size(); j ++) {
                ret[i] += mat_[i][j] * a[_ids[j]];
            }
        }
    }
    /*! Return (S * a) * b */
    ValueType mult(const Vector < ValueType > & a,
                   const Vector < ValueType > & b){
        ValueType ret = 0;
        for (Index i = 0; i < size(); i ++) {
            ValueType t = 0;
            for (Index j = 0; j < size(); j ++) {
                t += mat_[i][j] * a[_ids[j]];
            }
            ret += t * b[_ids[i]];
        }
        return ret;
    }

    /*! Return (S * (a-b)) * (m-n) TODO Check if its the same like mult(a-b, m-n)*/
    template < class Val > Val mult_(const Vector < Val > & a,
                                     const Vector < Val > & b,
                                     const Vector < Val > & m,
                                     const Vector < Val > & n){
        Val ret = 0;
        for (Index i = 0; i < size(); i ++) {
            Val t = 0;
            for (Index j = 0; j < size(); j ++) {
                t += mat_[i][j] * (a[_ids[j]]-b[_ids[j]]);
            }
            ret += t * (m[_ids[i]]-n[_ids[i]]);
        }
        return ret;
    }

    double mult(const RVector & a, const RVector & b,
                const RVector & m, const RVector & n){
        return mult_(a, b, m, n);
    }
    Complex mult(const CVector & a, const CVector & b,
                const CVector & m, const CVector & n){
        return mult_(a, b, m, n);
    }

    //** new interface starts here **//
    ElementMatrix(Index nCoeff, Index dofPerCoeff, Index dofOffset);

    ElementMatrix(const ElementMatrix < ValueType > & E);

    void copyFrom(const ElementMatrix < ValueType > & E, bool withMat=true);

    void init(Index nCoeff, Index dofPerCoeff, Index dofOffset);

    /*! Fill this ElementMatrix with value (u for scalar, v for vector values) basis. Cache the matrix in entity*/
    ElementMatrix < ValueType > & pot(const MeshEntity & ent, Index order,
                                      bool sum,
                                      Index nCoeff, Index dof, Index dofOffset);

    /*! Fill this ElementMatrix with value (u for scalar, v for vector values) basis.*/
    ElementMatrix < ValueType > & pot(const MeshEntity & ent, Index order,
                                      bool sum=false);

    /*! Fill this ElementMatrix with value (u for scalar, v for vector values) basis. Cache the matrix in entity*/
    ElementMatrix < ValueType > & grad(const MeshEntity & ent, Index order,
                                       bool elastic, bool sum, bool div,
                                      Index nCoeff, Index dof, Index dofOffset, bool kelvin=false
                                      );

    /*! Fill this ElementMatrix with gradient of ent.*/
    ElementMatrix < ValueType > & grad(const MeshEntity & ent, Index order,
                                       bool elastic=false, bool sum=false,
                                       bool div=false, bool kelvin=false);

    /*! Integrate, i.e., sum over quadrature matrices.*/
    void integrate() const;

    /*! Return reference to all matrices per quadrature point.*/
    const std::vector < Matrix < ValueType > > & matX() const { return _matX; }

    std::vector < Matrix < ValueType > > * pMatX() { return &_matX; }

    /*! Return const reference to the last active entity.*/
    const MeshEntity & entity() const { ASSERT_PTR(_ent); return *_ent; }

    /*! Return const reference to the last quadrature points.*/
    const PosVector & x() const { ASSERT_PTR(_x); return *_x; }

    /*! Return const reference to the last quadrature weights.*/
    const RVector & w() const { ASSERT_PTR(_w); return *_w; }

    Index order() const { return _order; }
    Index nCoeff() const { return _nCoeff; }
    Index dofPerCoeff() const { return _dofPerCoeff; }
    Index dofOffset() const { return _dofOffset; }

    void setDiv(bool div){ _div = true;}
    bool isDiv() const { return _div;}

    bool isIntegrated() const { return _integrated; }
    void integrated(bool i) { _integrated = i; }

    bool valid() const { return _valid; }
    void setValid(bool v) { _valid = v; }

    bool elastic() const { return _elastic;}

    bool oldStyle() const { return !this->_newStyle; }
protected:
    mutable Matrix < ValueType > mat_;
    IndexArray _ids;
    IndexArray _idsC;
    IndexArray _idsR;

    std::map< uint, RVector > uCache_;
    std::map< uint, Matrix < ValueType > > u2Cache_;

    std::vector< Matrix < ValueType > > _B;
    Matrix < ValueType > _grad;

    // number of single dof
    Index _nDof;
    // number of coefficients: 1, 2, 3 for scalar(dim), 1, 3, 6 for vector(dim)
    // Index _nC;

    RMatrix dNdr_;
    RMatrix dNds_;
    RMatrix dNdt_;

    RMatrix dNdx_; // (nRules, nVerts)
    RMatrix dNdy_; // (nRules, nVerts)
    RMatrix dNdz_; // (nRules, nVerts)

    RMatrix _abaTmp; // temp workspace

    //** new interface starts here **//
    // const Mesh * _mesh;
    Index _order;
    Index _nCoeff;
    Index _dofPerCoeff;
    Index _dofOffset;

    const MeshEntity * _ent;
    const RVector * _w;
    const PosVector * _x;
    // matrices per quadrature point
    std::vector < Matrix < ValueType > > _matX;

    bool _newStyle;
    bool _div;
    bool _valid;
    bool _elastic;
    mutable bool _integrated;

private:
    /*! No copy operator. */

    /*! No assignment operator. */
    ElementMatrix < ValueType > & operator = (const ElementMatrix < ValueType > & E) {
        std::cout << "ElementMatrix::operator = (" << std::endl;
        THROW_TO_IMPL
        if (this != & E){
//             this->resize(E.size());
//             for (uint i = 0; i < E.size(); i ++) mat_[i] = E.row(i);
//             _ids = E.idx();
        } return *this;
    }
};

template < > DLLEXPORT
void ElementMatrix < double >::copyFrom(const ElementMatrix < double > & E,
                                        bool withMat);

template < > DLLEXPORT
void ElementMatrix < double >::init(Index nCoeff, Index dofPerCoeff,
                                    Index dofOffset);

DLLEXPORT void dot(const ElementMatrix < double > & A,
                   const ElementMatrix < double > & B,
                   double c, ElementMatrix < double > & ret);

DLLEXPORT void dot(const ElementMatrix < double > & A,
                   const ElementMatrix < double > & B,
                   const Pos & c, ElementMatrix < double > & ret);

DLLEXPORT void dot(const ElementMatrix < double > & A,
                   const ElementMatrix < double > & B,
                   const RMatrix & c, ElementMatrix < double > & ret);

DLLEXPORT void dot(const ElementMatrix < double > & A,
                   const ElementMatrix < double > & B,
                   const FEAFunction & c, ElementMatrix < double > & ret);

// DLLEXPORT void dot(const ElementMatrix < double > & A,
//                    const ElementMatrix < double > & B,
//                    A_TYPE c, ElementMatrix < double > & C);

#define DEFINE_DOT_MULT(A_TYPE) \
DLLEXPORT const ElementMatrix < double > dot( \
                                        const ElementMatrix < double > & A, \
                                        const ElementMatrix < double > & B, \
                                        A_TYPE c); \
DLLEXPORT void mult(const ElementMatrix < double > & A, A_TYPE b, \
                    ElementMatrix < double > & C); \
DLLEXPORT const ElementMatrix < double > mult( \
                    const ElementMatrix < double > & A, A_TYPE b); \

DEFINE_DOT_MULT(double)
DEFINE_DOT_MULT(const Pos &)
DEFINE_DOT_MULT(const RMatrix &)
DEFINE_DOT_MULT(const FEAFunction &)

#undef DEFINE_DOT_MULT

DLLEXPORT const ElementMatrix < double > dot(
                                        const ElementMatrix < double > & A,
                                        const ElementMatrix < double > & B);
// return dot(A, B, 1.0);}

DLLEXPORT void dot(const ElementMatrix < double > & A,
                   const ElementMatrix < double > & B,
                   ElementMatrix < double > & ret);
// return dot(A, B, 1.0, ret);
// }


/*! scalar per quadrature point */
DLLEXPORT void mult(const ElementMatrix < double > & A, const RVector & b,
                    ElementMatrix < double > & C);
/*! vector per quadrature point */
DLLEXPORT void mult(const ElementMatrix < double > & A, const PosVector & b,
                    ElementMatrix < double > & C);
/*! Matrix per quadrature point */
DLLEXPORT void mult(const ElementMatrix < double > & A,
                    const std::vector < RMatrix > & b,
                    ElementMatrix < double > & C);

/*!Evaluate scalars per cell.*/
DLLEXPORT void evaluateQuadraturePoints(const MeshEntity & ent,
                                        const PosVector & x,
                                        const FEAFunction & f, RVector & ret);
/*!Evaluate vectors per cell.*/
DLLEXPORT void evaluateQuadraturePoints(const MeshEntity & ent,
                                        const PosVector & x,
                                        const FEAFunction & f, PosVector & ret);
/*!Evaluate matrices per cell.*/
DLLEXPORT void evaluateQuadraturePoints(const MeshEntity & ent,
                                        const PosVector & x,
                                        const FEAFunction & f,
                                        std::vector < RMatrix > & ret);
/*!Evaluate scalar for each cell.*/
DLLEXPORT void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                                        const FEAFunction & f,
                                        std::vector< RVector > & ret);
/*!Evaluate vectors for each cell.*/
DLLEXPORT void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                                        const FEAFunction & f,
                                        std::vector< PosVector > & ret);
/*!Evaluate matrices for each cell.*/
DLLEXPORT void evaluateQuadraturePoints(const Mesh & mesh, Index order,
                                        const FEAFunction & f,
                              std::vector< std::vector< RMatrix > > & ret);


#define DEFINE_CREATE_FORCE_VECTOR(A_TYPE) \
DLLEXPORT void createForceVector(const Mesh & mesh, Index order, \
                                 RVector & ret, A_TYPE a, \
                                 Index nCoeff, Index dofOffset); \
DLLEXPORT void createMassMatrix(const Mesh & mesh, Index order, \
                                RSparseMapMatrix & ret, A_TYPE a, \
                                Index nCoeff, Index dofOffset); \
DLLEXPORT void createStiffnessMatrix(const Mesh & mesh, Index order, \
                                     RSparseMapMatrix & ret, A_TYPE a, \
                                     Index nCoeff, Index dofOffset, \
                                     bool elastic=false, bool kelvin=false); \

DEFINE_CREATE_FORCE_VECTOR(double)             // const scalar for all cells
DEFINE_CREATE_FORCE_VECTOR(const Pos &)   // const vector for all cells
DEFINE_CREATE_FORCE_VECTOR(const RVector &)    // scalar for each cell
DEFINE_CREATE_FORCE_VECTOR(const PosVector &)   // vector for each cell
DEFINE_CREATE_FORCE_VECTOR(const RMatrix &)    // const matrix for all cells
// matrix for each cell
DEFINE_CREATE_FORCE_VECTOR(const std::vector< RVector > &)
// vector at each quadrature point for each cell
DEFINE_CREATE_FORCE_VECTOR(const std::vector< PosVector > &)
// matrix at each quadrature point for each cell
DEFINE_CREATE_FORCE_VECTOR(const std::vector < RMatrix > &)
// scalar at each quadrature point for each cell
DEFINE_CREATE_FORCE_VECTOR(const std::vector< std::vector < RMatrix > > &)
// generic function for each point
DEFINE_CREATE_FORCE_VECTOR(const FEAFunction &)

#undef DEFINE_CREATE_FORCE_VECTOR

/*!Interface to function q=f(p, ent) with q, p = Pos() in R1, R2, R 3
and ent assiated mesh entity.*/
class DLLEXPORT FEAFunction {
public:
    FEAFunction(Index valueSize): _valueSize(valueSize){ }

    virtual ~FEAFunction() { }

    virtual double evalR1(const Pos & arg, const MeshEntity * ent=0) const{
        log(Warning, "FEAFunction.eval should be overloaded.");
        return 0.0;
    }
    virtual Pos evalR3(const Pos & arg, const MeshEntity * ent=0) const{
        log(Warning, "FEAFunction.eval should be overloaded.");
        return Pos(0.0, 0.0, 0.0);
    }
    virtual RMatrix evalRM(const Pos & arg, const MeshEntity * ent=0) const{
        log(Warning, "FEAFunction.eval should be overloaded.");
        return RMatrix(0, 0);
    }

    Index valueSize() const { return _valueSize; }

protected:
    Index _valueSize;
};

class DLLEXPORT ElementMatrixMap {
public:

    void add(Index row, const ElementMatrix < double > & Ai);

    //TODO .. check if its the same like mult(a-b, m-n))
    RVector mult(const RVector & a, const RVector & b,
                 const RVector & m, const RVector & n) const;

    /*! Return (S_i * a) * b for all i*/
    RVector mult(const RVector & a, const RVector & b) const;

    Index rows() const { return rows_; }

    Index cols() const { return cols_; }


protected:
    std::vector< RMatrix > mat_;
    std::vector< IndexArray > _ids;
    std::vector< Index > row_;

    Index rows_;
    Index cols_;
};

template < class ValueType > std::ostream & operator << (std::ostream & str,
                                                         const ElementMatrix< ValueType > & e);

template < > DLLEXPORT std::ostream & operator << (std::ostream & str,
                                            const ElementMatrix< double > & e);



} // namespace GIMLI{

#endif // _GIMLI_ELEMENTMATRIX__H
