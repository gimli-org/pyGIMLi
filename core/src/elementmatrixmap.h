/******************************************************************************
 *   Copyright (C) 2006-2025 by the GIMLi development team                    *
 *   Carsten Rücker carsten@gimli.org
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

#pragma once

#include "gimli.h"
#include "elementmatrix.h"
#include "vector.h"
#include "matrix.h"

namespace GIMLI{

class DLLEXPORT ElementMatrixMap {
public:
    ElementMatrixMap(){
        rows_ = 0; // think to remove
        cols_ = 0; // think to remove
        dofA_ = 0;
        dofB_ = 0;
        nCoeff_ = 0;
        dofPerCoeff_ = 0;
        dofOffset_ = 0;
    }

    void push_back(const ElementMatrix < double > & Ai);

    void resize(Index size);

    void clear();

    #define DEFINE_INTEGRATOR_LF(A_TYPE, varname) \
        /*! R = \int_mesh this * f \d d mesh and R = RVector(dof) and \
            f = A_TYPE \
        */ \
        void integrate(const A_TYPE & varname, RVector & R, const double & alpha=1.0) const; \
        /*! Integrate into linear form R = alpha * \int_mesh this * f * R \d d mesh and \
        R = RVector(dof) and f = A_TYPE \
        */ \
        void integrate(const A_TYPE & varname, RVector & R, const RVector & alpha) const; \
        /*! Integrate into linear form R = alpha[cellID]*\int_mesh this * f * R \d d mesh and \
        R = RVector(dof) and f = A_TYPE \
        */ \
        RVector integrate(const A_TYPE & varname, const double & alpha=1.0) const; \
        /*! R = this * f */ \
        void mult(const A_TYPE & varname, ElementMatrixMap & ret) const;
    DEFINE_INTEGRATOR_LF(double, d)   // const scalar for all cells
    DEFINE_INTEGRATOR_LF(RSmallMatrix, rm)  // const Matrix for all cells
    DEFINE_INTEGRATOR_LF(RVector, rv)      // const scalar for each cells | nodes
    DEFINE_INTEGRATOR_LF(Pos, p)      // const vector for all cells //!calling order!
    DEFINE_INTEGRATOR_LF(PosVector, pv)    // const vector for each cells
    DEFINE_INTEGRATOR_LF(std::vector< RSmallMatrix >, vrm)// const matrix for each cells
    DEFINE_INTEGRATOR_LF(std::vector< RVector >, vrv)// scalar for quadr. on each cells
    DEFINE_INTEGRATOR_LF(std::vector< PosVector >, vpv)// vector for quadr. on each cells
    DEFINE_INTEGRATOR_LF(std::vector< std::vector< RSmallMatrix > >, vvrm)// mat for quadr. on each cells
    #undef DEFINE_INTEGRATOR_LF

    #define DEFINE_INTEGRATOR_LF_NODE(A_TYPE, varname)                         \
        /*! R = \int_mesh this * f \d d mesh and R = RVector(dof) and          \
            f = A_TYPE                                                         \
        */                                                                     \
        void integrate_n(const A_TYPE & varname, RVector & R,                  \
                         const double & alpha=1.0) const;                      \
        /*! Integrate into linear form R = alpha * \int_mesh                   \
        this * f * R \d d mesh and R = RVector(dof) and f = A_TYPE             \
        */                                                                     \
        void integrate_n(const A_TYPE & varname, RVector & R,                  \
                         const RVector & alpha) const;                         \
        /*! Integrate into linear form R = alpha[cellID]*\int_mesh             \
        this * f * R \d d mesh and R = RVector(dof) and f = A_TYPE             \
        */
    DEFINE_INTEGRATOR_LF_NODE(RVector, rv)      // const scalar for each node
    DEFINE_INTEGRATOR_LF_NODE(PosVector, pv)    // const vector for each node
    DEFINE_INTEGRATOR_LF_NODE(std::vector< RSmallMatrix >, vrm)// const matrix for each node
    DEFINE_INTEGRATOR_LF_NODE(std::vector< RVector >, vrv)// scalar for quadr. on each node
    DEFINE_INTEGRATOR_LF_NODE(std::vector< PosVector >, vpv)// vector for quadr. on each node
    DEFINE_INTEGRATOR_LF_NODE(std::vector< std::vector< RSmallMatrix > >, vvrm)// mat for quadr. on each node
    #undef DEFINE_INTEGRATOR_LF_NODE

    #define DEFINE_INTEGRATOR(A_TYPE) \
        /*! Integrate into bilinear form R = \int_mesh this * f * R \d mesh an\
        R = RSparseMapMatrix(dof, dof) and \
            f = A_TYPE \
        */ \
        void integrate(const ElementMatrixMap & R, const A_TYPE & f, \
                       SparseMatrixBase & A, const double & scale=1.0) const; \
        RSparseMapMatrix integrate(const ElementMatrixMap & R, \
                                   const A_TYPE & f, const double & scale=1.0) const; \

    DEFINE_INTEGRATOR(double)   // const scalar for all cells
    DEFINE_INTEGRATOR(RSmallMatrix)  // const Matrix for all cells
    DEFINE_INTEGRATOR(RVector)      // const scalar for each cells | nodes
    DEFINE_INTEGRATOR(std::vector< RSmallMatrix >)// const matrix for each cells
    DEFINE_INTEGRATOR(std::vector< RVector >)// scalar for quadr. on each cells
    DEFINE_INTEGRATOR(std::vector< PosVector >)// vector for quadr. on each cells
    DEFINE_INTEGRATOR(std::vector< std::vector< RSmallMatrix > >)// mat for quadr. on each cells

    #undef DEFINE_INTEGRATOR

    // To avoid ambiguity between Pos|RVector and PosVector|Matrix we define a slightly different interface with differnt keyword names and python have to be decide what is called
    #define DEFINE_INTEGRATOR(A_TYPE) \
        /*! Integrate into bilinear form R = \int_mesh this * f * R \d mesh an\
        R = RSparseMapMatrix(dof, dof) and \
            f = A_TYPE \
        */ \
        void integrate(const ElementMatrixMap & R, const A_TYPE & v, \
                       SparseMatrixBase & A, const double & scale=1.0) const; \
        RSparseMapMatrix integrate(const ElementMatrixMap & R, \
                                   const A_TYPE & v, const double & scale=1.0) const; \

    DEFINE_INTEGRATOR(Pos)      // const scalar for each cells (u * (pos*v))
    DEFINE_INTEGRATOR(PosVector) // const Pos (u * (pos*v)) for each cells
    #undef DEFINE_INTEGRATOR

    /*! Create generic bilinear form */
    void dot(const ElementMatrixMap & B, ElementMatrixMap & ret) const;

    #define DEFINE_ASSEMBLER(A_TYPE) \
        /*! Assemble linear form with non continuous properties. */ \
        void assemble(const A_TYPE & f, RVector & R, \
                      const double & scale=1.0) const; \
        /*! Assemble bilinear form with non continuous properties. */ \
        void assemble(const A_TYPE & f, SparseMatrixBase & A, \
                      const double & scale=1.0) const; \

    DEFINE_ASSEMBLER(double)   // const scalar for all cells
    DEFINE_ASSEMBLER(RSmallMatrix)  // const Matrix for all cells
    DEFINE_ASSEMBLER(RVector)  // const scalar for each cell
    DEFINE_ASSEMBLER(std::vector< RSmallMatrix >)// const matrix for each cell
    #undef DEFINE_ASSEMBLER

    // To avoid ambiguity between Pos|RVector and PosVector|Matrix we define a slightly different interface with differnt keyword names and python have to be decide what is called
    #define DEFINE_ASSEMBLER(A_TYPE) \
        /*! Assemble linear form with non continuous properties. */ \
        /*! Assemble bilinear form with non continuous properties. */ \
        void assemble(const A_TYPE & v, RVector & R, \
                      const double & scale=1.0) const; \
        void assemble(const A_TYPE & v, SparseMatrixBase & A, \
                      const double & scale=1.0) const; \

    DEFINE_ASSEMBLER(Pos)  // const Pos for all cells
    DEFINE_ASSEMBLER(std::vector< Pos >)  // const Pos for each cell

    #undef DEFINE_ASSEMBLER

    /*!Performance issues here. Don't use in python until
    return_value_policy< bp::copy_const_reference > has been changed. */
    const std::vector< ElementMatrix < double > > & mats() const;

    ElementMatrix < double > * pMat(Index i){ return & mats_[i]; }

    void collectQuadraturePoints() const;

    std::vector < PosVector > & quadraturePoints() const;

    std::vector < PosVector > * quadraturePointsRef(){
        return &quadrPnts_;
    }
    PosVector entityCenters() const;

    void fillSparsityPattern(RSparseMatrix & R) const ;

    void fillSparsityPattern(RSparseMatrix & R, const ElementMatrixMap & A) const;

    /*!Calculate copy of this + B, depending on requested dim. */
    void add(const ElementMatrixMap & B, ElementMatrixMap & ret, Index dim=1, double b=1.0) const;

    void add(Index row, const ElementMatrix < double > & Ai);

    //TODO .. check if its the same like mult(a-b, m-n))
    RVector mult(const RVector & a, const RVector & b,
                 const RVector & m, const RVector & n) const;

    /*! Return (S_i * a) * b for all i*/
    RVector mult(const RVector & a, const RVector & b) const;

    Index size() const { return mats_.size();}

    Index rows() const { return rows_; }
    Index cols() const { return cols_; }


    Index dof() const { return this->dofA_; }
    Index dofA() const { return this->dofA_; }
    Index dofB() const { return this->dofB_; }
    void setDof(Index a) { this->dofA_ = a; }
    void setDof(Index a, Index b) { this->dofA_ = a; this->dofB_ = b;}
    inline void setDofs(Index nCoeff, Index dofPerCoeff, Index dofOffset){
        nCoeff_ = nCoeff;
        dofPerCoeff_ = dofPerCoeff;
        dofOffset_ = dofOffset;
        this->setDof(dofPerCoeff_ * nCoeff_ + dofOffset_);
    }

    inline Index nCoeff() const { return nCoeff_;}
    inline Index dofPerCoeff() const { return dofPerCoeff_;}
    inline Index dofOffset() const { return dofOffset_;}

protected:
    std::vector< ElementMatrix < double > > mats_;
    Index nCoeff_;
    Index dofPerCoeff_;
    Index dofOffset_;

    mutable std::vector < PosVector > quadrPnts_; // cache Q for mats_

    std::vector< RSmallMatrix > mat_;
    std::vector< IndexArray > _ids;
    std::vector< Index > row_;

    Index rows_;
    Index cols_;

    Index dofA_;
    Index dofB_;
};


DLLEXPORT void
createUMap(const Mesh & mesh, Index order,
           ElementMatrixMap & ret,
           Index nCoeff=1, Index dofOffset=0);

DLLEXPORT ElementMatrixMap
createUMap(const Mesh & mesh, Index order,
           Index nCoeff=1, Index dofOffset=0);

DLLEXPORT void
createdUMap(const Mesh & mesh, Index order,
            ElementMatrixMap & ret,
            bool elastic, bool div, bool kelvin,
            Index nCoeff=1, Index dofOffset=0);

DLLEXPORT ElementMatrixMap
createdUMap(const Mesh & mesh, Index order,
            bool elastic, bool div, bool kelvin,
            Index nCoeff=1, Index dofOffset=0);

DLLEXPORT void
createIdentityMap(const Mesh & mesh, Index order,
                  ElementMatrixMap & ret,
                  Index nCoeff=1, Index dofOffset=0);

DLLEXPORT ElementMatrixMap
createIdentityMap(const Mesh & mesh, Index order,
                  Index nCoeff, Index dofOffset);

/*! Return symmetrized copy of A as 0.5*(A + A.T).
ATM. Only for gradients without Voigt or Kelvin notation.
*/
DLLEXPORT ElementMatrixMap sym(const ElementMatrixMap & A);
/*! Fill ret with symmetrized copy of A as 0.5*(A + A.T).
ATM. Only for gradients without Voigt or Kelvin notation.
*/
DLLEXPORT void sym(const ElementMatrixMap & A, ElementMatrixMap & ret);
/*! Return trace of A.
ATM. Only for gradients without Voigt or Kelvin notation.
*/
DLLEXPORT ElementMatrixMap trace(const ElementMatrixMap & A);
/*! Fill ret with trace of A.
*/
DLLEXPORT void trace(const ElementMatrixMap & A, ElementMatrixMap & ret);

DLLEXPORT void createUMap_(const Mesh & mesh, Index order, ElementMatrixMap & ret,
                Index nCoeff, Index dofOffset);

DLLEXPORT void createUMap0_(const Mesh & mesh, Index order, ElementMatrixMap & ret,
                Index nCoeff, Index dofOffset);
DLLEXPORT void createUMap1_(const Mesh & mesh, Index order, ElementMatrixMap & ret,
                Index nCoeff, Index dofOffset);
DLLEXPORT void createUMap2_(const Mesh & mesh, Index order, ElementMatrixMap & ret,
                Index nCoeff, Index dofOffset);

DLLEXPORT void testEvalEmap(const ElementMatrixMap & eMap,
                            const stdVectorR3Vector & v, stdVectorRVector & ret);

} // namespace GIMLI{
