//cp ../gimli/core/python/explicit_instances.cpp core/python/generated/explicit_instances.cpp && make pg
//cp ../gimli/core/python/explicit_instances.cpp core/python/generated/explicit_instances.cpp && make pg && python -c 'import pygimli as pg' | grep "undefined symbol"

#include <typeinfo>

#include "tuples.hpp"
#include <Python.h>

#include "baseentity.h"
#include "blockmatrix.h"
#include "curvefitting.h"
#include "cholmodWrapper.h"
#include "datacontainer.h"
#include "gimli.h"
#include "gravimetry.h"
#include "dc1dmodelling.h"
#include "elementmatrix.h"
#include "em1dmodelling.h"
#include "exitcodes.h"
#include "expressions.h"
#include "interpolate.h"
#include "integration.h"
#include "inversion.h"
#include "ldlWrapper.h"
#include "line.h"
#include "linSolver.h"
#include "matrix.h"
#include "memwatch.h"
#include "mesh.h"
#include "meshentities.h"
#include "meshgenerators.h"
#include "modellingbase.h"
#include "node.h"
#include "numericbase.h"
#include "plane.h"
#include "platform.h"
#include "pos.h"
#include "polynomial.h"
#include "quaternion.h"
#include "regionManager.h"
#include "shape.h"
#include "solver.h"
#include "solverWrapper.h"
#include "sparsematrix.h"
#include "spline.h"
#include "stopwatch.h"
#include "trans.h"
#include "triangleWrapper.h"
#include "ttdijkstramodelling.h"
#include "vector.h"
#include "vectortemplates.h"

namespace GIMLI{

#define DEFINE_COMPARE_OPERATOR__(OP) \
template BVector operator OP (const GIMLI::RVector & v1, const GIMLI::RVector & v2); \
template BVector operator OP (const GIMLI::IVector & v1, const GIMLI::IVector & v2); \

DEFINE_COMPARE_OPERATOR__(<)
DEFINE_COMPARE_OPERATOR__(<=)
DEFINE_COMPARE_OPERATOR__(>=)
DEFINE_COMPARE_OPERATOR__(>)
#undef DEFINE_COMPARE_OPERATOR__

#define DEFINE_COMPARE_OPERATOR__(OP) \
template BVector operator OP (const std::vector < GIMLI::SIndex > & vec, const GIMLI::SIndex & v); \
template BVector operator OP (const std::vector < GIMLI::Index > & vec, const GIMLI::Index & v); \
template BVector operator OP (const GIMLI::RVector & vec, const double & v); \
template BVector operator OP (const GIMLI::RVector & vec, int v); \
template BVector operator OP (const GIMLI::IVector & vec, int v); \

DEFINE_COMPARE_OPERATOR__(<)
DEFINE_COMPARE_OPERATOR__(<=)
DEFINE_COMPARE_OPERATOR__(>=)
DEFINE_COMPARE_OPERATOR__(==)
DEFINE_COMPARE_OPERATOR__(!=)
DEFINE_COMPARE_OPERATOR__(>)
#undef DEFINE_COMPARE_OPERATOR__



    template class Vector< bool >;
    template class Vector< double >;
    template class Vector< GIMLI::RVector3 >;
    template class Vector< GIMLI::SIndex >;
    template class Vector< GIMLI::Complex >;

    template class VectorIterator< bool >;
    template class VectorIterator< double >;
    template class VectorIterator< GIMLI::RVector3 >;
    template class VectorIterator< GIMLI::SIndex >;
    template class VectorIterator< GIMLI::Index >;
    template class VectorIterator< GIMLI::Complex >;

    template class Matrix< double >;
    template class Matrix< std::complex< double > >;
    template class Matrix3< double >;

    template class BlockMatrix< double >;
    template class Quaternion< double >;

    template class PolynomialElement< double >;
    template class PolynomialFunction< double >;

    template std::ostream & operator << (std::ostream & str, const PolynomialFunction < double > & p);
    template PolynomialFunction < double > operator - (const PolynomialFunction < double > & f);
    template PolynomialFunction < double > operator - (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    template PolynomialFunction < double > operator + (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    template PolynomialFunction < double > operator * (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    template PolynomialFunction < double > operator * (const PolynomialFunction < double > & f, const double & val);
    template PolynomialFunction < double > operator * (const double & val, const PolynomialFunction < double > & f);
    template PolynomialFunction < double > operator + (const PolynomialFunction < double > & f, const double & val);
    template PolynomialFunction < double > operator + (const double & val, const PolynomialFunction < double > & f);

    template std::vector < PolynomialFunction < double > >
        createPolynomialShapeFunctions(const Shape & ent, uint nCoeff, bool pascale, bool serendipity, const RVector &);
    template std::vector < PolynomialFunction < double > >
        createPolynomialShapeFunctions(const MeshEntity & ent, uint nCoeff, bool pascale, bool serendipity, const RVector &);

    template double besselI0< double >(const double & x);
    template double besselI1< double >(const double & x);
    template double besselK0< double >(const double & x);
    template double besselK1< double >(const double & x);

    template class Trans< RVector >;
    template class TransLinear< RVector >;
    template class TransLin< RVector >;
    template class TransPower< RVector >;
    template class TransLog< RVector >;
    template class TransLogLU< RVector >;
    template class TransCotLU< RVector >;

#define DEFINE_XVECTOR_STUFF__(VEC) \
template bool haveInfNaN(const VEC & v); \
template BVector isInf(const VEC & vec); \
template BVector isNaN(const VEC & vec); \
template BVector isInfNaN(const VEC & vec); \

DEFINE_XVECTOR_STUFF__(CVector)
DEFINE_XVECTOR_STUFF__(BVector)
DEFINE_XVECTOR_STUFF__(IVector)
DEFINE_XVECTOR_STUFF__(RVector) //RVector last since auto rhs conversion will fail else
#undef DEFINE_XVECTOR_STUFF__

    template RVector fliplr(const RVector & a);
    template RVector round(const RVector & a, double tol);
    template RVector increasingRange(const double & first, const double & last, Index n);

    template Pos & Pos::transform(const Matrix < double > & mat);

    template class SparseMatrix< double >;
    template class SparseMatrix< GIMLI::Complex >;

    template RSparseMatrix operator + (const RSparseMatrix & A, const RSparseMatrix & B);
    template RSparseMatrix operator - (const RSparseMatrix & A, const RSparseMatrix & B);
    template RSparseMatrix operator * (const RSparseMatrix & A, const double & b);
    template RSparseMatrix operator * (const double & b, const RSparseMatrix & A);

    template CSparseMatrix operator + (const CSparseMatrix & A, const CSparseMatrix & B);
    template CSparseMatrix operator - (const CSparseMatrix & A, const CSparseMatrix & B);
    template CSparseMatrix operator * (const GIMLI::Complex & b, const CSparseMatrix & A);
    template CSparseMatrix operator * (const CSparseMatrix & A, const GIMLI::Complex & b);

    template class ElementMatrix< double >;
    // template std::ostream & operator << (std::ostream & str,
    //                                      const ElementMatrix< double > & p);

    template std::vector< Index > unique(const std::vector < Index > & a);
    template std::vector< SIndex > unique(const std::vector < SIndex > & a);
    template RVector unique(const RVector & a);
    // template IndexArray unique(const IndexArray & a);
    template IVector unique(const IVector & a);

    template RVector sort(const RVector & a);
    // template IndexArray sort(const IndexArray & a);
    template IVector sort(const IVector & a);
    template std::vector< Index > sort(const std::vector < Index > & a);
    template std::vector< SIndex > sort(const std::vector < SIndex > & a);

    template RVector pow(const RVector & a, double power);
    template RVector pow(const RVector & a, int power);
    template RVector pow(const RVector & a, const RVector & power);
    template RVector cat(const RVector & a, const RVector & b);

    template RVector real(const CVector & a);
    template RVector imag(const CVector & a);
    template RVector angle(const CVector & a);
    template RVector phase(const CVector & a);
    template RVector abs(const CVector & a);
    template CVector conj(const CVector & a);

    template RMatrix real(const CMatrix & a);
    template RMatrix imag(const CMatrix & a);

    template double det(const RMatrix & a);
    template double det(const RMatrix3 & a);

    template double min(const RVector & v);
    template double max(const RVector & v);

    template double rms(const RVector & a);
    template double rms(const RVector & a, const RVector & b);
    template double rrms(const RVector & a, const RVector & b);

    template double norm(const RVector & a);
    template double normlp(const RVector & a, int p);
    template double norml1(const RVector & a);
    template double norml2(const RVector & a);
    template double normlInfinity(const RVector & a);

    template double norm(const CVector & a);
    template double normlp(const CVector & a, int p);
    template double norml1(const CVector & a);
    template double norml2(const CVector & a);
    template double normlInfinity(const CVector & a);

    template double euclideanNorm(const RVector & a);
    template double stdDev(const RVector & v);
    template double median(const RVector & a);
    template double mean(const RVector & a);
    template double arithmeticMean(const RVector & v);
    template double geometricMean(const RVector & v);
    template double harmonicMean(const RVector & a);

    template double dot(const RVector & v1, const RVector & v2);

    template void rand(RVector & vec, double min = 0.0, double max = 1.0);
    template void randn(RVector & vec);

    template double degToRad(const double & deg);
    template double radToDeg(const double & rad);

    template RVector3 degToRad(const RVector3 & deg);
    template RVector3 radToDeg(const RVector3 & rad);

    template void sort(const RVector & a, RVector & b, IndexArray & idx);
    template IndexArray sortIdx(const RVector & a);

    template std::string save(const RVector &v,
                              const std::string & fname,
                              IOFormat format=Ascii);
    template bool load(RVector &v, const std::string & fname, IOFormat format = Ascii,
                        bool verbose = true);
    template bool load(RMatrix & A, const std::string & filename);

    template bool saveMatrixCol(const RMatrix & A, const std::string & filename);
    template bool saveMatrixCol(const RMatrix & A, const std::string & filename,
                                 const std::string & comments);

    template bool loadMatrixCol(RMatrix & A, const std::string & filename);
    template bool loadMatrixCol(RMatrix & A, const std::string & filename,
                                 std::vector < std::string > & comments);

    template bool saveMatrixRow(const RMatrix & A, const std::string & filename);
    template bool saveMatrixRow(const RMatrix & A, const std::string & filename,
                                 const std::string & comments);

    template bool loadMatrixRow(RMatrix & A, const std::string & filename);
    template bool loadMatrixRow(RMatrix & A, const std::string & filename,
                                std::vector < std::string > & comments);

    template std::set< Node * > commonNodes(const std::set < Boundary * > &);
    template std::set< Node * > commonNodes(const std::set < Cell * > &);

    void mult_Mv_add_sv(RVector & rhs, const RSparseMapMatrix & M, const RVector & v1,
                   const double s, const RVector & v2){
        mult(M, v1, rhs);
        rhs += (s*v2);
    }
    void mult_Mv_add_sv(RVector & rhs, const RSparseMatrix & M, const RVector & v1,
                        const double s, const RVector & v2){

        // mult(M, v1, rhs);
        // rhs += (s*v2);
        rhs.resize(M.rows());
        // 60% faster, TODO: check why
        double si = 0.0;
        int dataIdx = M.vecColPtr()[0];
        for (Index i = 0, iMax = M.rows(); i < iMax; i++){
            si = 0.0;
            for (; dataIdx < M.vecColPtr()[i+1]; ++dataIdx){
                si += M.vecVals()[dataIdx] * v1[M.vecRowIdx()[dataIdx]];
            }
            rhs[i] = si + s*v2[i];
        }

//#pragma omp parallel
{
        // for (Index i = 0, iMax = M.rows(); i < iMax; i++){
        //     double si = 0.0;
        //     for (int j = M.vecColPtr()[i], jMax=M.vecColPtr()[i+1];
        //         j < jMax; j ++){
        //         si += M.vecVals()[j] * v1[M.vecRowIdx()[j]];
        //     }
        //     rhs[i] = si + s*v2[i];
        // }
}

    }

// template class PolynomialElement< double >;
// template class BlockMatrix< double >;
// template class TransLinear< RVector >;
// template class TransLin< RVector >;
// template class TransPower< RVector >;
// template class Matrix< std::complex< double > >;

// //template std::vector< Index > unique(const std::vector < Index > & a);
// template std::vector< SIndex > unique(const std::vector < SIndex > & a);

// template bool loadMatrixCol(RMatrix & A, const std::string & filename);
// template class VectorIterator< double >;
// template class SparseMatrix< double >;
// template class SparseMatrix< GIMLI::Complex >;
// template class Quaternion< double >;
// template class Vector< bool >;
// template class Vector< double >;
// template class Vector< Complex >;
// template class Vector< SIndex >;
// template class Vector< GIMLI::RVector3 >;

// #define DEFINE_XVECTOR_STUFF__(VEC) \
// template bool haveInfNaN(const VEC & v); \
// template BVector isInf(const VEC & vec); \
// template BVector isNaN(const VEC & vec); \
// template BVector isInfNaN(const VEC & vec); \

// DEFINE_XVECTOR_STUFF__(CVector)
// DEFINE_XVECTOR_STUFF__(BVector)
// DEFINE_XVECTOR_STUFF__(IVector)
// DEFINE_XVECTOR_STUFF__(RVector) //RVector last since auto rhs conversion will fail else
// #undef DEFINE_XVECTOR_STUFF__

}