#ifndef PYTHON_PYGIMLI__H
#define PYTHON_PYGIMLI__H

//See best practices section in Py++ documentation

#ifdef __clang__
// probably clang bug in combination with gcc > 7
#ifndef __STRICT_ANSI__
    #define __STRICT_ANSI__
#endif
//#define __float128 long double
//typedef struct { long double x, y; } __float128;
#endif

//#define PYTEST
#ifdef PYTEST

 #include <iostream>
 #include <fstream>
//  #include "stopwatch.h"
//  #include "vector.h"

#ifndef PYGIMLI_CAST
    #include "pos.h"
    #include "matrix.h"
    #include "elementmatrix.h"
#endif

namespace GIMLI{
    inline void tmp(){
         std::cout << "tmp"<< std::endl;  
    }
    // class TMP {
    //         public:
    //         void a(){std::cout << "TMP"<< std::endl;}  
    // };

    // extern template class Vector< double >;
//     template class VectorIterator< double >;
//     template class Vector< Complex >;
//     template class Vector< int >;
//     template class BlockMatrix< double >;

    inline void ___instantiation___(){
        // sizeof(::GIMLI::Index *);
		// sizeof(::GIMLI::Index &);
		// sizeof(::GIMLI::Index);
		// sizeof(::GIMLI::IndexArray *);
		// sizeof(::GIMLI::IndexArray &);
		// sizeof(::GIMLI::IndexArray);
		// sizeof(int *);
        // sizeof(int &);
        // sizeof(int);
        // sizeof(long unsigned int *);
        // sizeof(long unsigned int &);
        // sizeof(long unsigned int);
		// sizeof(long long unsigned int *);
        // sizeof(long long unsigned int &);
        // sizeof(long long unsigned int);
		// sizeof(long long int *);
        // sizeof(long long int &);
        // sizeof(long long int);
        // sizeof(double *);
        // sizeof(double);
        // sizeof(bool);
        // sizeof(double &);
     }

} // namespace GIMLI

//     typedef std::complex< double >                  Complex;
//     typedef GIMLI::Vector< double >                 RVector;
//     // typedef GIMLI::Vector< bool >                 BVector;
//     // typedef GIMLI::Vector< long >                 IVector;
    
//     typedef GIMLI::VectorIterator< double >          RVectorIter;// avoid doubs
//     typedef GIMLI::VectorIterator< bool >            BVectorIter;// avoid doubs
//     typedef GIMLI::VectorIterator< long >            IVectorIter;// avoid doubs
    
//     // typedef GIMLI::VectorIterator< Complex >              CVectorIter;
//     // typedef GIMLI::BlockMatrix< double >                 RBlockMatrix;

//     }
// } //pyplusplus::aliases

#else // if not PYTEST

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

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
#include "bert/bert.h"
#include "bert/bertMisc.h"
#include "bert/bertJacobian.h"
#include "bert/datamap.h"
#include "bert/dcfemmodelling.h"
#include "bert/bertDataContainer.h"
#include "bert/electrode.h"
// #include "inversionRollalong.h"
// #include "eameshwrapper.h"
// #include "matrixTemplates.h"

namespace GIMLI{

#define DEFINE_PY_VEC_OPERATOR__(OP) \
    inline IndexArray operator OP (const IndexArray & a, const IndexArray & b){\
        IndexArray ret(a);   ret OP##=b; return ret; } \
    inline IndexArray operator OP (const IndexArray & a, Index b){ \
        IndexArray ret(a);   ret OP##=b; return ret; } \
    inline IndexArray operator OP (const IndexArray & a, int b){ \
        IndexArray ret(a);   ret OP##=(Index)abs(b); return ret; } \
    inline IndexArray operator OP (Index a, const IndexArray & b){ \
        IndexArray ret(b.size()); for (Index i = 0; i < b.size(); i ++) ret[i] = a OP b[i]; return ret; } \
    inline IndexArray operator OP (int a, const IndexArray & b){ \
        IndexArray ret(b.size()); for (Index i = 0; i < b.size(); i ++) ret[i] = (Index)abs(a) OP b[i]; return ret; } \
    inline RVector operator OP (const RVector & a, const RVector & b){ \
        RVector ret(a);   ret OP##=b; return ret; }                           \
    inline RVector operator OP (const double & a, const RVector & b){ \
        RVector ret(b.size()); for (Index i = 0; i < b.size(); i ++) ret[i] = a OP b[i]; return ret; } \
    inline RVector operator OP (const RVector & a, const double & b){ \
        RVector ret(a);   ret OP##=b; return ret; }                           \
    inline IVector operator OP (const IVector & a, const IVector & b){ \
        IVector ret(a);   ret OP##=b; return ret; }                           \
    inline IVector operator OP (const IVector & a, SIndex b){ \
        IVector ret(a);   ret OP##=b; return ret; }                           \
    inline IVector operator OP (SIndex a, const IVector & b){ \
        IVector ret(b.size()); \
        for (Index i = 0; i < b.size(); i ++) ret[i] = a OP b[i]; return ret; } \
    inline CVector operator OP (const CVector & a, const CVector & b){ \
        CVector c(a);   c OP##=b; return c; }                           \
    inline CVector operator OP (const CVector & a, const Complex & b){ \
        CVector c(a);   c OP##=b; return c; }                           \
    inline CVector operator OP (const CVector & a, const RVector & b){ \
        CVector c(a);   c OP##=toComplex(b); return c; }                        \
    inline CVector operator OP (const CVector & a, const double  & b){ \
        CVector c(a);   c OP##=Complex(b); return c; }                         \
    inline CVector operator OP (const Complex & v, const CVector & b){ \
        CVector c(b.size());   for (Index i = 0; i < c.size(); i ++) c[i] = v OP b[i]; return c; } \
    inline CVector operator OP (const double  & v, const CVector & b){ \
        CVector c(b.size());   for (Index i = 0; i < c.size(); i ++) c[i] = v OP b[i]; return c; } \
    inline CVector operator OP (const RVector & a, const CVector & b){ \
        CVector c(b.size());   for (Index i = 0; i < c.size(); i ++) c[i] = a[i] OP b[i]; return c; } \
    inline RMatrix operator OP (const RMatrix & a, const double & b){ \
        RMatrix ret(a);   ret OP##=b; return ret; } \
    inline RMatrix operator OP (const double & a, const RMatrix & b){ \
        RMatrix ret(b.rows()); for (Index i = 0; i < b.rows(); i ++) ret[i] = a OP b[i]; return ret; } \
    inline RMatrix operator OP (const RMatrix & a, const RMatrix & b){ \
         RMatrix ret(a);   ret OP##=b; return ret; } \

DEFINE_PY_VEC_OPERATOR__(+)
DEFINE_PY_VEC_OPERATOR__(-)
DEFINE_PY_VEC_OPERATOR__(*)
DEFINE_PY_VEC_OPERATOR__(/)
#undef DEFINE_PY_VEC_OPERATOR__

#define DEFINE_PY_VEC_UNARY_OPERATOR__(OP, FUNCT)                      \
    inline RVector OP (const RVector & a) { \
        RVector tmp(a.size()); for (uint i = 0; i < a.size(); i ++) tmp[i] = FUNCT()(a[i]); return tmp; } \

DEFINE_PY_VEC_UNARY_OPERATOR__(abs,   ABS_)
DEFINE_PY_VEC_UNARY_OPERATOR__(acot,  ACOT)
DEFINE_PY_VEC_UNARY_OPERATOR__(atan,  ATAN)
DEFINE_PY_VEC_UNARY_OPERATOR__(cos,   COS)
DEFINE_PY_VEC_UNARY_OPERATOR__(cot,   COT)
DEFINE_PY_VEC_UNARY_OPERATOR__(exp,   EXP)
DEFINE_PY_VEC_UNARY_OPERATOR__(exp10, EXP10)
DEFINE_PY_VEC_UNARY_OPERATOR__(fabs,  ABS_)
DEFINE_PY_VEC_UNARY_OPERATOR__(log,   LOG)
DEFINE_PY_VEC_UNARY_OPERATOR__(log10, LOG10)
DEFINE_PY_VEC_UNARY_OPERATOR__(sign,  SIGN)
DEFINE_PY_VEC_UNARY_OPERATOR__(sin,   SIN)
DEFINE_PY_VEC_UNARY_OPERATOR__(sqrt,  SQRT)
DEFINE_PY_VEC_UNARY_OPERATOR__(square, SQR)
DEFINE_PY_VEC_UNARY_OPERATOR__(tan,   TAN)
DEFINE_PY_VEC_UNARY_OPERATOR__(tanh,  TANH)
#undef DEFINE_PY_VEC_UNARY_OPERATOR__

    inline RVector mag(const CVector & a){return abs(a);}


    //** maybe better to move these instantiation into libgimli, but why should the lib have unused template symbols??
    // explicit instantiaion leads to duplicate symbols for architecture arm64 
    // if the symbols allready exists, extern to mark them

#define DEFINE_COMPARE_OPERATOR__(OP) \
extern template BVector operator OP (const std::vector < GIMLI::SIndex > & vec, const GIMLI::SIndex & v); \
extern template BVector operator OP (const std::vector < GIMLI::Index > & vec, const GIMLI::Index & v); \

DEFINE_COMPARE_OPERATOR__(<)
DEFINE_COMPARE_OPERATOR__(<=)
DEFINE_COMPARE_OPERATOR__(>=)
DEFINE_COMPARE_OPERATOR__(==)
DEFINE_COMPARE_OPERATOR__(!=)
DEFINE_COMPARE_OPERATOR__(>)

    extern template class Vector< bool >;
    extern template class Vector< double >;
    extern template class Vector< GIMLI::RVector3 >;
    extern template class Vector< GIMLI::SIndex >;
    extern template class Vector< GIMLI::Complex >;

    extern template class VectorIterator< bool >;
    extern template class VectorIterator< double >;
    extern template class VectorIterator< GIMLI::RVector3 >;
    extern template class VectorIterator< GIMLI::SIndex >;
    extern template class VectorIterator< GIMLI::Index >;
    extern template class VectorIterator< GIMLI::Complex >;

    extern template class Matrix< double >;
    extern template class Matrix< std::complex< double > >;
    extern template class Matrix3< double >;

    extern template class BlockMatrix< double >;
    extern template class Quaternion< double >;

    extern template class PolynomialElement< double >;
    extern template class PolynomialFunction< double >;

    extern template std::ostream & operator << (std::ostream & str, const PolynomialFunction < double > & p);
    extern template PolynomialFunction < double > operator - (const PolynomialFunction < double > & f);
    extern template PolynomialFunction < double > operator - (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    extern template PolynomialFunction < double > operator + (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    extern template PolynomialFunction < double > operator * (const PolynomialFunction < double > & f, const PolynomialFunction < double > & g);
    extern template PolynomialFunction < double > operator * (const PolynomialFunction < double > & f, const double & val);
    extern template PolynomialFunction < double > operator * (const double & val, const PolynomialFunction < double > & f);
    extern template PolynomialFunction < double > operator + (const PolynomialFunction < double > & f, const double & val);
    extern template PolynomialFunction < double > operator + (const double & val, const PolynomialFunction < double > & f);

    extern template std::vector < PolynomialFunction < double > >
        createPolynomialShapeFunctions(const Shape & ent, uint nCoeff, bool pascale, bool serendipity, const RVector &);
    extern template std::vector < PolynomialFunction < double > >
        createPolynomialShapeFunctions(const MeshEntity & ent, uint nCoeff, bool pascale, bool serendipity, const RVector &);

    extern template double besselI0< double >(const double & x);
    extern template double besselI1< double >(const double & x);
    extern template double besselK0< double >(const double & x);
    extern template double besselK1< double >(const double & x);

    extern template class Trans< RVector >;
    extern template class TransLinear< RVector >;
    extern template class TransLin< RVector >;
    extern template class TransPower< RVector >;
    extern template class TransLog< RVector >;
    extern template class TransLogLU< RVector >;
    extern template class TransCotLU< RVector >;

#define DEFINE_XVECTOR_STUFF__(VEC) \
extern template bool haveInfNaN(const VEC & v); \
extern template BVector isInf(const VEC & vec); \
extern template BVector isNaN(const VEC & vec); \
extern template BVector isInfNaN(const VEC & vec); \

DEFINE_XVECTOR_STUFF__(CVector)
DEFINE_XVECTOR_STUFF__(BVector)
DEFINE_XVECTOR_STUFF__(IVector)
DEFINE_XVECTOR_STUFF__(RVector) //RVector last since auto rhs conversion will fail else
#undef DEFINE_XVECTOR_STUFF__

    extern template RVector fliplr(const RVector & a);
    extern template RVector round(const RVector & a, double tol);
    extern template RVector increasingRange(const double & first, const double & last, Index n);

    extern template Pos & Pos::transform(const Matrix < double > & mat);

    extern template class SparseMatrix< double >;
    extern template class SparseMatrix< GIMLI::Complex >;

    extern template RSparseMatrix operator + (const RSparseMatrix & A, const RSparseMatrix & B);
    extern template RSparseMatrix operator - (const RSparseMatrix & A, const RSparseMatrix & B);
    extern template RSparseMatrix operator * (const RSparseMatrix & A, const double & b);
    extern template RSparseMatrix operator * (const double & b, const RSparseMatrix & A);

    extern template CSparseMatrix operator + (const CSparseMatrix & A, const CSparseMatrix & B);
    extern template CSparseMatrix operator - (const CSparseMatrix & A, const CSparseMatrix & B);
    extern template CSparseMatrix operator * (const GIMLI::Complex & b, const CSparseMatrix & A);
    extern template CSparseMatrix operator * (const CSparseMatrix & A, const GIMLI::Complex & b);

    extern template class ElementMatrix< double >;
    // extern template std::ostream & operator << (std::ostream & str,
    //                                      const ElementMatrix< double > & p);

    extern template std::vector< Index > unique(const std::vector < Index > & a);
    extern template std::vector< SIndex > unique(const std::vector < SIndex > & a);
    extern template RVector unique(const RVector & a);
    // extern template IndexArray unique(const IndexArray & a);
    extern template IVector unique(const IVector & a);

    extern template RVector sort(const RVector & a);
    // extern template IndexArray sort(const IndexArray & a);
    extern template IVector sort(const IVector & a);
    extern template std::vector< Index > sort(const std::vector < Index > & a);
    extern template std::vector< SIndex > sort(const std::vector < SIndex > & a);

    extern template RVector pow(const RVector & a, double power);
    extern template RVector pow(const RVector & a, int power);
    extern template RVector pow(const RVector & a, const RVector & power);
    extern template RVector cat(const RVector & a, const RVector & b);

    extern template RVector real(const CVector & a);
    extern template RVector imag(const CVector & a);
    extern template RVector angle(const CVector & a);
    extern template RVector phase(const CVector & a);
    extern template RVector abs(const CVector & a);
    extern template CVector conj(const CVector & a);

    extern template RMatrix real(const CMatrix & a);
    extern template RMatrix imag(const CMatrix & a);

    extern template double det(const RMatrix & a);
    extern template double det(const RMatrix3 & a);

    extern template double min(const RVector & v);
    extern template double max(const RVector & v);

    extern template double rms(const RVector & a);
    extern template double rms(const RVector & a, const RVector & b);
    extern template double rrms(const RVector & a, const RVector & b);

    extern template double norm(const RVector & a);
    extern template double normlp(const RVector & a, int p);
    extern template double norml1(const RVector & a);
    extern template double norml2(const RVector & a);
    extern template double normlInfinity(const RVector & a);

    extern template double norm(const CVector & a);
    extern template double normlp(const CVector & a, int p);
    extern template double norml1(const CVector & a);
    extern template double norml2(const CVector & a);
    extern template double normlInfinity(const CVector & a);

    extern template double euclideanNorm(const RVector & a);
    extern template double stdDev(const RVector & v);
    extern template double median(const RVector & a);
    extern template double mean(const RVector & a);
    extern template double arithmeticMean(const RVector & v);
    extern template double geometricMean(const RVector & v);
    extern template double harmonicMean(const RVector & a);

    extern template double dot(const RVector & v1, const RVector & v2);

    extern template void rand(RVector & vec, double min = 0.0, double max = 1.0);
    extern template void randn(RVector & vec);

    extern template double degToRad(const double & deg);
    extern template double radToDeg(const double & rad);

    extern template RVector3 degToRad(const RVector3 & deg);
    extern template RVector3 radToDeg(const RVector3 & rad);

    extern template void sort(const RVector & a, RVector & b, IndexArray & idx);
    extern template IndexArray sortIdx(const RVector & a);

    extern template bool save(const RVector &v, const std::string & fname, IOFormat format = Ascii);
    extern template bool load(RVector &v, const std::string & fname, IOFormat format = Ascii,
                        bool verbose = true);
    extern template bool load(RMatrix & A, const std::string & filename);

    extern template bool saveMatrixCol(const RMatrix & A, const std::string & filename);
    extern template bool saveMatrixCol(const RMatrix & A, const std::string & filename,
                                 const std::string & comments);

    extern template bool loadMatrixCol(RMatrix & A, const std::string & filename);
    extern template bool loadMatrixCol(RMatrix & A, const std::string & filename,
                                 std::vector < std::string > & comments);

    extern template bool saveMatrixRow(const RMatrix & A, const std::string & filename);
    extern template bool saveMatrixRow(const RMatrix & A, const std::string & filename,
                                 const std::string & comments);

    extern template bool loadMatrixRow(RMatrix & A, const std::string & filename);
    extern template bool loadMatrixRow(RMatrix & A, const std::string & filename,
                                std::vector < std::string > & comments);

    extern template std::set< Node * > commonNodes(const std::set < Boundary * > &);
    extern template std::set< Node * > commonNodes(const std::set < Cell * > &);


//     inline void ___instantiation___(){
//         sizeof(bool);
//         sizeof(bool *);
//         sizeof(bool &);
//         sizeof(char);
//         sizeof(char *);
//         sizeof(char &);
//         sizeof(false);
//         sizeof(true);
//         sizeof(int *);
//         sizeof(int &);
//         sizeof(int);
//         sizeof(uint8 *);
//         sizeof(uint8 &);
//         sizeof(uint8);
//         sizeof(unsigned char);
//         sizeof(long unsigned int *);
//         sizeof(long unsigned int &);
//         sizeof(long unsigned int);
//         sizeof(long *);
//         sizeof(long &);
//         sizeof(long);
// 		sizeof(long long unsigned int *);
// 		sizeof(long long unsigned int &);
// 		sizeof(long long unsigned int);
//         sizeof(double *);
//         sizeof(double &);
//         sizeof(double  );
//         sizeof(GIMLI::Complex *);
//         sizeof(GIMLI::Complex &);
//         sizeof(GIMLI::Complex );
//         sizeof(GIMLI::SparseMatrix<double> &);
//         sizeof(GIMLI::SparseMatrix<double> *);
//         sizeof(GIMLI::SparseMatrix<double>);
// 		sizeof(::GIMLI::Index *);
// 		sizeof(::GIMLI::Index &);
// 		sizeof(::GIMLI::Index);
//         sizeof(::GIMLI::SIndex *);
//         sizeof(::GIMLI::SIndex &);
//         sizeof(::GIMLI::SIndex);
// 		sizeof(::GIMLI::IndexArray *);
// 		sizeof(::GIMLI::IndexArray &);
// 		sizeof(::GIMLI::IndexArray);
//     }

    /*! Temporary workaround until there is as solution for
        http://www.mail-archive.com/cplusplus-sig@python.org/msg00333.html/

        declare and init new function, inject them to original class within __init__.py
    */
    // template void Quaternion< double >::rotMatrix(Matrix < double > & mat) const;
    inline void getRotMatrix__(const Quaternion < double > & q, Matrix < double > & mat){
        q.rotMatrix(mat);
    }

// DEPRECATED
//     inline void interpolate_GILsave__(const Mesh & mesh, const RMatrix & data,
//                                       const R3Vector & pos, RMatrix & iData,
//                                       bool verbose=false){
//         __MS("interpolate_GILsave__ 1")
//         ALLOW_PYTHON_THREADS
//         return interpolate(mesh, data, pos, iData, verbose);
//     }
//     inline RVector interpolate_GILsave__(const Mesh & mesh, const RVector & data,
//                                           const RVector & x, const RVector & y, bool verbose = false){
//         ALLOW_PYTHON_THREADS
//         std::cout << "interpolate_GILsave__ 2" << std::endl;
//         return interpolate(mesh, data, x, y, verbose);
//     }
//
//     inline void interpolate_GILsave__(const Mesh & mesh, const RVector & data,
//                                        const Mesh & pos, RVector & iData, bool verbose = false){
//         ALLOW_PYTHON_THREADS
//         std::cout << "interpolate_GILsave__ 3" << std::endl;
//         return interpolate(mesh, data, pos, iData, verbose);
//     }

} // namespace GIMLI

// needed for castxml caster .. need to check .. the name Complex seems not enrolled in CMatrix generate code
// typedef std::complex< double >                       Complex;

//** define some aliases to avoid insane method names
namespace pyplusplus{ namespace aliases{
    typedef std::complex< double >                       Complex;
    typedef GIMLI::Pos                                   RVector3;

    typedef GIMLI::Vector< bool >                        BVector;
    typedef GIMLI::Vector< double >                      RVector;
    typedef GIMLI::Vector< RVector3 >                    R3Vector;
    typedef GIMLI::Vector< GIMLI::SIndex >               IVector;
    typedef GIMLI::Vector< GIMLI::Index >                IndexArray;
    typedef GIMLI::Vector< Complex >                     CVector;

    typedef GIMLI::VectorIterator< bool >                BVectorIter;
    typedef GIMLI::VectorIterator< double >              RVectorIter;
    typedef GIMLI::VectorIterator< RVector3 >            R3VectorIter;
    typedef GIMLI::VectorIterator< GIMLI::SIndex >       SIVectorIter;
    typedef GIMLI::VectorIterator< GIMLI::Index >        IVectorIter;
    typedef GIMLI::VectorIterator< Complex >             CVectorIter;

    typedef GIMLI::Matrix< double >                      RMatrix;
    typedef GIMLI::Matrix3< double >                     RMatrix3;
    typedef GIMLI::Matrix< std::complex< double > >      CMatrix;

    typedef GIMLI::BlockMatrix< double >                 RBlockMatrix;

    typedef GIMLI::Quaternion< double >                  RQuaternion;

    typedef GIMLI::PolynomialFunction< double >          RPolynomialFunction;
    typedef GIMLI::PolynomialElement< double >           RPolynomialElement;

//     typedef GIMLI::RollalongInSpace< double >            RRollalongInSpace;

    typedef GIMLI::ElementMatrix< double >                   ElementMatrix;
    typedef GIMLI::SparseMapMatrix< int, GIMLI::Index >      ISparseMapMatrix;
    typedef GIMLI::SparseMapMatrix< double, GIMLI::Index >   RSparseMapMatrix;
    typedef GIMLI::SparseMapMatrix< Complex, GIMLI::Index >  CSparseMapMatrix;
    typedef GIMLI::SparseMatrix< int >                       ISparseMatrix;
    typedef GIMLI::SparseMatrix< double >                    RSparseMatrix;
    typedef GIMLI::SparseMatrix< Complex >                   CSparseMatrix;

    typedef GIMLI::Trans< GIMLI::RVector >        RTrans;
    typedef GIMLI::TransLinear< GIMLI::RVector >  RTransLinear;
    typedef GIMLI::TransLin< GIMLI::RVector >     RTransLin;
    typedef GIMLI::TransPower< GIMLI::RVector >   RTransPower;
    typedef GIMLI::TransLog< GIMLI::RVector >     RTransLog;
    typedef GIMLI::TransLogLU< GIMLI::RVector >   RTransLogLU;
    typedef GIMLI::TransCotLU< GIMLI::RVector >   RTransCotLU;
    typedef GIMLI::TransCumulative < GIMLI::RVector > RTransCumulative;

    typedef GIMLI::Singleton< GIMLI::ShapeFunctionCache >   SingletonShapeFunction;
    typedef GIMLI::Singleton< GIMLI::IntegrationRules >     SingletonIntegrationsRules;
    typedef GIMLI::Singleton< GIMLI::MemWatch >             SingletonMemWatch;
    typedef GIMLI::Function< double, double >               FunctionDD;
//     typedef GIMLI::Variable< GIMLI::XAxis >                 VariableXAxis;
//     typedef GIMLI::Variable< GIMLI::YAxis >                 VariableYAxis;
//     typedef GIMLI::Variable< GIMLI::ZAxis >                 VariableZAxis;

    typedef GIMLI::Expr< GIMLI::ExprIdentity >              ExprExprIdent;

    typedef std::vector< std::string >                  stdVectorString;
    typedef std::vector< int >                          stdVectorI;
    typedef std::vector< GIMLI::SIndex >                stdVectorSIndex;
    typedef std::vector< GIMLI::Index >                 stdVectorIndex;

    //typedef std::vector< long int >                     stdVectorLI;
    //typedef std::vector< double >                       stdVectorR;
    //typedef std::vector< std::complex < double > >      stdVectorC;
    //typedef std::vector< unsigned int >                 stdVectorUI;

    typedef std::vector< std::pair< GIMLI::Index, GIMLI::Index > > stdVectorPairLongLong;
    typedef std::vector< std::pair< unsigned int, unsigned int> >      stdVectorPairUintUint;

    typedef std::vector< GIMLI::Vector< double > >      stdVectorRVector;
    typedef std::vector< GIMLI::Vector< bool > >        stdVectorBVector;
    typedef std::vector< GIMLI::Vector< GIMLI::SIndex > >      stdVectorSIndexVector;
    typedef std::vector< GIMLI::Matrix< double > >      stdVectorRMatrix;
    typedef std::vector< GIMLI::MatrixBase * >          stdpMatrixBase;
    typedef std::vector< GIMLI::RVector3 >              stdVectorRVector3;
    typedef std::vector< GIMLI::R3Vector >              stdVectorR3Vector;
    typedef std::vector< GIMLI::RMatrix3 >              stdVectorMatrix3;
    
    typedef std::vector< GIMLI::PolynomialElement<double> > stdVectorPolynomialElementR;
    typedef std::vector< GIMLI::PolynomialFunction< double > > stdVectorPolynomialFunctionR;

    typedef std::vector< GIMLI::Trans < GIMLI::Vector < double > > * > stdVectorTrans;
    typedef std::vector< GIMLI::Vector< std::complex< double > > > stdVectorCVector;
    typedef std::vector< GIMLI::Vector< std::complex< double > > > stdVectorCVector;
    typedef std::vector< std::vector< GIMLI::Matrix< double > > > stdVectorMatrixVector;

    typedef std::vector< GIMLI::Boundary * >         stdVectorBounds;
    typedef std::vector< GIMLI::Cell * >             stdVectorCells;
    typedef std::vector< GIMLI::Node * >             stdVectorNodes;
    typedef std::vector< GIMLI::MeshEntity * >       stdVectorMeshEntity;
    typedef std::vector< GIMLI::CubicFunct >         stdVectorCubicFunct;
    typedef std::vector< GIMLI::RegionMarker >       stdVectorRegionMarker;
    typedef std::vector< GIMLI::BlockMatrixEntry >   stdVectorBlockMatrixEntry;

    typedef std::map< std::string, GIMLI::Vector< double > > stdMapStringRVector;
    typedef std::map< GIMLI::SIndex, GIMLI::Region * >         stdMapRegion;
    typedef std::map< float, float >                    stdMapF_F;
    typedef std::map< float, Complex >                  stdMapF_C;
    typedef std::map< int, double >                     stdMapI_D;
    typedef std::map< long, unsigned long >             stdMapL_UL;
    typedef std::map< unsigned long, double >           stdMapUL_D;
    typedef std::map< int, int >                        stdMapI_I;
    typedef std::map< std::string, std::string >        stdMapS_S;

#ifdef WIN32
    typedef std::map< long long, double > stdMapL_D;
    typedef std::map< unsigned long long, unsigned long long > stdMapL_L;
    typedef std::set< unsigned long long  > stdSetL;
    typedef std::map< std::pair< long long, long long > , double > stdMapL_L_D;
#else
    typedef std::map< long, double > stdMapL_D;
    typedef std::map<std::pair<unsigned long, unsigned long>, double > stdMapL_L_D;
    typedef std::set< long int >                        stdSetL;
    typedef std::map< long, long >                      stdMapL_L;
#endif

    typedef std::map<std::pair<unsigned long, unsigned long>, std::complex< double > > stdMapL_L_C;

    typedef std::map< GIMLI::Index, GIMLI::GraphDistInfo > NodeDistMap;

    // std::map<long long, double>
    // std::map<long long, long long>
    // std::map<unsigned long long, unsigned long long>
    // std::set<unsigned long long

    #ifdef _WIN64
        typedef std::map<std::pair<unsigned long long, unsigned long long>, double > stdMapLL_LL_D;
        typedef std::map<std::pair<unsigned long long, unsigned long long>, std::complex< double > > stdMapLL_LL_C;
    #elif _WIN32
        typedef std::map<std::pair<unsigned int, unsigned int>, double > stdMapI_I_D;
        typedef std::map<std::pair<unsigned int, unsigned int>, std::complex< double > > stdMapI_I_C;
    #endif

    typedef std::list< long unsigned int >              stdListUL;

    typedef std::set< std::string >                     stdSetS;
    typedef std::set< GIMLI::Node * >                 stdSetNodes;
    typedef std::set< GIMLI::Boundary * >             stdSetBoundaries;
    typedef std::set< GIMLI::Cell * >                 stdSetCells;

}} //pyplusplus::aliases

#endif // else not PYTEST

#endif // PYTHON_PYGIMLI__H
