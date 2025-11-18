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
#ifndef _GIMLI_GIMLI__H
#define _GIMLI_GIMLI__H

#ifdef HAVE_CONFIG_CMAKE_H
    #include <config.cmake.h>
#else
    #if defined(HAVE_CONFIG_H)
        #include <config.h>
        #define PACKAGE_AUTHORS "carsten@gimli.org, thomas@gimli.org, florian@gimli.org"
    #endif
#endif

#ifndef TRUE
    #define TRUE 1
#endif

#ifndef ON
    #define ON 1
#endif

#if BOOST_THREAD_FOUND || defined(HAVE_BOOST_THREAD_HPP)
    #define USE_BOOST_THREAD TRUE
#endif

#if BOOST_BIND_FOUND || defined(HAVE_BOOST_BIND_HPP)
    #define USE_BOOST_BIND TRUE
#endif

#ifndef PACKAGE_NAME
        #define PACKAGE_NAME "libgimli"
        #define PACKAGE_VERSION "untagt-win"

        #define PACKAGE_BUGREPORT "carsten@gimli.org"
        #define PACKAGE_AUTHORS "carsten@gimli.org, thomas@gimli.org, florian@gimli.org"
#endif // PACKAGE_NAME

#ifdef _MSC_VER
	#pragma warning(disable: 4251)
#endif

//#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cstdlib>
#include <cerrno>
#include <stdint.h>
#include <complex>
#include <algorithm>
#include <functional>

#include "platform.h"
#include "exitcodes.h"

//! GIMLi main namespace for the Geophyiscal Inversion and Modelling Library
namespace GIMLI{

#ifndef __USE_MISC
    typedef unsigned int uint;
#endif

typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

#ifdef _WIN64
    // for some magic reasons gccxml ignores the size_t = __int64 typedef under win64, so this helps let the python bindings through
//#warning WIN64
    typedef uint64 Index;
    typedef int64 SIndex;
#elif defined(_WIN32)
//#warning WIN32
    typedef uint32 Index;
    typedef int32 SIndex;
#elif defined(_MSC_VER)
//#warning WINMSC
	typedef unsigned __int32 Index;
    typedef __int32 SIndex;
#else
//    #warning DEFAULT
    #include <sys/types.h>
    typedef size_t Index;
    // typedef signed long SIndex;
    typedef ssize_t SIndex; // strange pygimli conversion into int on jenkins
#endif

#ifndef __ASSERT_FUNCTION
	#if defined(_MSC_VER)
		#define __ASSERT_FUNCTION __FUNCTION__
	#else
		#define __ASSERT_FUNCTION __PRETTY_FUNCTION__
	#endif
#endif


/*!Replace from with to inside str and return the result*/
DLLEXPORT  std::string replace(const std::string & str,
                               const std::string & from,
                               const std::string & to);

#ifndef SRC_DIR
    #define SRC_DIR "./"
#endif

#ifdef _WIN64__
    /// testme and refactor with below # 1.3.1
    #define __FILENAME__ GIMLI::replace(__FILE__, SRC_DIR, ".")
    // #define __FILENAME__ std::max<const char*>(__FILE__,\
    //     std::max(strrchr(__FILE__, '\\')+1, strrchr(__FILE__, '/')+1))
#else
    #define __FILENAME__ GIMLI::replace(__FILE__, SRC_DIR, ".")
#endif


#ifndef PYGIMLI_CAST // castxml complains on older gcc/clang
inline std::string str(){ return "";}
#endif
//! General template for conversion to string, should supersede all sprintf etc.
template< typename T > inline std::string str(const T & v){
    std::ostringstream os;
    os << v;
    return os.str();
}
enum LogType {Verbose, Info, Warning, Error, Debug, Critical};
DLLEXPORT void log(LogType type, const std::string & msg);

template<typename Value, typename... Values>
std::string str(Value v, Values... vs){
    std::ostringstream os;
    using expander = int[];
    os << v; // first
    (void) expander{ 0, (os << " " << vs, void(), 0)... };
    return os.str();
}

template<typename... Values>
void log(LogType type, Values... vs){
    return log(type, str(vs...));
}

// simple python like std out
inline void print() {
    std::cout << std::endl;
}

template<class Head>
void print(std::ostream& s, Head&& head) {
    s << std::forward<Head>(head) << std::endl;
}

template<class Head, class... Tail>
void print(std::ostream& s, Head&& head, Tail&&... tail) {
    s << std::forward<Head>(head) << " ";
    print(s, std::forward<Tail>(tail)...);
}

template<class... Args>
void print(Args&&... args) {
    print(std::cout, std::forward<Args>(args)...);
}

#define WHERE GIMLI::str(__FILENAME__) + ":" + GIMLI::str(__LINE__) + "\t"
#define WHERE_AM_I WHERE + "\t" + GIMLI::str(__ASSERT_FUNCTION) + " "
#define TO_IMPL WHERE_AM_I + " not yet implemented\n " + GIMLI::versionStr() + "\nPlease send the messages above, the commandline and all necessary data to the author."

#define __M std::cout << "*** " << WHERE << std::endl;
#define __MS(str) std::cout << "*** " <<str << " " << WHERE << std::endl;
// temporary
#define __MSP(...) GIMLI::print("+++", WHERE, "\n", __VA_ARGS__, "\n---");
#define __D if (debug()) std::cout << "Debug: " << WHERE << std::endl;
#define __DS(str) if (debug()) std::cout << "Debug: " << str << " " << WHERE << std::endl;

#define THROW_TO_IMPL throwToImplement(TO_IMPL);
#define CERR_TO_IMPL std::cerr << TO_IMPL << std::endl;
#define DEPRECATED std::cerr << WHERE_AM_I << " is deprecated " << std::endl;
#define DEPR_STR(s) std::cerr << WHERE_AM_I << " is deprecated. Use: " << s " instead."<< std::endl;
#define COUTMARKER std::cerr << WHERE_AM_I << std::endl;
#define UNTESTED std::cerr << "WARNING!" << WHERE_AM_I << " " << "this function is untested" << std::endl;

#define TOLERANCE 1e-12
#define TOUCH_TOLERANCE 1e-12
#define MAX_DOUBLE std::numeric_limits<double>::max()
#define MAX_INT std::numeric_limits< int >::max()

#define MESHBINSUFFIX   ".bms"
#define MATRIXBINSUFFIX ".bmat"
#define MATRIXASCSUFFIX ".matrix"
#define VECTORBINSUFFIX ".bvec"
#define VECTORASCSUFFIX ".vector"
#define NOT_DEFINED "notDefined"

#define ASSERT_EQUAL_SIZE(m, n) if (m.size() != n.size()) \
    throwLengthError(WHERE_AM_I + " " + str(m.size()) + " != " + str(n.size()));
#define ASSERT_THIS_SIZE(n) if (n < 0 || n >= this->size()) \
    throwLengthError(WHERE_AM_I + " " + str(this->size()) + " <= " + str(n));
#define ASSERT_VEC_SIZE(vec, n) if (n != vec.size()) \
    throwLengthError(WHERE_AM_I + " " + str(vec.size()) + " != " + str(n));
#define ASSERT_SIZE(vec, n) if (n < 0 || n >= vec.size()) \
    throwLengthError(WHERE_AM_I + " " + str(vec.size()) + " <= " + str(n));
#define ASSERT_EQUAL(m, n) if (m != n) \
    throwLengthError(WHERE_AM_I + " " + str(m) + " != " + str(n));
#define ASSERT_RANGE(i, start, end) if (i < start || i >= end) \
    throwRangeError(WHERE_AM_I, i, start, end);
#define ASSERT_EMPTY(v) if (v.size()==0) \
    throwLengthError(WHERE_AM_I + " array size is zero.");
#define ASSERT_PTR(v) if (v == 0) \
    throwError(WHERE_AM_I + " object not initialized.");

static const int MARKER_BOUND_HOMOGEN_NEUMANN = -1;
static const int MARKER_BOUND_MIXED = -2;
static const int MARKER_BOUND_HOMOGEN_DIRICHLET = -3;
static const int MARKER_BOUND_DIRICHLET = -4;
static const int MARKER_CELL_PARAMETER = 2;
static const int MARKER_NODE_SENSOR = -99;
static const int MARKER_FIXEDVALUE_REGION = -1000000;

static const uint8 MESH_BASEENTITY_RTTI      = 00;
static const uint8 MESH_MESHENTITY_RTTI      = 01;
static const uint8 MESH_NODE_RTTI            = 10;
static const uint8 MESH_BOUNDARY_RTTI        = 20;
static const uint8 MESH_BOUNDARY_NODE_RTTI   = 21;
static const uint8 MESH_EDGE_RTTI            = 22;
static const uint8 MESH_EDGE3_RTTI           = 23;
static const uint8 MESH_TRIANGLEFACE_RTTI    = 24;
static const uint8 MESH_TRIANGLEFACE6_RTTI   = 25;
static const uint8 MESH_QUADRANGLEFACE_RTTI  = 26;
static const uint8 MESH_QUADRANGLEFACE8_RTTI = 27;
static const uint8 MESH_POLYGON_FACE_RTTI    = 28;
static const uint8 MESH_CELL_RTTI            = 30;
static const uint8 MESH_EDGE_CELL_RTTI       = 31;
static const uint8 MESH_EDGE3_CELL_RTTI      = 32;
static const uint8 MESH_TRIANGLE_RTTI        = 33;
static const uint8 MESH_TRIANGLE6_RTTI       = 34;
static const uint8 MESH_QUADRANGLE_RTTI      = 35;
static const uint8 MESH_QUADRANGLE8_RTTI     = 36;
static const uint8 MESH_QUADRANGLE9_RTTI     = 37;
static const uint8 MESH_TETRAHEDRON_RTTI     = 41;
static const uint8 MESH_TETRAHEDRON10_RTTI   = 42;
static const uint8 MESH_HEXAHEDRON_RTTI      = 43;
static const uint8 MESH_HEXAHEDRON20_RTTI    = 44;
static const uint8 MESH_TRIPRISM_RTTI        = 45;
static const uint8 MESH_TRIPRISM15_RTTI      = 46;
static const uint8 MESH_PYRAMID_RTTI         = 47;
static const uint8 MESH_PYRAMID13_RTTI       = 48;

static const uint8 MESH_SHAPE_NODE_RTTI        = 210;
static const uint8 MESH_SHAPE_EDGE_RTTI        = 211;
static const uint8 MESH_SHAPE_TRIANGLE_RTTI    = 221;
static const uint8 MESH_SHAPE_QUADRANGLE_RTTI  = 222;
static const uint8 MESH_SHAPE_POLYGON_FACE_RTTI = 223;
static const uint8 MESH_SHAPE_TETRAHEDRON_RTTI = 231;
static const uint8 MESH_SHAPE_HEXAHEDRON_RTTI  = 232;
static const uint8 MESH_SHAPE_TRIPRISM_RTTI    = 233;
static const uint8 MESH_SHAPE_PYRAMID_RTTI     = 234;

static const uint8 GIMLI_MATRIXBASE_RTTI        = 0;
static const uint8 GIMLI_MATRIX_RTTI            = 1;
static const uint8 GIMLI_SPARSE_MAP_MATRIX_RTTI = 2;
static const uint8 GIMLI_SPARSE_CRS_MATRIX_RTTI = 3;
static const uint8 GIMLI_BLOCKMATRIX_RTTI       = 4;

/*! Flag load/save Ascii or binary */
enum IOFormat{Ascii, Binary};

//** start forward declaration
class Boundary;
class Cell;
class DataContainer;
class Line;
class MatrixBase;
class Mesh;
class MeshEntity;
class ModellingBase;
class Node;
class Plane;
class Region;
class RegionManager;
class Shape;
class Stopwatch;
class Pos;
class SolverWrapper;

typedef Pos RVector3;
typedef std::complex < double > Complex;

template < class ValueType >        class SparseMatrix;
typedef SparseMatrix< int >         ISparseMatrix;
typedef SparseMatrix< double >      RSparseMatrix;
typedef SparseMatrix< Complex >     CSparseMatrix;

template< class ValueType, class IndexType > class SparseMapMatrix;
typedef SparseMapMatrix< int, Index >     ISparseMapMatrix;
typedef SparseMapMatrix< double, Index >  RSparseMapMatrix;
typedef SparseMapMatrix< Complex, Index >  CSparseMapMatrix;

template < class ValueType > class Matrix;
template < class ValueType > class BlockMatrix;
template < class ValueType > class Matrix3;
template < class ValueType > class Vector;

typedef Vector < double > RVector;
typedef Vector < Complex > CVector;
typedef Vector < Pos > PosVector;
typedef PosVector R3Vector;
typedef Vector < bool > BVector;
typedef Vector < SIndex > IVector;
typedef Vector < Index > IndexArray;

typedef std::vector < SIndex > SIndexArray;

typedef Matrix < double > RMatrix;
typedef Matrix3< double > RMatrix3;
typedef Matrix < Complex > CMatrix;
typedef BlockMatrix < double > RBlockMatrix;


//#typedef Vector< unsigned char > BVector;

// typedef std::vector < RVector > RMatrix;
// typedef std::vector < CVector > CMatrix;

template < class ValueType > class PolynomialFunction;
typedef PolynomialFunction< double > RPolynomialFunction;

template < class ValueType > class ElementMatrix;
typedef ElementMatrix < double > RElementMatrix;

template < class Vec > class Trans;

/*! */
DLLEXPORT void savePythonGIL(bool s);
DLLEXPORT bool pythonGIL();

/*! Set global gimli debug flag on or off */
DLLEXPORT void setDebug(bool s);
DLLEXPORT bool debug();

/*! For several levels of deep debugging. Mainly used for python rvalue conversion. */
DLLEXPORT void setDeepDebug(int level);
DLLEXPORT int deepDebug();

/*! Set maximum amount of threads used by thirdparty software (e.g. openblas).
Default is number of CPU. */
DLLEXPORT void setThreadCount(Index nThreads);
DLLEXPORT Index threadCount();

DLLEXPORT void setUseOMP(bool o);
DLLEXPORT bool useOMP();

/*! For some debug purposes only */
DLLEXPORT void showSizes();

class DLLEXPORT PythonGILSave {
public:
    PythonGILSave() : saved_(false) { save(); }
    ~PythonGILSave() { restore(); }
    void save() ;
    void restore();
private:
    bool saved_;
#ifdef PYGIMLI
    PyThreadState *save_;
#endif
};

#define ALLOW_PYTHON_THREADS PythonGILSave __pygil_t__;
#define RESTORE_PYTHON_THREADS __pygil_t__.restore();
#define SAVE_PYTHON_THREADS __pygil_t__.save();

#ifndef PYTHON_FOUND
    #ifndef Py_BEGIN_ALLOW_THREADS
        #define Py_BEGIN_ALLOW_THREADS
        #define Py_END_ALLOW_THREADS
    #endif
#else
//     #include <Python.h>
#endif

DLLEXPORT std::string versionStr();

DLLEXPORT std::string authors();

template < class T, class U > T min(const T & a, const U & b){ return std::min(a, T(b)); }
template < class T, class U > T max(const T & a, const U & b){ return std::max(a, T(b)); }

DLLEXPORT int openFile(const std::string & fname, std::fstream * file,
                        std::ios_base::openmode farg, bool terminate);

inline int openInFileTerm(const std::string & fname, std::fstream * file){
    return openFile(fname, file, std::ios::in, true);
}
inline int openInFile(const std::string & fname, std::fstream * file,
                      bool terminate=true){
    return openFile(fname, file, std::ios::in, terminate);
}
inline int openOutFile(const std::string & fname, std::fstream * file,
                       bool terminate=true){
    return openFile(fname, file, std::ios::out, terminate);
}

DLLEXPORT bool fileExist(const std::string & filename);
DLLEXPORT uint countColumnsInFile(const std::string & fname, uint & columnCount);
DLLEXPORT uint countColumnsInFile(const std::string & fname);
DLLEXPORT uint countRowsInFile(const std::string & fname);
DLLEXPORT uint fileLength(std::fstream & file);
DLLEXPORT std::vector < std::string > getRowSubstrings(std::fstream & file, char comment='#');

inline std::vector < std::string > getRow(std::fstream & file, char comment='#'){
    return getRowSubstrings(file, comment);
}

DLLEXPORT std::vector < std::string > getNonEmptyRow(std::fstream & file, char comment='#');
DLLEXPORT std::vector < std::string > getCommentLine(std::fstream & file, char comment='#');

DLLEXPORT std::vector < std::string > getSubstrings(const std::string & str);
DLLEXPORT std::vector < std::string > split(const std::string & str, char delimiter=':');

DLLEXPORT std::map < float, Complex > loadCFloatMap(const std::string & filename);
DLLEXPORT std::map < float, float > loadFloatMap(const std::string & filename);
DLLEXPORT std::map < int, int > loadIntMap(const std::string & filename);

inline void convert(bool          & var, char * opt) { var = true; }
inline void convert(int           & var, char * opt) { if (!opt) var ++;    else var = atoi(opt); }
inline void convert(uint          & var, char * opt) { if (!opt) var ++;    else var = atoi(opt); }
inline void convert(Index         & var, char * opt) { if (!opt) var ++;    else var = atoi(opt); }
inline void convert(float         & var, char * opt) { if (!opt) var = 0.0f; else var = (float)atof(opt); }
inline void convert(double        & var, char * opt) { if (!opt) var = 0.0; else var = atof(opt); }
inline void convert(std::string   & var, char * opt) { if (!opt) var = "";  else var = opt ; }
inline void convert(std::vector < std::string >  & var, char * opt) { if (opt) var.push_back(opt); }

inline std::string type(const bool          & var) { return "bool"; }
inline std::string type(const int32         & var) { return "int32"; }
inline std::string type(const int64         & var) { return "int64"; }
inline std::string type(const uint32        & var) { return "uint32"; }
inline std::string type(const uint64        & var) { return "uint64"; }
// inline std::string type(const Index         & var) { return "Index"; }
// inline std::string type(const SIndex        & var) { return "SIndex"; }
inline std::string type(const float         & var) { return "float"; }
inline std::string type(const double        & var) { return "double"; }
inline std::string type(const Complex       & var) { return "complex"; }
inline std::string type(const std::string   & var) { return "string"; }
inline std::string type(const std::vector < std::string >  & var) { return "string"; }

inline std::string type(const RVector  & var) { return "RVector"; }
inline std::string type(const RVector3 & var) { return "RVector3"; }
inline std::string type(const R3Vector & var) { return "R3Vector"; }
inline std::string type(const CVector  & var) { return "CVector"; }
inline std::string type(const RMatrix  & var) { return "RMatrix"; }
inline std::string type(const CMatrix  & var) { return "CMatrix"; }

inline int       toInt(const std::string & str){ return std::atoi(str.c_str()); }
inline float   toFloat(const std::string & str){ return (float)std::atof(str.c_str()); }
inline double toDouble(const std::string & str){ return std::strtod(str.c_str(), NULL); }

/*! Read value from environment variable.
 * Return default value if environment not set.
 * Environment var can be set in sh via:
 * export name=val, or simple passing name=val in front of executable. */
template < typename ValueType > ValueType getEnvironment(const std::string & name,
                                                         ValueType def,
                                                         bool verbose=false){
    ValueType var = def;

    char * cVar = getenv(name.c_str());
    if (cVar != NULL){
        convert(var, cVar);
        if (verbose) std::cout << "Found: export " << name << "=" << cVar << std::endl;
    }
    return var;
}

/*! Set environment variable. Probably only for internal use and maybe only
 * for posix systems*/
template < typename ValueType > void setEnvironment(const std::string & name,
                                                    ValueType val,
                                                    bool verbose=false){
    int ret = setenv(name.c_str(), str(val).c_str(), 1);
    switch(ret){
        case EINVAL:
            __MS(name << " " << val)
            throwError("name is NULL, points to a string of length 0, or contains an '=' character.");
        case ENOMEM:
            __MS(name << " " << val)
            throwError("name is NULL, Insufficient memory to add a new variable to the environment.");
    }
    if (verbose) std::cout << "set: export " << name << "=" << val << std::endl;
}

/*!Replace from with to inside str and return the result*/
DLLEXPORT std::string replace(const std::string & str, const char from, const char to);

/*!convert all chars in str to lower and return the result*/
DLLEXPORT std::string lower(const std::string & str);

// template < typename T > inline void swapVal(T & a, T & m){
//     T tmp(a); a = m; m = tmp;
// }

/*! General template for deleting an object. This is not exception-safe unless you use some kind of smart pointer.\n
Example: Delete all objects in a container.
vector < ptr * > vP;
// ... // something that fills vP with the new operator.
for_each(vP.begin(), vP.end(), deletePtr());
*/
struct DLLEXPORT deletePtr{
    template < typename T > void operator()(T * p) { delete p; }
};
struct DLLEXPORT cerrPtr{
  template < typename T > void operator() (const T * p) const { std::cerr << p << " " << std::endl; }
};
struct DLLEXPORT cerrPtrObject{
  template < typename T > void operator() (const T * p) const { std::cerr << *p << " " << std::endl; }
};
struct DLLEXPORT coutPtr{
  template < typename T > void operator() (const T * p) const { std::cout << p << " " << std::endl; }
};
struct DLLEXPORT coutPtrObject{
  template < typename T > void operator() (const T * p) const { std::cout << *p << " " << std::endl; }
};

template < typename Set > inline void intersectionSet(Set & dest,
                                                      const Set & a,
                                                      const Set & b){
    dest.clear();
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(dest, dest.begin()));
}
template < typename Set > inline void intersectionSet(Set & dest,
                                                      const Set & a,
                                                      const Set & b,
                                                      const Set & c){
    dest.clear();
    set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                     std::inserter(dest, dest.begin()));
    Set tmp(dest);
    dest.clear();
    set_intersection(tmp.begin(), tmp.end(), c.begin(), c.end(),
                     std::inserter(dest, dest.begin()));
}

template < typename Set > inline void intersectionSet(Set & dest,
                                                      const Set & a,
                                                      const Set & b,
                                                      const Set & c,
                                                      const Set & d){
    dest.clear();
    set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                     std::inserter(dest, dest.begin()));
    Set tmp(dest);
    dest.clear();
    set_intersection(tmp.begin(), tmp.end(), c.begin(), c.end(),
                     std::inserter(dest, dest.begin()));
    tmp = dest;
    dest.clear();
    set_intersection(tmp.begin(), tmp.end(), d.begin(), d.end(),
                     std::inserter(dest, dest.begin()));
}

template < typename Set > inline void intersectionSet(Set & dest,
                                                      const std::vector < Set > & a){
    if (a.size() > 1) {
        intersectionSet(dest, a[0], a[1]);
        for (size_t i = 2; i < a.size(); i ++){
            Set tmp(dest);
            dest.clear();
            set_intersection(tmp.begin(), tmp.end(), a[i].begin(), a[i].end(),
                             std::inserter(dest, dest.begin()));
        }
    } else if (a.size() == 1){
        dest = a[0];
    } else {
        dest.clear();
    }
}

template < class T > class IncrementSequence {
public:
    IncrementSequence(T initialValue = 0) : value_(initialValue) {   }
    inline T operator() () { return value_++; }
private:
    T value_;
};

/*! Template class for singleton instances */
template < typename Classname > class DLLEXPORT Singleton {
public:

    virtual ~Singleton() {
    #ifndef PYGIMLI_CAST
        delete pInstance_; pInstance_ = NULL;
    #endif
    }

    /*! This call create one instance of the class and return a pointer to it. */
    static Classname * pInstance() {
    #ifndef PYGIMLI_CAST
        return pInstance_ ? pInstance_ : (pInstance_ = new Classname());
    #endif
        return 0;
    }

    /*! This call create one instance of the class and return a reference to it. */
    static Classname & instance() {
        return *pInstance();
    }

protected:
    /*! Protected so it can only be called from derived classes */
    Singleton(){ }
    Singleton(const Singleton &){__M};
private:
    /*! Private so that it can not be called */

    /*! Copy constructor is private, so don't use it */
#ifndef PYGIMLI_CAST
    static Classname * pInstance_;
#endif
};

template < typename T > Index hash_(T v){
    #ifndef PYGIMLI_CAST
        return std::hash<T>()(v);
    #endif
    __M
    return 0;
}
/*! Combine
https://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine */
template <typename T>
void hashCombine(Index & seed, const T& val){
    seed ^= hash_(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}
template <typename T, typename... Types>
void hashCombine (Index & seed, const T & val, const Types&... args){
    hashCombine(seed, val);
    hashCombine(seed, args...);
}
inline void hashCombine (Index & seed){}
template <typename... Types>
Index hash(const Types&... args){
    Index seed = 0;
    hashCombine(seed, args...);
    return seed;
}
template void hashCombine(Index & seed, const Index & hash);
// template void hashCombine(Index & seed, const PosVector & val);
// template void hashCombine(Index & seed, const Pos & val);
// template void hashCombine(Index & seed, const DataContainer & val);

} // namespace GIMLI

#endif // _GIMLI_GIMLI__H
