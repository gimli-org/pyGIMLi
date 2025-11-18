#include <typeinfo>

#include "boost/python/object.hpp"  //len function
#include "boost/python/ssize_t.hpp" //ssize_t type definition
#include "boost/python/detail/none.hpp"
#include <boost/mpl/int.hpp>
#include <boost/mpl/next.hpp>
#include "tuples.hpp"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include <numpy/npy_math.h>

#include <Python.h>

#include "gimli.h"
#include "pos.h"
#include "vector.h"
#include "matrix.h"

// cp ../gimli/core/python/custom_rvalue.cpp core/python/generated/custom_rvalue.cpp && make pg
static int __numpy_initialized = 0;

// this fails for win since py39
// const static int __numpy_initialized = initNumpy();

int initNumpy(){
    // Needed or py* checks will segfault
	if (__numpy_initialized == 0){
        import_array2("Cannot import numpy.core.multiarray c-api for rvalue converters.", 0);
        __numpy_initialized = 1;
    }
    return 0;
}

namespace bp = boost::python;

#if defined ( __APPLE__ ) || ( defined (__SVR4) && defined (__sun) )
    #define __DC(str) ;
#else
    #define __DC(str) if (GIMLI::deepDebug() > 0) __MS(str)
#endif

// ** TODO check if we need a delete somewhere for all the new stuff

#include <fstream>
std::ostream & operator << (std::ostream & os, const bp::object& o){
    return os << bp::extract<std::string>(bp::str(o))();
}
namespace r_values_impl{

template < class ValueType > void * checkConvertibleSequenz(PyObject * obj){
    initNumpy();

    if (!obj){
        __DC("\t" << obj << "\t abborting .. !Object")
        return NULL;
    }
    __DC(obj << "(" << obj->ob_type->tp_name << ") -> sequenz of type ("
             << GIMLI::type(ValueType(0)) << ")" << typeid(ValueType).name())
         // FW: Caused problems during Mac build // still?

    // is obj is a sequence
    if(!PySequence_Check(obj)){
        __DC("\t" << obj << "\t abborting no sequence")
        return NULL;
    }

    // has the obj a len method
    if (!PyObject_HasAttrString(obj, "__len__")){
        __DC("\t" << obj << "\t abborting no len")
        return NULL;
    }

    if (PyObject_TypeCheck(obj, &PyArray_Type)){
        PyArrayObject *arr = (PyArrayObject *)obj;

        __DC(obj << "\t numpy.ndarray to " << typeid(ValueType).name() << " "
        << typeid(GIMLI::Index).name() << " " << typeid(float).name()<< "... okay" << " ndim: " << PyArray_NDIM(arr))

        if (typeid(ValueType) == typeid(GIMLI::Pos)){
            if (PyArray_NDIM(arr) == 1){
                __DC("\t" << obj << "\t abborting .. numpy array.ndim == 1")
                return NULL;
            }
        }
        if (typeid(ValueType) == typeid(GIMLI::RMatrix)){
            if (PyArray_NDIM(arr) != 2){
                __DC("\t" << obj << "\t ndarray.ndim != 2 and is non convertable to std::vector(GIMLI::RMatrix)")
                return NULL;
            }
        }

        if (typeid(ValueType) == typeid(GIMLI::Index)){
            if (PyArray_TYPE(arr) == NPY_BOOL){
                __DC("\t" << obj << "\t Object is nd.array with dtype == bool*"
                     "non convertable to GIMLI::IVector")
                return NULL;
            }
        } else if (typeid(ValueType) == typeid(bool)){
           if (PyArray_TYPE(arr) != NPY_BOOL){
                __DC("\t" << obj << "\t Object is nd.array with dtype != bool"
                     ".. non convertable to GIMLI::BVector")
                return NULL;
            }
        }

        return obj;
    }

    bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

    if (len(py_sequence) > 0) {
        // FW: Causes problems on Mac build.
        // __DC(obj << "\t len: " << len(py_sequence) << " type: " << GIMLI::type(ValueType(0)) << ": " << typeid(ValueType).name())
        bp::object element = py_sequence[0];
        __DC(obj << "\t seq[0]: " << element << " is of type: " << element.ptr()->ob_type->tp_name)

        //** do not convert [long] || [ulong]  > [bool]
        if (typeid(ValueType) == typeid(GIMLI::Index) || typeid(ValueType) == typeid(GIMLI::SIndex)){
            if (strcmp(element.ptr()->ob_type->tp_name, "bool") == 0) {
                __DC(obj << "\t abborting: Index requested but sequence of "<< element.ptr()->ob_type->tp_name)
                return NULL;
            }
        } else if (typeid(ValueType) == typeid(bool)){
            // special check for BVector .. we oonly want to convert sequence of bool objects
            if (strcmp(element.ptr()->ob_type->tp_name, "bool") != 0) {
                __DC(obj << "\t abborting: bools requested but sequence of "<< element.ptr()->ob_type->tp_name)
                return NULL;
            }
        }

        bp::extract< ValueType > type_checker(element);
        if (type_checker.check()){
            __DC(obj << "\t ->construct: " << len(py_sequence))
            return obj;
        } else {
            __DC(obj << "\t cannot convert: " << type_checker.check())

            PyTypeObject* type = obj->ob_type;
            const char* p = type->tp_name;
            __DC("type is " << p)

            //use PyObject_IsInstance(obj, ndarray())
            if (strcmp(p, "numpy.ndarray")==0){
                return obj;
            }

            __DC(WHERE_AM_I << "element cannot converted ")
        }

    } else {
        __DC(obj << " len == 0")
        return NULL;
    }

    __DC(obj << " fail")
    return NULL;
}

template < class ValueType > void * checkConvertibleNumpyScalar(PyObject * obj){
    initNumpy();

    // will be used for every convert of numpy scalars here e.g. for list conversion
    if (!obj){
        __DC("\t" << obj << "\t abort check .. !Object")
        return NULL;
    }
    if (GIMLI::deepDebug() > 0){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> " +
            GIMLI::type(ValueType(0))) // FW: Caused problems during Mac build
        __DC("\tType:" << Py_TYPE(obj))
        __DC("\tArray:" << PyObject_TypeCheck(obj, &PyArray_Type))
        __DC("\tPyGenericArrType_Type:" << PyObject_TypeCheck(obj, &PyGenericArrType_Type))
        __DC("\tPyIntegerArrType_Type:" << PyObject_TypeCheck(obj, &PyIntegerArrType_Type))
        __DC("\tPySignedIntegerArrType_Type:" << PyObject_TypeCheck(obj, &PySignedIntegerArrType_Type))
        __DC("\tPyUnsignedIntegerArrType_Type:" << PyObject_TypeCheck(obj, &PyUnsignedIntegerArrType_Type))
        __DC("\tPyIntArrType_Type:" << PyObject_TypeCheck(obj, &PyIntArrType_Type))
        __DC("\tPyLongArrType_Type:" << PyObject_TypeCheck(obj, &PyLongArrType_Type))
        __DC("\tPyUIntArrType_Type:" << PyObject_TypeCheck(obj, &PyUIntArrType_Type))
        __DC("\tPyULongArrType_Type:" << PyObject_TypeCheck(obj, &PyULongArrType_Type))
        __DC("\tPyFloatArrType_Type:" << PyObject_TypeCheck(obj, &PyFloatArrType_Type))
        __DC("\tPyDoubleArrType_Type:" << PyObject_TypeCheck(obj, &PyDoubleArrType_Type))
    }

    if (PyObject_TypeCheck(obj, &PyGenericArrType_Type)){
        if (typeid(ValueType) == typeid(GIMLI::Index)){
            if (!PyObject_TypeCheck(obj, &PyIntegerArrType_Type)){
                __DC("\t" << obj << "\t abort check .. Object cannot convert to GIMLI::Index")
                return NULL;
            }
        }
        __DC("\t" << obj << "\t ok.")
        return obj;
    }
    __DC("\t" << obj << "\t abort: no numpy scalar.")
    return NULL;
}


struct PyTuple2RVector3{

    typedef boost::tuples::tuple< double, double > xy_type;
    typedef boost::tuples::tuple< double, double, double > xyz_type;
    typedef bp::from_py_sequence< xy_type > xy_converter_type;
    typedef bp::from_py_sequence< xyz_type > xyz_converter_type;

    typedef GIMLI::RVector3 xyz_t;

    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> RVector3")
        if (xy_converter_type::convertible(obj) ||
             xyz_converter_type::convertible(obj)){
            return obj;
        } else{
            return NULL;
        }
    }

    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){

        typedef bp::converter::rvalue_from_python_storage< xyz_t > xyz_storage_t;

        xyz_storage_t * the_storage = reinterpret_cast< xyz_storage_t * >(data);
        void * memory_chunk = the_storage->storage.bytes;

        double x(0.0), y(0.0), z(0.0);

        bp::tuple py_tuple(bp::handle<>(bp::borrowed(obj)));

        if (3 == bp::len(py_tuple)){
            boost::tuples::tie(x, y, z) = xyz_converter_type::to_c_tuple(obj);
        } else if (2 == bp::len(py_tuple)){
            boost::tuples::tie(x, y) = xy_converter_type::to_c_tuple(obj);
        }

        //** don't know where this will be deleted but it is necessary
        new (memory_chunk) xyz_t(x, y, z);
        data->convertible = memory_chunk;
    }
};

struct PySequence2RVector{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> RVector")
        return checkConvertibleSequenz< double >(obj);
    }

    /*! Convert List[] or ndarray into RVector */
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
       __DC(obj << " constructing RVector:")

        typedef bp::converter::rvalue_from_python_storage< GIMLI::Vector< double > > storage_t;
        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));
        GIMLI::Vector< double > * vec = new (memory_chunk) GIMLI::Vector< double >(len(py_sequence));
        data->convertible = memory_chunk;

        if (strcmp(obj->ob_type->tp_name, "numpy.ndarray") == 0){
            PyArrayObject *arr = (PyArrayObject *)obj;
            __DC("type is " << obj->ob_type->tp_name << " " << PyArray_TYPE(arr))

            if (PyArray_TYPE(arr) == 12 && PyArray_ISONESEGMENT(arr)){
                void * arrData = PyArray_DATA(arr);
                std::memcpy(&(*vec)[0], arrData, vec->size() * sizeof(double));
                return;
            } else if (PyArray_TYPE(arr) == 7 && PyArray_ISONESEGMENT(arr)){ //numpy.int64
                    __DC(arr << " ** from array of type " << PyArray_TYPE(arr))

                    bp::object element;

                    for (GIMLI::Index i = 0; i < vec->size(); i ++){
                        element = py_sequence[i];
                        (*vec)[i] = PyArrayScalar_VAL(element.ptr(), Int64);

    //                     __DC(i << " a " << element);
    //                     __DC(i << " a " << element.ptr()->ob_type->tp_name)
    //                     __DC(i << " d " <<  PyArrayScalar_VAL(element.ptr(), Int64));

                    }
                    return;
            } else {
                    __DC("fixme: type=" << PyArray_TYPE(arr))
            }
        }

        // convert from list
        __DC(obj << " ** from sequence ")

        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = bp::extract< double >(py_sequence[i]);
        }
    }
private:
};

struct PySequence2CVector{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> CVector")
        return checkConvertibleSequenz<GIMLI::Complex>(obj);
    }

    /*! Convert List[] or ndarray into RVector */
    static void construct(PyObject* obj,
                          bp::converter::rvalue_from_python_stage1_data * data){
       __DC(obj << " constructing CVector:")

        typedef bp::converter::rvalue_from_python_storage< GIMLI::Vector< GIMLI::Complex > > storage_t;
        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));
        GIMLI::Vector< GIMLI::Complex > * vec =
            new (memory_chunk) GIMLI::Vector< GIMLI::Complex >(len(py_sequence));
        data->convertible = memory_chunk;

        if (strcmp(obj->ob_type->tp_name, "numpy.ndarray") == 0){
            PyArrayObject *arr = (PyArrayObject *)obj;
            __DC("type is " << obj->ob_type->tp_name << " " << PyArray_TYPE(arr))

            if (PyArray_TYPE(arr) == 12 && PyArray_ISONESEGMENT(arr)){
                void * arrData = PyArray_DATA(arr);
                std::memcpy(&(*vec)[0], arrData, vec->size() * sizeof(double));
                return;
            } else if (PyArray_TYPE(arr) == 7 && PyArray_ISONESEGMENT(arr)){ //numpy.int64
                    __DC(arr << " ** from array of type " << PyArray_TYPE(arr))

                    bp::object element;

                    for (GIMLI::Index i = 0; i < vec->size(); i ++){
                        element = py_sequence[i];
                        (*vec)[i] = PyArrayScalar_VAL(element.ptr(), Int64);
    //                     __DC(i << " a " << element);
    //                     __DC(i << " a " << element.ptr()->ob_type->tp_name)
    //                     __DC(i << " d " <<  PyArrayScalar_VAL(element.ptr(), Int64));
                    }
                    return;
            } else if (PyArray_TYPE(arr) == 15 && PyArray_ISONESEGMENT(arr)){ //numpy.complex
                __DC(arr << " ** from array of type " << PyArray_TYPE(arr))

                bp::object element;

                for (GIMLI::Index i = 0; i < vec->size(); i ++){
                    element = py_sequence[i];

#ifndef NPY_2_0_API_VERSION
                    (*vec)[i] = GIMLI::Complex(
                            PyArrayScalar_VAL(element.ptr(), Complex128).real,
                            PyArrayScalar_VAL(element.ptr(), Complex128).imag);
#else
                    //see and fix me! https://numpy.org/devdocs/release/2.0.0-notes.html#complex-types-underlying-c-type-changes
                    (*vec)[i] = GIMLI::Complex(
                            npy_creal(PyArrayScalar_VAL(element.ptr(), Complex128)),
                            npy_cimag(PyArrayScalar_VAL(element.ptr(), Complex128)));
#endif
//                     (*vec)[i] = PyArrayScalar_VAL(element.ptr(),
//                     __DC(i << " a " << element);
//                     __DC(i << " a "
//                         << bp::extract< double >(element.attr('real')));
//                     __DC(i << " a "
//                     << bp::extract< double >(bp::extract<bp::tuple>(element)));
//                     __DC(i << " a " << element.ptr()->ob_type->tp_name);
//                     PyArrayScalar_VAL(element.ptr(), Complex128).real;
//                     PyArrayScalar_VAL(element.ptr(), Complex128).imag;
//                     __DC(i << " d " <<  PyArrayScalar_VAL(element.ptr(),
//                                                           Complex128));
                }
                return;
            } else {
                 __DC("fixme: type="
                 << PyArray_TYPE(arr) << " " << PyArray_ISONESEGMENT(arr))
            }
        }

        // convert from list
        __DC(obj << " ** from sequence not implemented")
        GIMLI::throwToImplement("implementme: PySequence2CVector");
//         for (GIMLI::Index i = 0; i < vec->size(); i ++){
//             (*vec)[i] = bp::extract< double >(py_sequence[i]);
//         }
     }
private:
};

struct PySequence2IndexArray{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> IndexArray")
        return checkConvertibleSequenz<GIMLI::Index>(obj);
    }

    /*! Convert obj into IndexArray */
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        __DC(obj << "\t constructing IndexArray")
        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

        typedef bp::converter::rvalue_from_python_storage< GIMLI::IndexArray > storage_t;

        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        GIMLI::IndexArray * vec = new (memory_chunk) GIMLI::IndexArray(len(py_sequence));
        data->convertible = memory_chunk;
        __DC(obj << "\t from list")
        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = bp::extract< GIMLI::Index >(py_sequence[i]);
        }
    }
private:
};

struct PySequence2IVector{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> IVector")
        return checkConvertibleSequenz<GIMLI::SIndex>(obj);
    }

    /*! Convert obj into IndexArray */
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        __DC(obj << "\t constructing IVector")
        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

        typedef bp::converter::rvalue_from_python_storage< GIMLI::IVector> storage_t;

        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        GIMLI::IVector * vec = new (memory_chunk) GIMLI::IVector(len(py_sequence));
        data->convertible = memory_chunk;
        __DC(obj << "\t from list")
        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            //__DC(obj << " " << i << " " << bp::extract< long >(py_sequence[i]))
            (*vec)[i] = bp::extract< GIMLI::SIndex >(py_sequence[i]);
        }
    }
private:
};

struct PySequence2BVector{
    static void * convertible(PyObject * obj){
        __DC(obj  << "(" << obj->ob_type->tp_name << ") -> BVector")
        return checkConvertibleSequenz< bool >(obj);
    }
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        __DC(obj << "\t constructing BVector")
        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

        typedef bp::converter::rvalue_from_python_storage< GIMLI::BVector> storage_t;

        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        GIMLI::BVector * vec = new (memory_chunk) GIMLI::BVector(len(py_sequence));
        data->convertible = memory_chunk;
        __DC(obj << "\t from list")
        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = PyArrayScalar_VAL(bp::object(py_sequence[i]).ptr(), Bool);

            // __DC(i << " " << bp::object(py_sequence[i]) << " " << (*vec)[i])
            //(*vec)[i] = bp::extract< bool >(py_sequence[i]);
        }
    }
};

struct PySequence2StdVectorRVector3{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
         __DC(obj << "(" << obj->ob_type->tp_name << ") -> StdVectorRVector3")
        return checkConvertibleSequenz< GIMLI::Pos >(obj);
    }

    /*! Convert obj into RVector */
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        __DC(obj << "\t constructing RVector3")
        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

        typedef bp::converter::rvalue_from_python_storage< std::vector < GIMLI::Pos > > storage_t;

        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        std::vector < GIMLI::Pos > *vec = new (memory_chunk) std::vector < GIMLI::Pos >(len(py_sequence));
        data->convertible = memory_chunk;

        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = bp::extract< GIMLI::Pos >(py_sequence[i]);
        }
    }
private:
};

struct PySequence2R3Vector{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> R3Vector")
        return checkConvertibleSequenz< GIMLI::Pos >(obj);
    }

    /*! Convert obj into RVector */
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        __DC(obj << "\t constructing R3Vector")
        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

        typedef bp::converter::rvalue_from_python_storage< GIMLI::R3Vector > storage_t;

        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        GIMLI::R3Vector *vec = new (memory_chunk) GIMLI::R3Vector(len(py_sequence));
        data->convertible = memory_chunk;

        for (GIMLI::Index i = 0; i < vec->size(); i ++){
            (*vec)[i] = bp::extract< GIMLI::Pos >(py_sequence[i]);
        }
    }
private:
};

struct Numpy2RMatrix{

    /*! Check if the object is convertible */
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> RMatrix")
        return checkConvertibleSequenz< GIMLI::Matrix< double > >(obj);
    }

    /*! Convert obj into RVector */
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        __DC(obj << "\t constructing RMatrix")

        bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

        typedef bp::converter::rvalue_from_python_storage< GIMLI::Matrix < double > > storage_t;

        storage_t* the_storage = reinterpret_cast<storage_t*>(data);
        void* memory_chunk = the_storage->storage.bytes;

        PyArrayObject *arr = (PyArrayObject *)obj;

        if (PyArray_TYPE(arr) == 12) {
            __DC("\ttype=" << PyArray_TYPE(arr)
                << " ISONESEGMENT:" << PyArray_ISONESEGMENT(arr)
                << " IS_C_CONTIGUOUS:" << PyArray_IS_C_CONTIGUOUS(arr)
                << " IS_F_CONTIGUOUS:" << PyArray_IS_F_CONTIGUOUS(arr)
                )
            int nDim = PyArray_NDIM(arr);
            if (nDim != 2){
                __DC("nDim=" << nDim)
                GIMLI::throwToImplement("Only numpy.ndarray with ndim == 2 can be converted to GIMLI::RMatrix");
            }
            GIMLI::Index rows = PyArray_DIM(arr, 0);
            GIMLI::Index cols = PyArray_DIM(arr, 1);
            GIMLI::Matrix < double > *mat = new (memory_chunk)
                                        GIMLI::Matrix < double >(rows, cols);
            data->convertible = memory_chunk;
            __DC("rows=" << rows << " cols=" << cols)

            if (PyArray_ISONESEGMENT(arr)){
                double * arrData = (double*)PyArray_DATA(arr);
                if (PyArray_IS_C_CONTIGUOUS(arr)){
                    for (GIMLI::Index i = 0; i < mat->rows(); i ++ ){
                        std::memcpy(&(*mat)[i][0],
                                    arrData + (i * mat->cols()),
                                    mat->cols() * sizeof(double));
                    }
                } else {
                    // assume Fortran like column orientated mem
                    // slow elementwise copy needed until someone knows
                    // a better way
                    for (GIMLI::Index i = 0; i < cols; i ++ ){
                        for (GIMLI::Index j = 0; j < rows; j ++ ){
                            mat->setVal(j, i, (double)arrData[j + i*rows]);
                        }
                    }
                }
            } else {
                GIMLI::throwToImplement("numpy.ndarray is not one segment .. not yet implemented.");
            }

            return;
        } else {
            __DC("implementme: type=" << PyArray_TYPE(arr))
        }

        GIMLI::throwToImplement("Unknown rvalue type conversion from numpy.ndarray of type " + GIMLI::str(PyArray_TYPE(arr)) +
         + " to GIMLI::RMatrix");

//         for (GIMLI::Index i = 0; i < vec->size(); i ++){
//             //(*mat)[i] = bp::extract< GIMLI::Pos >(py_sequence[i]);
//         }
    }
private:
};

template < class ValueType > void convertFromNumpyScalar(PyObject* obj,
                        bp::converter::rvalue_from_python_stage1_data * data){
    // __DC(obj << "\tNumpyScalar -> " + GIMLI::type(ValueType(0)) + " check OK: ") // FW: Mac problems
    //bp::object py_sequence(bp::handle<>(bp::borrowed(obj)));

    typedef bp::converter::rvalue_from_python_storage< ValueType > storage_t;

    storage_t* the_storage = reinterpret_cast<storage_t*>(data);
    void* memory_chunk = the_storage->storage.bytes;

    ValueType * val = new (memory_chunk) ValueType[1];
    data->convertible = memory_chunk;

    if (PyObject_TypeCheck(obj, &PyLongArrType_Type)){
        *val = PyArrayScalar_VAL(obj, Int32);
        __DC(obj << "\tnumpy.int32 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyLongLongArrType_Type)){
        *val = PyArrayScalar_VAL(obj, Int64);
        __DC(obj << "\tnumpy.int64 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyULongArrType_Type)){
        *val = PyArrayScalar_VAL(obj, UInt64);
        __DC(obj << "\tnumpy.uint32 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyULongLongArrType_Type)){
        *val = PyArrayScalar_VAL(obj, UInt64);
        __DC(obj << "\tnumpy.uint64 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyIntArrType_Type)){
        *val = PyArrayScalar_VAL(obj, Int32);
        __DC(obj << "\tnumpy.int32 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyUIntArrType_Type)){
        *val = PyArrayScalar_VAL(obj, UInt32);
        __DC(obj << "\tnumpy.uint32 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyFloatArrType_Type)){
        *val = PyArrayScalar_VAL(obj, Float32);
        __DC(obj << "\tnumpy.float32 = " << *val)
    } else if (PyObject_TypeCheck(obj, &PyDoubleArrType_Type)){
        *val = PyArrayScalar_VAL(obj, Float64);
        __DC(obj << "\tnumpy.float64 = " << *val)
    } else {
        __MS(obj << "\tconvertFromNumpyScalar -> unhandled dtype")
        __MS(obj << "\tconvertFromNumpyScalar -> name: " << obj->ob_type->tp_name)
        __MS("\tPyGenericArrType_Type:" << PyObject_TypeCheck(obj, &PyGenericArrType_Type))
        __MS("\tPyIntegerArrType_Type:" << PyObject_TypeCheck(obj, &PyIntegerArrType_Type))
        __MS("\tPySignedIntegerArrType_Type:" << PyObject_TypeCheck(obj, &PySignedIntegerArrType_Type))
        __MS("\tPyUnsignedIntegerArrType_Type:" << PyObject_TypeCheck(obj, &PyUnsignedIntegerArrType_Type))
        __MS("\tPyIntArrType_Type:" << PyObject_TypeCheck(obj, &PyIntArrType_Type))
        __MS("\tPyLongArrType_Type:" << PyObject_TypeCheck(obj, &PyLongArrType_Type))
        __MS("\tPyUIntArrType_Type:" << PyObject_TypeCheck(obj, &PyUIntArrType_Type))
        __MS("\tPyULongArrType_Type:" << PyObject_TypeCheck(obj, &PyULongArrType_Type))
        __MS("\tPyFloatArrType_Type:" << PyObject_TypeCheck(obj, &PyFloatArrType_Type))
        __MS("\tPyDoubleArrType_Type:" << PyObject_TypeCheck(obj, &PyDoubleArrType_Type))
    }
}


//template <typename T, NPY_TYPES NumPyScalarType>
struct Numpy2Long{
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> SIndex")
        return checkConvertibleNumpyScalar< GIMLI::SIndex >(obj);
    }
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        return convertFromNumpyScalar< GIMLI::SIndex >(obj, data);
    }
private:
};

//template <typename T, NPY_TYPES NumPyScalarType>
struct Numpy2ULong{
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> Index(Numpy2ULong)")
        return checkConvertibleNumpyScalar< GIMLI::Index >(obj);
    }
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        return convertFromNumpyScalar< GIMLI::Index >(obj, data);
    }
private:
};

//template <typename T, NPY_TYPES NumPyScalarType>
struct Numpy2Int{
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> int(Numpy2Int)")
        return checkConvertibleNumpyScalar< GIMLI::int32 >(obj);
    }
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        return convertFromNumpyScalar< GIMLI::int32 >(obj, data);
    }
private:
};

//template <typename T, NPY_TYPES NumPyScalarType>
struct Numpy2UInt{
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> uint(Numpy2UInt)")
        return checkConvertibleNumpyScalar< GIMLI::uint32 >(obj);
    }
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        return convertFromNumpyScalar< GIMLI::uint32 >(obj, data);
    }
private:
};

//template <typename T, NPY_TYPES NumPyScalarType>
struct Numpy2Double{
    static void * convertible(PyObject * obj){
        __DC(obj << "(" << obj->ob_type->tp_name << ") -> double(Numpy2Double)")
        return checkConvertibleNumpyScalar< double >(obj);
    }
    static void construct(PyObject* obj, bp::converter::rvalue_from_python_stage1_data * data){
        return convertFromNumpyScalar< double >(obj, data);
    }
private:
};

} //r_values_impl

void register_numpy_to_int64_conversion(){
    bp::converter::registry::push_back(& r_values_impl::Numpy2Long::convertible,
                                        & r_values_impl::Numpy2Long::construct,
                                        bp::type_id< GIMLI::SIndex >());
}
void register_numpy_to_uint64_conversion(){
    bp::converter::registry::push_back(& r_values_impl::Numpy2ULong::convertible,
                                        & r_values_impl::Numpy2ULong::construct,
                                        bp::type_id< GIMLI::Index >());
}
void register_numpy_to_int32_conversion(){
    bp::converter::registry::push_back(& r_values_impl::Numpy2Int::convertible,
                                        & r_values_impl::Numpy2Int::construct,
                                        bp::type_id< GIMLI::int32 >());
}
void register_numpy_to_uint32_conversion(){
    bp::converter::registry::push_back(& r_values_impl::Numpy2UInt::convertible,
                                        & r_values_impl::Numpy2UInt::construct,
                                        bp::type_id< GIMLI::uint32 >());
}
void register_numpy_to_double_conversion(){
    bp::converter::registry::push_back(& r_values_impl::Numpy2Double::convertible,
                                        & r_values_impl::Numpy2Double::construct,
                                        bp::type_id< double >());
}
void register_pysequence_to_indexvector_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2IndexArray::convertible,
                                        & r_values_impl::PySequence2IndexArray::construct,
                                        bp::type_id< GIMLI::IndexArray >());
}

void register_pysequence_to_ivector_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2IVector::convertible,
                                        & r_values_impl::PySequence2IVector::construct,
                                        bp::type_id< GIMLI::IVector >());
}

void register_pysequence_to_rvector_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2RVector::convertible,
                                        & r_values_impl::PySequence2RVector::construct,
                                        bp::type_id< GIMLI::Vector< double > >());
}
void register_pysequence_to_cvector_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2CVector::convertible,
                                        & r_values_impl::PySequence2CVector::construct,
                                        bp::type_id< GIMLI::Vector< GIMLI::Complex > >());
}

void register_pysequence_to_bvector_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2BVector::convertible,
                                        & r_values_impl::PySequence2BVector::construct,
                                        bp::type_id< GIMLI::Vector< bool > >());
}
void register_pysequence_to_StdVectorRVector3_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2StdVectorRVector3::convertible,
                                        & r_values_impl::PySequence2StdVectorRVector3::construct,
                                        bp::type_id< std::vector< GIMLI::Pos > >());
}
void register_pysequence_to_r3vector_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PySequence2R3Vector::convertible,
                                        & r_values_impl::PySequence2R3Vector::construct,
                                        bp::type_id< GIMLI::R3Vector >());
}
void register_pytuple_to_rvector3_conversion(){
    bp::converter::registry::push_back(& r_values_impl::PyTuple2RVector3::convertible,
                                        & r_values_impl::PyTuple2RVector3::construct,
                                        bp::type_id< GIMLI::Pos >());
}
void register_numpy_to_rmatrix_conversion(){
    bp::converter::registry::push_back(& r_values_impl::Numpy2RMatrix::convertible,
                                        & r_values_impl::Numpy2RMatrix::construct,
                                        bp::type_id< GIMLI::Matrix< double > >());
}
