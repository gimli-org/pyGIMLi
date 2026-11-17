#ifndef PYGIMLI__NUMPY__HPP
#define PYGIMLI__NUMPY__HPP

#include <numpy/arrayobject.h>

#include "gimli.h"


static int __numpy_initialized = 0;

// this fails for win since py39
// const static int __numpy_initialized = initNumpy();

inline int initNumpy(){
    // Needed or py* checks will segfault
    if (__numpy_initialized == 0){
        import_array2("Cannot import numpy.core.multiarray c-api for value converters.", 0);
        __numpy_initialized = 1;
    }
    return 0;
}

namespace GIMLI{


template < class ValueType > int npType(){
    __MS(sizeof(ValueType))
    __MS(typeid(ValueType).name())

    THROW_TO_IMPL
    return 0;
}

template < > inline int npType<bool>(){
    return NPY_BOOL;
}
template < > inline int npType<double>(){
    return NPY_DOUBLE;
}
template < > inline int npType<GIMLI::Complex>(){
    return NPY_COMPLEX128;
}
template < > inline int npType<GIMLI::Index>(){
    return NPY_UINT64;
}
template < > inline int npType<int>(){
    // int might by ambiguous, check here .. MAC vs. WIN cs. LIN
    return NPY_INT32;
}
template < > inline int npType<uint8>(){
    return NPY_UINT8;
}

template < class ValueType, template < class > class Vec >
PyObject * toNumpyArray(Vec< ValueType > & vec){
    initNumpy();

    npy_intp length = (ssize_t)vec.size();

    //** this fails with typename==int!
    // PyObject * ret = PyArray_SimpleNewFromData(1, &length,
    //                                            npType<ValueType>(),
    //                                            (void *)(&vec[0]));
    //** this fails with typename==int!

    PyObject * ret = PyArray_SimpleNew(1, &length, npType<ValueType>());
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(ValueType));
    return ret;
}

// special case for std::vector < ValueType, ....>
template < class ValueType, class Iter, template < class, class > class Vec >
PyObject * toNumpyArray(Vec< ValueType, Iter > & vec){
    // __MS("++++toNumpyArray")
    initNumpy();

    npy_intp length = (ssize_t)vec.size();

    PyObject * ret = PyArray_SimpleNew(1, &length, npType<ValueType>());
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject*>(ret)),
                (void *)(&vec[0]), length * sizeof(ValueType));
    // __MS("----toNumpyArray")
    return ret;
}

template < class BufferClass >
void fromNumpyArrayIntoBuffer(BufferClass & cls, PyObject * obj){
    // __MS("++++fromNumpyArrayIntoBuffer")
    cls.buf().clear();

    // __MS(obj, "\t is array:", obj->ob_type->tp_name)
    // __MS(obj, "\t is array", PyObject_TypeCheck(obj, &PyGenericArrType_Type))

    if (strcmp(obj->ob_type->tp_name, "numpy.ndarray") == 0){

        PyArrayObject *arr = (PyArrayObject *)obj;

        if (PyArray_ISONESEGMENT(arr) and PyArray_TYPE(arr) == NPY_UINT8){

            // __MS("\t", obj, "\t ndarray.length: ", PyArray_DIM(arr, 0))
            cls.buf().resize(PyArray_DIM(arr, 0));
            void * arrData = PyArray_DATA(arr);
            std::memcpy(&cls.buf()[0], arrData, cls.size()*sizeof(uint8));

        } else {

            __MS("\t", obj, "\t ndarray.ndim: ", PyArray_NDIM(arr))
            __MS("\t", obj, "\t ndarray.dtype: ", PyArray_TYPE(arr))
            __MS("\t", obj, "\t ndarray.length: ", PyArray_DIM(arr, 0))
            __MS("\t", obj, "\t NPY_UINT64: ", NPY_UINT64)
            __MS("\t", obj, "\t NPY_INT64: ", NPY_INT64)
            __MS("\t", obj, "\t NPY_UINT8: ", NPY_UINT8)
            __MS("\t", obj, "\t ndarray.onsegment: ", PyArray_ISONESEGMENT(arr))
            THROW_TO_IMPL
        }

    } else {
        log(Error, "BufferClass can only be filled with numpy.ndarray.");
    }
    // __MS("----fromNumpyArrayIntoBuffer")
}


} // namespace GIMLI
#endif // PYGIMLI__NUMPY__HPP^
