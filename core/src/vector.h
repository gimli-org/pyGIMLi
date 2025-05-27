/******************************************************************************
 *   Copyright (C) 2007-2025 by the GIMLi development team                    *
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
//** Idea taken from
//
// @Article{Veldhuizen95b,
//   author =       "Todd L. Veldhuizen",
//   title =        "Expression templates",
//   journal =      "C++ Report",
//   volume =       "7",
//   number =       "5",
//   pages =        "26--31",
//   month =        jun,
//   year =         "1995",
//   note =         "Reprinted in C++ Gems, ed. Stanley Lippman"
// }

#pragma once

#define EXPRVEC_USE_TEMPORARY_EXPRESSION

//#define EXPRVEC_USE_BOOST_THREAD // case 1
//#define EXPRVEC_USE_STD_ALGORITHM // case 2
#define EXPRVEC_USE_INDIRECTION // case 3 #current default

#include "gimli.h"
#include "expressions.h"
#include "utils.h"

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <initializer_list>
#include <numeric>
#include <cmath>
#include <cstring>
#include <fstream>
#include <cerrno>
#include <iterator>
#include <memory>

#ifdef USE_BOOST_BIND
    // deprecated
    #include <boost/bind.hpp>
#else
    #include <functional>
#endif

//#define __DDS(...) __MS(__VA_ARGS__)
#define __DDS(...)

namespace GIMLI{

template < class ValueType, class A > class __VectorExpr;

// forwarding inline
IndexArray find(const BVector & v);

DLLEXPORT IndexArray range(Index start, Index stop, Index step=1);
DLLEXPORT IndexArray range(Index stop);


#ifndef PYGIMLI_CAST
inline void Dump(const void * mem, unsigned int n) {
    const char * p = reinterpret_cast< const char *>(mem);
    for (unsigned int i = 0; i < n; i++) {
        std::cout << std::hex << int(p[i]) << " ";
    }
    std::cout << std::endl;
}
#endif

template < class ValueType > class DLLEXPORT VectorIterator {
public:
    typedef ValueType value_type;

    #ifdef _MSC_VER
    // basic iterator traits needed by msvc
    typedef Index difference_type;
    typedef value_type & pointer;
    typedef value_type & reference;
    typedef std::random_access_iterator_tag iterator_category;
    bool operator < (const VectorIterator< ValueType > & a) const { return val_ < a.val_; }
    #endif

    VectorIterator()
        : val_(0), maxSize_(0), end_(0), data_(0), _hasBorrowedData(false){
    }

    VectorIterator(ValueType * v, Index size, const Vector < ValueType > * owner)
        : val_(v), maxSize_(size), end_(v + size),
        data_(owner->data_), _hasBorrowedData(owner->_hasBorrowedData){
    }

    VectorIterator(const VectorIterator < ValueType > & iter)
        : val_(iter.val_), maxSize_(iter.maxSize_), end_(iter.val_ + iter.maxSize_), data_(iter.data_),
        _hasBorrowedData(iter._hasBorrowedData){
    }

    VectorIterator < ValueType > & operator = (const VectorIterator < ValueType > & iter){
        if (this != & iter){
//             __M
            val_ = iter.val_;
            maxSize_ = iter.maxSize_ ;
            end_ = iter.end_;
            data_ = iter.data_;
            _hasBorrowedData = iter._hasBorrowedData;
        }
        return *this;
    }

    inline const ValueType & operator * () const { return * val_; }

    inline ValueType & operator * () { return * val_; }

    inline const ValueType & operator [] (const Index i) const { //__MS(this) __MS(val_) __MS(*val_) Dump(val_, sizeof(ValueType));
        return val_[i]; }

    inline ValueType & operator [] (const Index i) {
        return val_[i]; }

    inline VectorIterator< ValueType > & operator ++ () { // prefix ++A
        ++val_; return *this;
    }

    inline VectorIterator< ValueType > & operator -- () { // prefix ++A
        --val_; return *this;
    }

    inline VectorIterator< ValueType > operator ++ (int) { // postfix A++
         VectorIterator< ValueType > old(*this); // make a copy for result
//          TO_IMPL
         ++(*this);            // Now use the prefix version to do the work
         return old;        // return the copy (the old) value.
    }

    inline VectorIterator< ValueType > operator -- (int) { // postfix A++
         VectorIterator< ValueType > old(*this); // make a copy for result
//          TO_IMPL
         --(*this);                // Now use the prefix version to do the work
         return old;           // return the copy (the old) value.
    }

    inline bool operator == (const VectorIterator< ValueType > & a) const { return val_ == a.val_; }
    inline bool operator != (const VectorIterator< ValueType > & a) const { return val_ != a.val_; }

    inline Index size() const { return maxSize_; }

    inline ValueType * ptr() const { return val_; }
    inline ValueType * ptr() { return val_; }

    inline bool hasMore() const { return val_ != end_; }

    /*!Return the current and advance the iterator to the next.*/
    inline ValueType nextVal(){
        // __MS(this) __MS(val_) __MS(*val_)
        // __MS(this->hasMore())
        return *val_++; }

    inline ValueType nextForPy(){// __MS(val_) __MS(*val_)
        if(!hasMore()){
#if PYGIMLI
#include "boost/python.hpp"
            boost::python::objects::stop_iteration_error();
#endif            //will not come here
        }
        return nextVal();
    }

    ValueType * val_;
    Index maxSize_;
    ValueType * end_;
    ValueType * data_; // begin()
    bool _hasBorrowedData;

};


template < class VEC1, class VEC2 > void _assertDifferentBorrowedView(const VEC1 & a,
                                                         const VEC2 & b){
//    if (a._hasBorrowedData && a.data_ == b.data_) {
//        log(Error, "vector a and b have same borrowed data and are the same single rowview of RDensematrix. This is probably not what you want." );
//    }
}

//! One dimensional array aka Vector of limited size.
/*!
One dimensional array aka Vector of limited size. Size limit depends on platform (32bit system maxsize = 2^32, 64bit system, maxsize=2^64)
*/
template< class ValueType > class DLLEXPORT Vector {
public:
    typedef ValueType ValType;
//    typedef VectorIterator< ValueType const > const_iterator;
    typedef VectorIterator< ValueType > iterator;

    /*!
     * Construct one-dimensional array of size n. The vector is cleaned (filled with zero)
     */
// #ifndef PYGIMLI_CAST
// this constructor is dangerous for IndexArray in pygimli ..
// there is an autocast from int -> IndexArray(int)
    Vector()
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(0)", this);
    // explicit Vector(Index n = 0) : data_(NULL), begin_(NULL), end_(NULL) {
        resize(0);
        clean();
    }
    Vector(Index n)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(int)", this);
    // explicit Vector(Index n = 0) : data_(NULL), begin_(NULL), end_(NULL) {
        resize(n);
        clean();
    }
// #endif
    /*!
     * Construct one-dimensional array of size n, and fill it with val
     */
    Vector(Index n, const ValueType & val)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(int, float)", this);
        resize(n);
        fill(val);
    }
    Vector(std::initializer_list< ValueType > l):
        _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec([])", this);
        resize(l.size());
        std::copy(l.begin(), l.end(), data_);
    }
    Vector(ValueType * buf, Index n):
        _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(*, n)", this);
        reserve(n);
        size_ = n;
        std::memcpy(data_, buf, sizeof(ValueType) * n);
    }

    /*!
     * Construct vector from file. Shortcut for Vector::load
     */
    Vector(const std::string & filename, IOFormat format=Ascii)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        this->load(filename, format);
    }
    /*!
     * Copy constructor. Create new vector as a deep copy of v.
     */
    Vector(const Vector< ValueType > & v)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(copy)", this, &v)
        // this->_hasBorrowedData = v._hasBorrowedData;
        // no shallow copy
        copy_(v);

        // if (v._hasBorrowedData){
            // // __MS(v.size(), v._borrowedDataOffset)
            // size_ = v.size();
            // this->_borrowedData = v._borrowedData;
            // this->_borrowedDataOffset = v._borrowedDataOffset;
            // this->data_ = &v._borrowedData[v._borrowedDataOffset];
        // } else {
        // }
        __DDS("vec copied:", this)
    }
    /*!
     * Copy constructor. Create new vector as a deep copy of the slice v[start, end)
     */
    Vector(const Vector< ValueType > & v, Index start, Index end)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(copy, start, end)", this);
        resize(end - start);
        std::copy(&v[start], &v[end], data_);
    }
    /*!
     * Copy constructor. Create new vector from expression
     */
    template < class A > Vector(const __VectorExpr< ValueType, A > & v)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(expr)", this);
        resize(v.size());
        assign_(v);
    }

#ifndef PYGIMLI_CAST
    /*!
     * Copy constructor. Create new vector as a deep copy of std::vector(Valuetype)
     */
    Vector(const std::vector< ValueType > & v)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("std::vec()", this);
        resize(v.size());
        for (Index i = 0; i < v.size(); i ++) data_[i] = v[i];
        //std::copy(&v[0], &v[v.size()], data_);
    }
    /*! Create Vector from borrowed data. */
    Vector(Index n, const std::shared_ptr< ValueType [] > & d, Index offset=0)
        : _hasBorrowedData(true), size_(n), data_(0), capacity_(0){
        __DDS("vec(boorowed)", this, n, d)
        this->_borrowedData = d;

        if (d == nullptr){
            log(Error, "borrowed data is empty.");
        } else {
            this->setBorrowedDataOffset(offset);
        }
        __DDS("vec(borrowed) created:", this, *this)
    }

    template < class ValueType2 > Vector(const Vector< ValueType2 > & v)
        : _hasBorrowedData(false), size_(0), data_(0), capacity_(0){
        __DDS("vec(vec2)", this)
        resize(v.size());
        for (Index i = 0; i < v.size(); i ++) data_[i] = ValueType(v[i]);
        //std::copy(&v[0], &v[v.size()], data_);
    }
#endif

    /*! Default destructor. */
    ~Vector() {
        __DDS("~vec", this);
        free_();
    }

    /*! Assignment operator. Creates a new vector as copy of v */
    Vector< ValueType > & operator = (const Vector< ValueType > & v) {
        __DDS("vec=", this);
        if (this != &v) {
            copy_(v);
        }
        return *this;
    }

#if USE_EIGEN3
    Vector< ValueType > & operator = (const Eigen::VectorXd & v);
    Vector< ValueType > & operator += (const Eigen::VectorXd & v);
    Vector< ValueType > & setVal(const Eigen::VectorXd & v, const IVector & ids);
    Vector< ValueType > & addVal(const Eigen::VectorXd & v, const IVector & ids);

#endif

    /*! Assignment operator. Creates a new vector as from expression. */
    template < class A > Vector< ValueType > & operator = (const __VectorExpr< ValueType, A > & v) {
        __DDS("vec = expr", this);
        assign_(v);
        return *this;
    }

    /*! Assignment operator. Fill the existing vector with val. Shortcut for fill. */
    Vector< ValueType > & operator = (const ValueType & val) {
        fill(val);
        return *this;
    }

    /*! Return a deep copy. For numpy compatibility.*/
    inline Vector < ValueType > copy() const {
        return *this;
    }

    inline const ValueType & operator[](const Index i) const {
        // ASSERT_THIS_SIZE(i)  // will not work for std::copy, std::sort etc.
        return data_[i];
    }

    inline ValueType & operator[](const Index i) {
        // ASSERT_THIS_SIZE(i)  // will not work for std::copy, std::sort etc.
        return data_[i];
    }

    inline const Vector < ValueType > operator[](const IndexArray & i) const { return this->get_(i); }

    inline Vector < ValueType > operator[](const IndexArray & i) {
        return this->get_(i); }

    inline Vector < ValueType > operator[](const BVector & b) {
        return this->get_(b); }

#ifndef PYGIMLI_CAST
    inline const ValueType & operator () (const Index i) const {
        // ASSERT_THIS_SIZE(i)
        return data_[i]; }

    inline ValueType & operator () (const Index i) {
        // ASSERT_THIS_SIZE(i)
        return data_[i]; }
#endif
     /*!
      * Return a new vector that match the slice [start, end).
      *  end == -1 or larger size() sets end = size.
      * Throws exception on violating boundaries.
      */
    inline Vector < ValueType > operator () (Index start, SIndex end) const {
        return getVal(start, end);
    }

    inline Vector < ValueType > operator () (const std::pair< Index, SIndex > & pair) const {
        return getVal(pair.first, pair.second);
    }

//     inline Vector < ValueType > operator () (const std::vector < int > & idx) const {
//         return get_(idx);
//     }
    /*!
     * Return a new vector that based on indices's.
     * Throws exception if indices's are out of bound
     */
    inline Vector < ValueType > operator () (const IndexArray & iArray) const {
        return get_(iArray);
    }
    /*!
     * Return a new vector that based on indices's.
     * Throws exception if indices's are out of bound.
     */
    inline Vector < ValueType > operator () (const SIndexArray & siArray) const {
        return get_(siArray);
    }
    inline Vector < ValueType > operator () (const IVector & iVec) const {
        return get_(iVec);
    }
    template < class IndexContainer > Vector < ValueType > get_(const IndexContainer & idx) const {
//         __MS(&idx)
//         __MS(idx.size())
//         //__MS(idx)
        Vector < ValueType > v(idx.size());
        Index id;
        for (Index i = 0; i < idx.size(); i ++){
           id = idx[i];
           if (id >= 0 && id < size_){
                v[i] = data_[(Index)id];
           } else {
                throwLengthError(WHERE_AM_I + " idx out of range " +
                                     str(id) + " [" + str(0) + " " + str(size_) + ")");
           }
        }
        return v;
    }

    /*! */
    Vector < ValueType > get_(const BVector & bv) const {
        return this->get_(GIMLI::find(bv));
    }
    Vector < ValueType > getVUI_(const IndexArray & iA) const {
        return this->get_(iA);
    }
    Vector < ValueType > getVSI_(const IVector & iA) const {
        return this->get_(iA);
    }

#ifndef PYGIMLI_CAST
    /*!
        Implicit converter for Vector< T > = Vector< ValueType >
    */
//     template < class T > operator Vector< T >() const {
//         //COUTMARKER
//         Vector< T > f(this->size());
//         for (Index i = 0; i < this->size(); i ++){ f[i] = T(data_[i]); }
//         return f;
//     }

    /*!
     *  Implicit converter for std::vector < T > = Vector< ValueType >
     */
    template < class T > operator std::vector< T >() {
        std::vector< T > f(this->size());
        for (Index i = 0; i < this->size(); i ++){ f[i] = T(data_[i]); }
        return f;
    }
#endif

    /*! Set the value to the whole array. Same as fill(val) */
    inline Vector< ValueType > & setVal(const ValueType & val) {
        this->fill(val);
        return *this;
    }

    /*! Set the val where bv is true.
     * Throws out of length exception if sizes dismatch.
     * Same as setVal(val, find(bv)) but faster. */
    inline Vector< ValueType > & setVal(const ValueType & val,
                                        const BVector & bv) {
        ASSERT_EQUAL(this->size(), bv.size())
        for (Index i = 0; i < bv.size(); i ++ ) if (bv[i]) data_[i] = val;
        return *this;
    }

    /*! Set the value val at index i.
     * Throws out of range exception if index is not in [0, size). */
    inline Vector< ValueType > & setVal(const ValueType & val,
                                        Index i) {
        ASSERT_RANGE(i, 0, this->size())
        data_[i] = val;
        return *this;
    }

    /*! Return valid slice pair [start, end) after size correction and range assertion.
    End will be set to max size for end > size or to size + end for end < 0.
    Start larger end throws range exception. */
    std::pair < Index, Index > slice(Index start, SIndex end) const {
        Index e = (Index) end;
        if (end < 0){
            e = max(start, size_ + end);
        }
        return slice(start, e);
    }
    /*! Return valid slice pair [start, end) after size correction and range assertion.
    End will be set to max size for end > size.
    Start larger end throws range exception. */
    std::pair < Index, Index > slice(Index start, Index end) const {
        Index e = min(end, size_);
        if (start != e) {
            ASSERT_RANGE(start, 0, e)
        }
        return std::pair < Index, Index >(start, e);
    }

    /*! Set values at slice. See \ref slice. */
    inline Vector< ValueType > & setVal(const ValueType & val,
                                        Index start, SIndex end) {
        std::pair < Index, Index > se(slice(start, end));
        std::fill(data_ + se.first, data_ + se.second, val);
        return *this;
    }

    /*!Set value val from pair.start to pair.end*/
    inline Vector< ValueType > & setVal(const ValueType & val,
                                        const std::pair< Index, SIndex > & pair) {
        return setVal(val, pair.first, pair.second);
    }

    /*! Set multiple values. Throws out of range exception if index check fails. */
    inline Vector< ValueType > & setVal(const ValueType & val,
                                        const IndexArray & ids) {
        for (Index i = 0; i < ids.size(); i ++) {
            setVal(val, ids[i]);
        }
        return *this;
    }

    /*! Set multiple values from vals at index position iArray.
     * Throws out of range exception if index check fails. */
    inline Vector< ValueType > & setVal(const Vector < ValueType > & vals,
                                        const IndexArray & ids) {
        ASSERT_EQUAL(vals.size(), ids.size())
        for (Index i = 0; i < ids.size(); i ++){
//            data_[iArray[i]] = vals[i];
           setVal(vals[i], ids[i]);
        }
        return *this;
    }

    /*! Insert vals from start index. Resize if necessary.*/
    inline Vector< ValueType > & setVal(const Vector < ValueType > & vals,
                                        Index start) {
        Index newS = start + vals.size();
        if (this->size() < newS) this->resize(newS);
        this->setVal(vals, start, newS);

        return *this;
    }

    /*! Set values at slice. See \ref slice. */
    inline Vector< ValueType > & setVal(const Vector < ValueType > & vals,
                                        Index start, SIndex end) {

        std::pair < Index, Index > se(slice(start, end));

        if (vals.size() < (se.second - se.first)){
            throwLengthError( WHERE_AM_I + " vals.size() < (end-start) " +
                              str(vals.size()) + " " + str(se.first) + " " + str(se.second)) ;
        }

        if (this->size() == vals.size()){
            std::copy(&vals[se.first], &vals[se.second], &data_[se.first]);
        } else {
            std::copy(&vals[0], &vals[se.second - se.first], &data_[se.first]);
        }
        return *this;
    }

    /*! Set all vals from pair.start to pair.end */
    inline Vector< ValueType > & setVal(const Vector < ValueType > & vals,
                                        const std::pair< Index, SIndex > & pair) {
        return setVal(vals, pair.first, pair.second);
    }

    inline Vector< ValueType > & push_back(const ValueType & v) {
        __DDS("vec.push_back", this, v);
        resize(size_ + 1);
        return setVal(v, size_ -1);
    }

    /*! Add values at slice. See \ref slice. */
    inline Vector< ValueType > & addVal(const Vector < ValueType > & vals,
                                        Index start, Index end) {
        std::pair < Index, Index > se(slice(start, end));

        if (vals.size() < se.second - se.first){
            throwLengthError(WHERE_AM_I + " vals.size() < (end-start) " +
                                str(vals.size()) + " " + str(start) + " " + str(end)) ;
        }

        if (this->size() == vals.size()){
            for (Index i = se.first; i < se.second; i ++) data_[i] += vals[i];
        } else {
            for (Index i = se.first; i < se.second; i ++) data_[i] += vals[i - se.first];
        }

        return *this;
    }
    /*! Subtract values at slice. See \ref slice. */
    inline Vector< ValueType > & subVal(const Vector < ValueType > & vals,
                                        Index start, Index end) {
        std::pair < Index, Index > se(slice(start, end));

        if (vals.size() < end - start){
            throwLengthError(WHERE_AM_I + " vals.size() < (end-start) " +
                                str(vals.size()) + " " + str(start) + " " + str(end)) ;
        }

        if (this->size() == vals.size()){
            for (Index i = start; i < end; i ++) data_[i] -= vals[i];
        } else {
            for (Index i = start; i < end; i ++) data_[i] -= vals[i - start];
        }

        return *this;
    }

    inline Vector< ValueType > & addVal(const Vector < ValueType > & vals,
                                    const std::pair< Index, SIndex > & pair) {
        return addVal(vals, pair.first, pair.second);
    }
    inline Vector< ValueType > & subVal(const Vector < ValueType > & vals,
                                    const std::pair< Index, SIndex > & pair) {
        return subVal(vals, pair.first, pair.second);
    }

    /*! Add values from vals at IndexArray idx.
     * Throws length exception if sizes of vals and idx mismatch. */
    inline Vector< ValueType > & addVal(const Vector < ValueType > & vals,
                                        const IndexArray & idx) {
        ASSERT_EQUAL(idx.size(), vals.size())
        //ASSERT_THIS_SIZE(max(idx))
        for (Index i = 0; i < idx.size(); i ++) data_[idx[i]] += vals[i];
        return *this;
    }
    /*! This[ids] += vals[:] * alpha[ids] */
    inline Vector< ValueType > & addVal(const Vector < ValueType > & vals,
                                        const IndexArray & ids,
                                        const Vector < ValueType > & alpha) {
        // assert(max(ids) < vals.size() && vals.size() == alpha.size() // expensive?
        ASSERT_EQUAL(this->size(), alpha.size())

        Index id;
        for (Index i = 0; i < ids.size(); i ++) {
            id = ids[i];
            data_[id] += vals[i] * alpha[id];
        }
        return *this;
    }
    /*! sub values from vals at IndexArray idx.
     * Throws length exception if sizes of vals and idx mismatch. */
    inline Vector< ValueType > & subVal(const Vector < ValueType > & vals,
                                        const IndexArray & idx) {
        ASSERT_EQUAL(idx.size(), vals.size())
        //ASSERT_THIS_SIZE(max(idx))
        for (Index i = 0; i < idx.size(); i ++) data_[idx[i]] -= vals[i];
        return *this;
    }

    /*! Add val to index idx.
     */
    inline Vector< ValueType > & addVal(const ValueType & val, Index i) {
        ASSERT_RANGE(i, 0, this->size())
        data_[i] += val;
        return *this;
    }
    /*! Subtract val to index idx.
     */
    inline Vector< ValueType > & subVal(const ValueType & val, Index i) {
        ASSERT_RANGE(i, 0, this->size())
        data_[i] -= val;
        return *this;
    }

    /*! Add Values from an ElementMatrix.*/
    void add(const ElementMatrix < double > & A, bool neg=false);

    /*! Add Values from an ElementMatrix. Optional scale with constant scalar. */
    void add(const ElementMatrix < double > & A, const double & scale, bool neg=false);

    /*! Add Values from an ElementMatrix. Optional scale constant Pos. */
    void add(const ElementMatrix < double > & A, const Pos & scale,
             bool neg=false);

    /*! Add Values from an ElementMatrix. Optional scale constant RMatrix. */
    void add(const ElementMatrix < double > & A, const RSmallMatrix & scale,
             bool neg=false);

    /*! DEPRECATED Bad design (Per node values need to be interpolated to quadrature points first.)*/
    void add(const ElementMatrix < double > & A,
             const Vector < double > & scale);

    /*! Get value for index i.
     * Throws out of range exception if index check fails. */
    inline const ValueType & getVal(Index i) const {
        ASSERT_RANGE(i, 0, this->size())
        return data_[i];
    }

    /*! Get ith. value from the back or return 0.0*/
    inline ValueType last(Index i=0) const {
        if (this->size() <= i) return 0.0;
        return data_[this->size()-(i+1)];
    }

    Vector < ValueType > getVal(Index start, SIndex end) const {
        std::pair < Index, Index > se(slice(start, end));

        Vector < ValueType > v(se.second - se.first);

        if (se.first == se.second) return v;
        std::copy(&data_[se.first], &data_[se.second], &v[0]);

        return v;
    }

    Vector < ValueType > getVal(const std::pair< Index, SIndex > & pair) const {
        return getVal(pair.first, pair.second);
    }

#ifdef PYGIMLI
//    needed for: /usr/include/boost/python/def_visitor.hpp
    bool operator < (const Vector< ValueType > & v) const { return false; }
#else
    BVector operator < (const Vector< ValueType > & v) const {
        ASSERT_EQUAL(this->size(), v.size())
        _assertDifferentBorrowedView(*this, v);

        BVector ret(this->size(), 0);

        std::less<ValueType> l;
        for (Index i = 0; i < v.size(); i ++) {
            ret[i] = l(data_[i], v[i]);
        }
        return ret;
    }
#endif

#define DEFINE_COMPARE_OPERATOR_VEC__(OP, FUNCT) \
    BVector operator OP (const Vector< ValueType > & v) const { \
        _assertDifferentBorrowedView(*this, v); \
        ASSERT_EQUAL(this->size(), v.size()) \
        BVector ret(this->size(), 0); \
        FUNCT<ValueType> f; \
        for (Index i = 0; i < v.size(); i ++){ret[i] = f(data_[i], v[i]);} \
        return ret; \
    } \

// DEFINE_COMPARE_OPERATOR_VEC__(<=, std::less_equal)
// DEFINE_COMPARE_OPERATOR_VEC__(>=, std::greater_equal)
// DEFINE_COMPARE_OPERATOR_VEC__(>, std::greater)

#undef DEFINE_COMPARE_OPERATOR_VEC__

#define DEFINE_COMPARE_OPERATOR__(OP, FUNCT) \
    inline BVector operator OP (const int & v) const { \
        BVector ret(this->size(), false); \
        FUNCT<ValueType> f; \
        for (Index i = 0; i < this->size(); i ++){ ret[i] = f(data_[i], ValueType(v)); } \
        return ret;\
    } \
    inline BVector operator OP (const uint & v) const { \
        BVector ret(this->size(), false); \
        FUNCT<ValueType> f; \
        for (Index i = 0; i < this->size(); i ++){ ret[i] = f(data_[i], ValueType(v)); } \
        return ret;\
    } \
    inline BVector operator OP (const ValueType & v) const { \
        BVector ret(this->size(), false); \
        FUNCT<ValueType> f; \
        for (Index i = 0; i < this->size(); i ++){ ret[i] = f(data_[i], v); } \
        return ret;\
    } \

// DEFINE_COMPARE_OPERATOR__(<, std::less)
// DEFINE_COMPARE_OPERATOR__(<=, std::less_equal)
// DEFINE_COMPARE_OPERATOR__(>=, std::greater_equal)
// DEFINE_COMPARE_OPERATOR__(==, std::equal_to)
// DEFINE_COMPARE_OPERATOR__(!=, std::not_equal_to)
// DEFINE_COMPARE_OPERATOR__(>, std::greater)

#undef DEFINE_COMPARE_OPERATOR__

//     inline BVector operator > (const ValueType & v) const {
//         BVector ret(this->size(), 0);
//         std::transform(data_, data_ + size_, &ret[0], boost::bind(isGreater< ValueType >, _1, v));
//         return ret;
//     }

#define DEFINE_UNARY_MOD_OPERATOR__(OP, FUNCT) \
  inline Vector< ValueType > & operator OP##= (const Vector < ValueType > & v) { \
        _assertDifferentBorrowedView(*this, v); \
        if (v.size() == 1) return *this OP##= v[0];\
        ASSERT_SIZES_GREATER_EQUAL((*this), v) \
        std::transform(data_, data_ + v.size(), &v[0], data_, FUNCT()); return *this; } \
  inline Vector< ValueType > & operator OP##= (const ValueType & val) { \
        for (Index i = 0; i < size_; i ++){data_[i] OP##= val;}return *this; } \

DEFINE_UNARY_MOD_OPERATOR__(+, PLUS)
DEFINE_UNARY_MOD_OPERATOR__(-, MINUS)
DEFINE_UNARY_MOD_OPERATOR__(/, DIVID)
DEFINE_UNARY_MOD_OPERATOR__(*, MULT)

#undef DEFINE_UNARY_MOD_OPERATOR__

    /*! Negation operator thats return a copy of this with negative values. */
    inline Vector < ValueType > operator - () const { return *this * -1.0; }

    /*! Return dot(this, b[off::]) */
    ValueType mult(const Vector< ValueType > & b, Index off=0) const {
        ASSERT_SIZE_GREATER_EQUAL(b, this->size() + off)
        ValueType s = ValueType(0.0);
        for (Index i = 0; i < this->size(); i ++) s+= data_[i] * b[i+off];
        return s;
    }


    /*! Resize if n differs size() and fill new with val. Old data are preserved. */
    void resize(Index n, ValueType fill){
        if (n != size_){
            if (this->_hasBorrowedData){
                throwError("Vector with borrowed data can't be resized.");
            }
            reserve(n);
            for (Index i = size_; i < n; i ++) data_[i]=fill;
            size_ = n;
        }
    }
    // default args lead win64 pygimli segfault .. WTF??
    void resize(Index n){
        resize(n, 0);
    }

    /*! Reserve memory. Old data are preserved*/
    void reserve(Index n){

        if (this->_hasBorrowedData){
            throwError("Vector with borrowed data can't be resized.");
        }

        Index newCapacity = max(1, n);
        if (capacity_ != 0){
            int exp;
            frexp(n, &exp);
            newCapacity = pow(2, exp);
        }
//         __MS(n << " " << capacity_ << " " << newCapacity)

        if (newCapacity != capacity_) {

            if (this->_hasBorrowedData){
                log(Error, "Cannot resize Vector with borrowed data");
                return;
            }

            ValueType * buffer = new ValueType[newCapacity];

            std::memcpy(buffer, data_,
                        sizeof(ValueType) * min(capacity_, newCapacity));
            // std::destroy_at(data_);
            if (data_) delete [] data_;
            data_  = buffer;
            capacity_ = newCapacity;
            //std::copy(&tmp[0], &tmp[min(tmp.size(), n)], data_);
        }
     }

    /*! Fill the whole vector from the pointer of val.
     *  CAUTION!! There is no boundary check.
     * Val must be properly ([val, val+this->size()))assigned.  */
    template< class V > Vector< ValueType > & fill(V * val) {
        for (Index i = 0; i < size_; i ++) data_[i] = ValueType(val[i]);
        //std::copy(val, val + size_, data_);
        return *this;
    }

    /*! Fill the whole vector with val. */
    Vector< ValueType > & fill(const ValueType & val) {
        std::fill(data_, data_ + size_, val); return *this; }

    /*! Fill the whole vector with function expr(i) */
    template< class Ex > Vector< ValueType > & fill(Expr< Ex > expr){
        for (Index i = 0; i < size_; i ++){
            data_[i] = expr((ValueType)i);
        }
        return *this;
    }

    inline void assign(const Vector< ValueType > & v){
        _assertDifferentBorrowedView(*this, v);
        this->copy_(v);
    }
    template < class ExprOP > inline void assign(const ExprOP & v){
        if (v.size()) {
            resize(v.size());
            v.assign(*this);
        }
    }
    //     /*! Fill the whole vector with function expr(i) */
//     template< class V > void fill(const V & val){
//         for (Index i = 0; i < size_; i ++) data_[i] = ValueType(val);
//     }

    /*! Fill Vector with 0.0. Don't change size.*/
    void clean(){
        if (size_ > 0) std::memset(data_, '\0', sizeof(ValueType) * size_);
    }

    /*! Empty the vector. Frees memory and resize to 0.*/
    void clear(){ free_(); }

    /*! Round all values of this array to a given tolerance. */
    Vector< ValueType > & round(const ValueType & tolerance){
        THROW_TO_IMPL
        return *this;
    }

    //VectorIterator< ValueType > * end() { return end_; }

    inline bool empty() const { return size() == 0; }

    inline Index size() const { return size_; }

    inline Index capacity() const { return capacity_; }

    inline Index nThreads() const { return nThreads_; }

    inline Index singleCalcCount() const { return singleCalcCount_; }

    ValueType * data() { return data_; }

    /*! Serialize vector into a byte buffer.
    Same structure like save a vector.*/
    ByteBuffer serialize() const {
        return GIMLI::serialize(*this);
    }

    /*! Deserialize the vector from a byte buffer.
    Same structure like loading a vector.*/
    void deserialize(const ByteBuffer & m) {
        return GIMLI::deserialize(m, *this);
    }

    /*! Write the vector to a stream.
    * Version =  0: [size(8), size*(valueType)]\n
    * Version = -1: [version(8), typeinfo(8), size(8), size*(valueType)]
    */
    void writeToStream(std::ostream & out, int version=-1) const {
        if (_hasBorrowedData){
            log(Error, "Vector with borrowed data can't "
                        "yet be written to stream.");
        }
        if (version == -1){
            writeStream(out, (int64)-1); // version info
            writeStream(out, (char)typeid(ValueType).name()[0]);
        }

        writeStream(out, (int64)size_);
        writeStream(out, data_[0], size_);
    }

    void readFromStream(std::istream & in) {
        this->clear();

        int64 version; readStream(in, version, 1);
        int64 size;
        if (version < 0) {
            // has version info for pg>=2
            char type; readStream(in, type, 1);
            if (type != typeid(ValueType).name()[0]) {
                throwError(WHERE_AM_I + " Type mismatch in readFromStream. "
                    "Expected: " + std::string(typeid(ValueType).name()) +
                    " but got: " + type);
            }
            readStream(in, size, 1);
        } else {
            size = version;
        }

        this->resize(size);
        readStream<ValType>(in, data_[0], size);
    }

    /*! Save the object to file. Returns true on success and in case of trouble an exception is thrown.
     * The IOFormat flag will be overwritten if the filename have a proper file suffix. Ascii format is forced if \ref VECTORASCSUFFIX is given.
     * Binary format is forced if \ref VECTORBINSUFFIX is set. If no suffix is provided \ref VECTORASCSUFFIX or \ref VECTORBINSUFFIX will be append. \n\n
    Binary format is: \n
    Unsigned int64 [8Byte] - length of the Array [0.. length) \n
    ValueType [sizeof(ValueType)] - 0th value \n
    ... \n
    ValueType [sizeof(ValueType)] - length -1 value \n\n
    Ascii format is a simple list of the values \n\n
    \param filename string of the file name
    \param IOFormat enum, either Ascii for human readable format, or Binary for fast save and load.
    */
    std::string save(const std::string & filename, IOFormat format=Ascii) const {

        if (filename.rfind(VECTORASCSUFFIX) != std::string::npos) {
            format = Ascii;
        } else if (filename.rfind(VECTORBINSUFFIX) != std::string::npos) {
            format = Binary;
        }
        std::string fname(filename);

        if (format == Ascii){
            if (fname.rfind(".") == std::string::npos) {
                fname += VECTORASCSUFFIX;
            }

            std::ofstream file; file.open(fname.c_str());
            if (!file) {
                throwError(filename + ": " + strerror(errno));
            }

            file.setf(std::ios::scientific, std::ios::floatfield);
            file.precision(14);

            for (Index i = 0, imax = size_; i < imax; i ++) file << data_[i] << std::endl;
            file.close();
        } else {
            return GIMLI::saveBinary(filename, *this, VECTORBINSUFFIX);
        }
        return fname;
    }

    /*!
     * Load the object from file. Returns true on success and in case of trouble an exception is thrown.
     * The IOFormat flag will be overwritten if the filename have a proper file suffix.
     * Ascii format is forced if \ref VECTORASCSUFFIX is given.
     * Binary format is forced if \ref VECTORBINSUFFIX is set.
     * See Vector< ValueType >::save for file format.
     */
    bool load(const std::string & filename, IOFormat format=Ascii){

        if (filename.rfind(VECTORASCSUFFIX) != std::string::npos) format = Ascii;
        else if (filename.rfind(VECTORBINSUFFIX) != std::string::npos) format = Binary;

        if (!fileExist(filename)){
            if (fileExist(filename + VECTORBINSUFFIX))
                return this->load(filename + VECTORBINSUFFIX, Binary);
            if (fileExist(filename + VECTORASCSUFFIX))
                return this->load(filename + VECTORASCSUFFIX, Ascii);
        }

        if (format == Ascii){
            std::vector < ValueType > tmp;

            std::fstream file; openInFile(filename.c_str(), &file);
            ValueType val; while(file >> val) {
                tmp.push_back(val);
            }

            this->resize(tmp.size());
            std::copy(tmp.begin(), tmp.end(), &data_[0]);
            file.close();

        } else {
            GIMLI::loadBinary(filename, *this, VECTORBINSUFFIX);
        }
        return true;
    }

    VectorIterator< ValueType > beginPyIter() const {
        return VectorIterator< ValueType >(data_, size_, this);
    }
    VectorIterator< ValueType > begin() const {
        return VectorIterator< ValueType >(data_, size_, this);
    }
    VectorIterator< ValueType > end() const {
        return VectorIterator< ValueType >(data_ + size_, 0, this);
    }

    Index hash() const {
        Index seed = 0;
        for (Index i = 0; i < this->size_; ++i) {
            hashCombine(seed, this->data_[i]);
        }
        return seed;
    }

    // std::shared_ptr< ValueType [] > & borrowedData() const {
    //     return _borrowedData; }

    inline void setBorrowedDataOffset(Index offset){
        data_ = &_borrowedData.get()[offset];
    }
    bool _hasBorrowedData; // we need a flag if borrowedData are used, it cannot be read from the shared_ptr in case an empty pointer was set.
    std::shared_ptr< ValueType [] > _borrowedData;

    Index size_;
    ValueType * data_;

protected:

    void free_(){
        size_ = 0;
        capacity_ = 0;
        if (this->_hasBorrowedData){
            // delete [] this->_borrowedData;
        } else {
            if (data_)  delete [] data_;
        }
        data_  = NULL;
    }

    void copy_(const Vector< ValueType > & v){
        if (v.size()) {
            // __MS(*this)
            resize(v.size());
                //"" check speed for memcpy here
                //std::memcpy(data_, v.data_, sizeof(ValType)*v.size());
                // memcpy is approx 0.3% faster but copy is extensively testet
                // cleanest solution needs iterator rewriting:
                // std::copy(v.begin(), v.end(), this->begin());
            // __MS(v)

            std::copy(&v[0], &v[v.size()], data_); // only works without bound check in subscription operator
            //  __MS(data_)
            //  __MS(*this)
        }
    }

    template < class ExprOP > inline void assign_(const ExprOP & v) {
        if (v.size()) {
            //std::cout << "assign_" << v.size() << std::endl;
            resize(v.size());
            v.assign(*this);
        }
    }

    Index capacity_;

    static const Index minSizePerThread = 10000;
    static const int maxThreads = 8;
    int nThreads_;
    Index singleCalcCount_;
};

template <> DLLEXPORT void Vector< Pos >::clean();

// /*! Implement specialized type traits in vector.cpp */
template <> DLLEXPORT void
Vector<double>::add(const ElementMatrix < double >& A, bool neg);
template <> DLLEXPORT void
Vector<double>::add(const ElementMatrix < double >& A, const double & a, bool neg);
template <> DLLEXPORT void
Vector<double>::add(const ElementMatrix < double >& A, const Pos & a, bool neg);
template <> DLLEXPORT void
Vector<double>::add(const ElementMatrix < double >& A, const RSmallMatrix & a,
                    bool neg);

template< typename ValueType > void
Vector< ValueType >::add(const ElementMatrix < double >& A,  bool neg){ THROW_TO_IMPL}
template< typename ValueType > void
Vector< ValueType >::add(const ElementMatrix < double >& A, const double & a, bool neg){THROW_TO_IMPL}
template< typename ValueType > void
Vector< ValueType >::add(const ElementMatrix < double >& A, const Pos & a, bool neg){THROW_TO_IMPL}
template< typename ValueType > void
Vector< ValueType >::add(const ElementMatrix < double >& A,
                         const RSmallMatrix & a,
                         bool neg){ THROW_TO_IMPL }

// removeme in V1.2, 20200727
template <> DLLEXPORT void Vector<double>::add(
    const ElementMatrix < double > & A, const RVector & a);
// removeme in V1.2, 20200727
template< typename ValueType > void Vector< ValueType >::add(
    const ElementMatrix < double >& A, const Vector< double> & a){THROW_TO_IMPL}


template < > inline
Vector< double > & Vector< double >::round(const double & tolerance){
#ifdef USE_BOOST_BIND
    std::transform(data_, data_ + size_, data_, boost::bind(roundTo< double >, _1, tolerance));
#else
    for (Index i = 0; i < size_; i ++) data_[i] = roundTo(data_[i], tolerance);
#endif
    return *this;
}

template< class ValueType, class Iter > class AssignResult{
public:
    AssignResult(Vector< ValueType > & a, const Iter & result, Index start, Index end)
     : a_(&a), iter_(result), start_(start), end_(end){
    }
    void operator()() {
        ValueType * iter = a_->begin().ptr();
        //std::cout << start_ << " " << end_ << std::endl;
        for (Index i = start_; i < end_; i++) iter[i] = iter_[i];
    }

    Vector< ValueType > * a_;
    Iter iter_;
    Index start_;
    Index end_;
};

struct BINASSIGN { template < class T > inline T operator()(const T & a, const T & b) const { return b; } };

template< class ValueType, class Iter > void assignResult(Vector< ValueType > & v, const Iter & result) {
#ifdef EXPRVEC_USE_TEMPORARY_EXPRESSION
  // Make a temporary copy of the iterator.  This is faster on segmented
  // architectures, since all the iterators are in the same segment.
    Iter result2 = result;
#else
    // Otherwise, cast away const (eek!).  No harmful side effects.
    Iter& result2 = (Iter&)result;
#endif


#ifdef EXPRVEC_USE_BOOST_THREAD
//DEPRECATED
#error DONT  USE THIS
    if (v.nThreads() == 1) {
      AssignResult< ValueType, Iter >(v, result, 0, v.size())();
    } else {
      boost::thread_group threads;
      for (Index i = 0; i < v.nThreads(); i ++){
	Index start = v.singleCalcCount() * i;
	Index end   = v.singleCalcCount() * (i + 1);
	if (i == v.nThreads() -1) end = v.size();
	threads.create_thread(AssignResult< ValueType, Iter >(v, result, start, end));
      }
      threads.join_all();
    }
#else // no boost thread
    #ifdef EXPRVEC_USE_STD_ALGORITHM

    std::transform(v.begin().ptr(),
                   v.end().ptr(),
                   result,
                   v.begin().ptr(), BINASSIGN());

    #else  // no std algo
        #ifdef EXPRVEC_USE_INDIRECTION
            ValueType * iter = v.begin().ptr();

            // Inlined expression
            for (Index i = v.size(); i--;) iter[i] = result2[i];
        #else
            ValueType * iter = v.begin().ptr();
            ValueType * end  = v.end().ptr();

            do {
                *iter = *result2;       // Inlined expression
                ++result2;
            } while (++iter != end);
        #endif
    #endif
#endif
}

template< class ValueType, class A > class __VectorExpr {
public:
    __VectorExpr(const A & a) : iter_(a) {
        this->data_ = a.data_;
        this->_hasBorrowedData = a._hasBorrowedData;
        // __MS("Expr", this->data_, this->_hasBorrowedData)
    }

    inline ValueType operator [] (Index i) const { return iter_[i]; }

    inline ValueType operator * () const { return *iter_; }

    inline void operator ++ () { ++iter_; }

    void assign(Vector< ValueType > & x) const { assignResult(x, *this); }

    inline Index size() const { return iter_.size(); }

    A * begin() { return iter_.begin(); }
    A * end() { return iter_.end(); }

// private:
    A iter_;
    ValueType * data_;
    bool _hasBorrowedData;
};

template< class ValueType, class A, class Op > class __VectorUnaryExprOp {
public:
    __VectorUnaryExprOp(const A & a) : iter_(a) {
        this->data_ = a.data_;
        this->_hasBorrowedData = a._hasBorrowedData;
        // __MS("V", this->data_, this->_hasBorrowedData)
    }

    inline ValueType operator [] (Index i) const { return Op()(iter_[i]); }

    inline ValueType operator * () const { return Op()(*iter_); }

    inline void operator ++ () { ++iter_;  }

    inline Index size() const { return iter_.size(); }

// private:
    A iter_;
    ValueType * data_;
    bool _hasBorrowedData;
};

template< class ValueType, class A, class B, class Op > class __VectorBinaryExprOp {
public:
    __VectorBinaryExprOp(const A & a, const B & b) : iter1_(a), iter2_(b) {
        _assertDifferentBorrowedView(a, b);
        this->data_ = a.data_;
        this->_hasBorrowedData = a._hasBorrowedData;
        // __MS("V*W", this->data_, this->_hasBorrowedData)
    }

    inline ValueType operator [] (Index i) const { return Op()(iter1_[i], iter2_[i]); }

    inline ValueType operator * () const { return Op()(*iter1_, *iter2_); }

    inline void operator ++ () { ++iter1_; ++iter2_; }

    inline Index size() const { return iter2_.size(); }

// private:
    A iter1_;
    B iter2_;
    ValueType * data_;
    bool _hasBorrowedData;
};

template< class ValueType, class A, class Op > class __VectorValExprOp {
public:
    __VectorValExprOp(const A & a, const ValueType & val) : iter_(a), val_(val) {
        this->data_ = a.data_;
        this->_hasBorrowedData = a._hasBorrowedData;
        // __MS("V*a", this->data_, this->_hasBorrowedData)
    }

    inline ValueType operator [] (Index i) const { return Op()(iter_[i], val_); }

    inline ValueType operator * () const { return Op()(*iter_, val_); }

    inline void operator ++ () { ++iter_; }

    inline Index size() const { return iter_.size(); }

// private:
    A iter_;
    ValueType val_;
    ValueType * data_;
    bool _hasBorrowedData;
};

template< class ValueType, class A, class Op > class __ValVectorExprOp {
public:
    __ValVectorExprOp(const ValueType & val, const A & a) : iter_(a), val_(val)
    {
        this->data_ = a.data_;
        this->_hasBorrowedData = a._hasBorrowedData;
        // __MS("a*V", this->data_, this->_hasBorrowedData)
    }

    inline ValueType operator [] (Index i) const { return Op()(val_, iter_[i]); }

    inline ValueType operator * () const { return Op()(val_, *iter_); }

    inline void operator ++ () { ++iter_; }

    inline Index size() const { return iter_.size(); }

// private:
    A iter_;
    ValueType val_;
    ValueType * data_;
    bool _hasBorrowedData;
};

#define DEFINE_UNARY_EXPR_OPERATOR__(OP, FUNCT)\
\
template < class T > \
__VectorExpr< T, __VectorUnaryExprOp< T, VectorIterator< T >, FUNCT > > \
OP(const Vector< T > & a){ \
    typedef __VectorUnaryExprOp< T, VectorIterator< T >, FUNCT > ExprT; \
    return __VectorExpr< T, ExprT >(ExprT(a.begin())); } \
\
template < class T, class A > \
__VectorExpr< T, __VectorUnaryExprOp< T, __VectorExpr< T, A >, FUNCT > > \
OP(const __VectorExpr< T, A > & a){ \
    typedef __VectorUnaryExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT; \
    return __VectorExpr< T, ExprT >(ExprT(a)); \
} \

DEFINE_UNARY_EXPR_OPERATOR__(abs,   ABS_)
DEFINE_UNARY_EXPR_OPERATOR__(acot,  ACOT)
DEFINE_UNARY_EXPR_OPERATOR__(atan,  ATAN)
DEFINE_UNARY_EXPR_OPERATOR__(cos,   COS)
DEFINE_UNARY_EXPR_OPERATOR__(cot,   COT)
DEFINE_UNARY_EXPR_OPERATOR__(exp,   EXP)
DEFINE_UNARY_EXPR_OPERATOR__(exp10, EXP10)
DEFINE_UNARY_EXPR_OPERATOR__(fabs,  ABS_)
DEFINE_UNARY_EXPR_OPERATOR__(log,   LOG)
DEFINE_UNARY_EXPR_OPERATOR__(log10, LOG10)
DEFINE_UNARY_EXPR_OPERATOR__(sign, SIGN)
DEFINE_UNARY_EXPR_OPERATOR__(sin,   SIN)
DEFINE_UNARY_EXPR_OPERATOR__(sqrt,  SQRT)
DEFINE_UNARY_EXPR_OPERATOR__(square, SQR)
DEFINE_UNARY_EXPR_OPERATOR__(tan,   TAN)
DEFINE_UNARY_EXPR_OPERATOR__(tanh,  TANH)

#undef DEFINE_UNARY_EXPR_OPERATOR__

#define DEFINE_EXPR_OPERATOR__(OP, FUNCT)				\
template < class T >							\
__VectorExpr< T, __VectorBinaryExprOp< T, VectorIterator< T >, VectorIterator< T >, FUNCT > > \
operator OP (const Vector< T > & a, const Vector< T > & b){		\
    typedef __VectorBinaryExprOp< T, VectorIterator< T >, VectorIterator< T >, FUNCT > ExprT; \
    _assertDifferentBorrowedView(a, b); \
    return __VectorExpr< T, ExprT >(ExprT(a.begin(), b.begin()));	\
}									\
                                                                        \
template < class T >							\
__VectorExpr< T, __VectorValExprOp< T, VectorIterator< T >, FUNCT > >	\
operator OP (const Vector< T > & a, const T & val){			\
  typedef __VectorValExprOp< T, VectorIterator< T >, FUNCT > ExprT;     \
  return __VectorExpr< T, ExprT >(ExprT(a.begin(), val));		\
}                                                                       \
                                                                        \
template < class T >							\
__VectorExpr< T, __ValVectorExprOp< T, VectorIterator< T >, FUNCT > >	\
operator OP (const T & val, const Vector< T > & a){			\
  typedef __ValVectorExprOp< T, VectorIterator< T >, FUNCT > ExprT;     \
  return __VectorExpr< T, ExprT >(ExprT(val, a.begin()));		\
}									\
    									\
template< class T, class A >						\
__VectorExpr< T, __VectorBinaryExprOp< T, __VectorExpr< T, A >, VectorIterator< T >, FUNCT > > \
operator OP (const __VectorExpr< T, A > & a, const Vector< T > & b){	\
  typedef __VectorBinaryExprOp< T, __VectorExpr< T, A >, VectorIterator< T >, FUNCT > ExprT; \
  _assertDifferentBorrowedView(a, b); \
  return __VectorExpr< T, ExprT >(ExprT(a, b.begin()));		\
}									\
    									\
template< class T, class A >					\
__VectorExpr< T, __VectorBinaryExprOp< T, VectorIterator< T >, __VectorExpr< T, A >, FUNCT > > \
operator OP (const Vector< T > & a, const __VectorExpr< T, A > & b){	\
  typedef __VectorBinaryExprOp< T, VectorIterator< T >, __VectorExpr< T, A >, FUNCT > ExprT; \
  _assertDifferentBorrowedView(a, b); \
  return __VectorExpr< T, ExprT >(ExprT(a.begin(), b));		\
}									\
									\
template< class T, class A >					\
__VectorExpr< T, __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
        operator OP (const __VectorExpr< T, A > & a, const T & val){	\
            typedef __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT; \
            return __VectorExpr< T, ExprT >(ExprT(a, val));	\
}								\
        \
template< class T, class A >				\
__VectorExpr< T, __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
operator OP (const T & val, const __VectorExpr< T, A > & a){	\
  typedef __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT; \
  return __VectorExpr< T, ExprT >(ExprT(val, a));		\
}								\
                                                                \
template< class T, class A, class B >				\
__VectorExpr< T, __VectorBinaryExprOp< T, __VectorExpr< T, A >, __VectorExpr< T, B >, FUNCT > > \
operator OP (const __VectorExpr< T, A > & a, const __VectorExpr< T, B > & b){ \
  typedef __VectorBinaryExprOp< T, __VectorExpr< T, A >, __VectorExpr< T, B >, FUNCT > ExprT; \
  _assertDifferentBorrowedView(a, b); \
  return __VectorExpr< T, ExprT >(ExprT(a, b));			\
}									\
\
template< class T, class T2, class A >					\
        __VectorExpr< T, __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
        operator OP (const __VectorExpr< T, A > & a, const T2 & val){	\
        typedef __VectorValExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT; \
        return __VectorExpr< T, ExprT >(ExprT(a, (T)val));		\
}								\
        \
template< class T, class T2, class A >				\
        __VectorExpr< T, __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > >	\
        operator OP (const T2 & val, const __VectorExpr< T, A > & a){	\
        typedef __ValVectorExprOp< T, __VectorExpr< T, A >, FUNCT > ExprT;	\
        return __VectorExpr< T, ExprT >(ExprT((T)val, a));		\
}								\
        \
template < class T, class T2 >					        \
        __VectorExpr< T, __ValVectorExprOp< T, VectorIterator< T >, FUNCT > >	\
        operator OP (const T2 & val, const Vector< T > & a){			\
        typedef __ValVectorExprOp< T, VectorIterator< T >, FUNCT > ExprT;  \
        return __VectorExpr< T, ExprT >(ExprT((T)val, a.begin()));		\
}									\
        \
template < class T, class T2 >						\
        __VectorExpr< T, __VectorValExprOp< T, VectorIterator< T >, FUNCT > >	\
        operator OP (const Vector< T > & a, const T2 & val){			\
        typedef __VectorValExprOp< T, VectorIterator< T >, FUNCT > ExprT;  \
        return __VectorExpr< T, ExprT >(ExprT(a.begin(), (T)val));	\
}                                                                       \
        \

DEFINE_EXPR_OPERATOR__(+, PLUS)
DEFINE_EXPR_OPERATOR__(-, MINUS)
DEFINE_EXPR_OPERATOR__(*, MULT)
DEFINE_EXPR_OPERATOR__(/, DIVID)

#undef DEFINE_EXPR_OPERATOR__

//********************************************************************************
//** define some utility functions

// inline bool operator < (const GIMLI::Vector<double>&a, const GIMLI::Vector<double> &b) {
//     return false;
// }


template < class ValueType >
bool operator == (const Vector< ValueType > & v1, const Vector< ValueType > & v2){
    if (v1.size() != v2.size()) return false;
    for (Index i = 0; i < v1.size(); i ++){
        if (!isEqual(v1[i], v2[i])) {
            // for (Index j = 0; j < v1.size(); j ++){
            //     __MS(j<< " " << v1[j] << "  "  << v2[j] <<  " " << v1[j]-v2[j] )
            // }
            return false;
        }
    }
    return true;
}

template < class ValueType, class A >
bool operator == (const Vector< ValueType > & v1, const __VectorExpr< ValueType, A > & v2){
    return v1 == Vector< ValueType >(v2);
}

/*! Return true if all values lower than TOLERANCE. */
template < class ValueType >
bool zero(const Vector< ValueType > & v){
    if (v.size() > 0) return (max(abs(v)) < TOLERANCE);
    return false;
}

/*! Return true if at least one value is greater than TOLERANCE. */
template < class ValueType >
bool nonZero(const Vector< ValueType > & v){
    if (v.size() == 0) return false;
    return !zero(v);
}

template < class ValueType >
bool operator != (const Vector< ValueType > & v1, const Vector< ValueType > & v2){
    return !(v1==v2);
}

/*! Return scalar product <v1, v2>. Redirect from \ref dot*/
template < class ValueType >
ValueType mult(const Vector< ValueType > & v1, const Vector< ValueType > & v2){
    return sum(v1*v2);
}

/*! Return scalar product <v1, v2>.*/
template < class ValueType >
ValueType dot(const Vector< ValueType > & v1, const Vector< ValueType > & v2){
    return mult(v1,v2);
}

//template double dot(const RVector & v1, const RVector & v2);
//template double mult(const RVector & v1, const RVector & v2);

/*! Find function. Return index vector of true values */
inline IndexArray find(const BVector & v){
    IndexArray idx;
    idx.reserve(v.size());
    for (Index i = 0; i < v.size(); i ++){
        if (v[i]) idx.push_back(i);
    }
    return idx;
}

/*! Refactor with expression templates */
inline BVector operator ~ (const BVector & a){
    BVector ret(a.size());
    for (Index i = 0; i < ret.size(); i ++) ret[i] = !a[i];
    return ret;
}

/*! For use in pygimli*/
inline BVector inv(const BVector & a){
    return ~a;
}

/*! Refactor with expression templates */
inline BVector operator & (const BVector & a, const BVector & b){
    BVector ret(a.size());
    for (Index i = 0; i < ret.size(); i ++) ret[i] = a[i] && b[i];
    return ret;
}

/*! Refactor with expression templates */
inline BVector operator | (const BVector & a, const BVector & b){
    BVector ret(a.size());
    for (Index i = 0; i < ret.size(); i ++) ret[i] = a[i] || b[i];
    return ret;
}

inline RVector operator * (const BVector & a, const RVector & b){
    RVector ret(a.size());
    for (Index i = 0; i < ret.size(); i ++) ret[i] = a[i] * b[i];
    return ret;
}

inline RVector operator * (const RVector & a, const BVector & b){
    RVector ret(a.size());
    for (Index i = 0; i < ret.size(); i ++) ret[i] = a[i] * b[i];
    return ret;
}

// /*! Refactor for speed */
// template < class ValueType, class A > BVector
// operator == (const __VectorExpr< ValueType, A > & vec, const ValueType & v){
//     return Vector< ValueType >(vec) == v;
// }

#define DEFINE_COMPARE_OPERATOR__(OP) \
template < class ValueType, class A > BVector \
operator OP (const __VectorExpr< ValueType, A > & vec, const ValueType & v){ \
    BVector ret(vec.size(), false); \
    for (Index i = 0; i < ret.size(); i ++) ret[i] = vec[i] OP v; \
    return ret;\
} \
template < class ValueType > BVector \
operator OP (const std::vector < ValueType > & vec, const ValueType & v){ \
    BVector ret(vec.size(), false); \
    for (Index i = 0; i < ret.size(); i ++) ret[i] = vec[i] OP v; \
    return ret;\
} \
template < class ValueType > BVector \
operator OP (const Vector < ValueType > & vec, const ValueType & v){ \
    BVector ret(vec.size(), false); \
    for (Index i = 0; i < ret.size(); i ++) { ret[i] = vec[i] OP v; }\
    return ret;\
} \
template < class ValueType > BVector \
operator OP (const Vector < ValueType > & vec, int v){ \
    BVector ret(vec.size(), false); \
    for (Index i = 0; i < ret.size(); i ++) ret[i] = vec[i] OP v; \
    return ret;\
} \

DEFINE_COMPARE_OPERATOR__(<)
DEFINE_COMPARE_OPERATOR__(<=)
DEFINE_COMPARE_OPERATOR__(>=)
DEFINE_COMPARE_OPERATOR__(==)
DEFINE_COMPARE_OPERATOR__(!=)
DEFINE_COMPARE_OPERATOR__(>)

#undef DEFINE_COMPARE_OPERATOR__


#define DEFINE_COMPARE_OPERATOR__(OP) \
template < class ValueType > BVector \
operator OP (const Vector < ValueType > & a, const Vector < ValueType > & b){ \
    BVector ret(a.size(), false); \
    ASSERT_EQUAL_SIZE(a, b) \
    for (Index i = 0; i < ret.size(); i ++) ret[i] = a[i] OP b[i]; \
    return ret;\
} \

DEFINE_COMPARE_OPERATOR__(<)
DEFINE_COMPARE_OPERATOR__(>)
DEFINE_COMPARE_OPERATOR__(<=)
DEFINE_COMPARE_OPERATOR__(>=)

#undef DEFINE_COMPARE_OPERATOR__




#define DEFINE_UNARY_COMPARE_OPERATOR__(OP, FUNCT) \
template < class ValueType, class A > \
BVector OP(const __VectorExpr< ValueType, A > & vec){ \
    BVector ret(vec.size(), 0); \
    for (Index i = 0; i < ret.size(); i ++) ret[i] = FUNCT()(vec[i]); \
    return ret; \
}\
template < class ValueType > \
BVector OP (const Vector< ValueType > & vec){ \
    BVector ret(vec.size(), 0); \
    std::transform(vec.begin().ptr(), vec.end().ptr(), ret.begin().ptr(), FUNCT()); \
    return ret; \
} \

DEFINE_UNARY_COMPARE_OPERATOR__(isInf, ISINF)
DEFINE_UNARY_COMPARE_OPERATOR__(isNaN, ISNAN)
DEFINE_UNARY_COMPARE_OPERATOR__(isInfNaN, ISINFNAN)

#undef DEFINE_UNARY_COMPARE_OPERATOR__



template < class T > Vector < T > cat(const Vector< T > & a, const Vector< T > & b){
    Vector < T > c (a.size() + b.size());
    std::copy(&a[0], &a[a.size()], &c[0]);
    std::copy(&b[0], &b[b.size()], &c[a.size()]);
    return c;
}

template < class T, class A > T sum(const __VectorExpr< T, A > & a){
    //std::cout << "sum(vectorExpr)" << std::endl;
    T tmp(0.0);
    for (Index i = 0, imax = a.size(); i < imax; i++) tmp += a[i];
    return tmp;

//     T tmp(0.0);
//     __VectorExpr< T, A > al = a;
//     for (Index i = 0; i < a.size(); i++, ++al) {
//         tmp += *al;
//     }
//     return tmp;

//     T tmp(0.0);
//     __VectorExpr< T, A > al = a;
//     for (; al != a.end(); ++al){
//         tmp += *al;
//     }
//     return tmp;


//return std::accumulate(a[0], a[a.size()], T());
}

//** Templates argue with python bindings
inline Complex sum(const CVector & c){
    return std::accumulate(c.begin(), c.end(), Complex(0));
}
inline double sum(const RVector & r){
    return std::accumulate(r.begin(), r.end(), double(0));
}
inline SIndex sum(const IVector & i){
    return std::accumulate(i.begin(), i.end(), SIndex(0));
}


template < class T, class A > T min(const __VectorExpr< T, A > & a){ return min(Vector< T >(a)); }
template < class T, class A > T max(const __VectorExpr< T, A > & a){ return max(Vector< T >(a)); }

inline Complex max(const CVector & v){
    ASSERT_NON_EMPTY(v)
    Complex ret=v[0];
    for (Index i = 1; i < v.size(); i ++ ) if (v[i] > ret) ret = v[i];
    return ret;
}

inline Complex min(const CVector & v){
    ASSERT_NON_EMPTY(v)
    Complex ret=v[0];
    for (Index i = 1; i < v.size(); i ++ ) if (v[i] < ret) ret = v[i];
    return ret;
}

template < class T > T min(const Vector < T > & v){
    ASSERT_NON_EMPTY(v)
    return *std::min_element(&v[0], &v[0] + v.size());
}
template < class T > T max(const Vector < T > & v){
    ASSERT_NON_EMPTY(v)
    return *std::max_element(&v[0], &v[0] + v.size());
}

template < class T > void capMax(Vector < T > & v, T max){
    ASSERT_NON_EMPTY(v)
    for (Index i = 0; i < v.size(); i ++ ) v[i] = min(v[i], max);
}

template < class T > void capMin(Vector < T > & v, T min){
    ASSERT_NON_EMPTY(v)
    for (Index i = 0; i < v.size(); i ++ ) v[i] = max(v[i], min);
}

template < class ValueType >
    ValueType mean(const Vector < ValueType > & a){
        return sum(a) / ValueType(a.size());
}

template < class ValueType, class A>
    ValueType mean(const __VectorExpr< ValueType, A > & a){
        return sum(a) / a.size();
}

template < class ValueType > ValueType stdDev(const Vector < ValueType > & a){
    return std::sqrt(sum(square(a - mean(a))) / (double)(a.size() - 1));
}

template < class ValueType > bool haveInfNaN(const Vector < ValueType > & v){
    for (VectorIterator < ValueType > it = v.begin(); it != v.end(); ++it){
        if (isInfNaN(*it)) return true;
    }
    return false;
}

template < class ValueType > Vector < ValueType > fixZero(const Vector < ValueType > & v, const ValueType tol = TOLERANCE){
    Vector < ValueType > ret(v);
    for (VectorIterator < ValueType > it = ret.begin(); it != ret.end(); ++it){
        if (::fabs(*it) < TOLERANCE) *it = tol;
    }
    return ret;
}

template < class ValueType > Vector < ValueType > round(const Vector < ValueType > & v, ValueType tol){
    return Vector< ValueType >(v).round(tol);
}

template < class T >
Vector < T > fliplr(const Vector < T > & v){
    Vector < T > n(v.size());
    for (Index i = 0; i < v.size(); i ++) n[i] = v[v.size() - 1 - i];
    return n;
}

template < class T, class A, class T2 > Vector < T > pow(const __VectorExpr< T, A > & a, T2 power){
    return pow(Vector< T >(a), power);
}

template < class T > Vector < T > pow(const Vector < T > & v, const Vector < T > & npower){
    ASSERT_EQUAL(v.size(), npower.size())

    Vector < T > r(v.size());
    for (Index i = 0; i < v.size(); i ++) r[i] = std::pow(v[i], T(npower[i]));
    return r;
}

template < class T > Vector < T > pow(const Vector < T > & v, double npower){
    Vector < T > r(v.size());
    for (Index i = 0; i < v.size(); i ++){
        r[i] = std::pow(v[i], T(npower));
    }
    return r;
}

// no template < int|double > since castxml interprets it as pow(vec,vec(int))
template < class T > Vector < T > pow(const Vector < T > & v, int npower){
    return pow(v, (double)npower);
}

template < class T > Vector< T > sort(const Vector < T > & a){
 #ifndef PYGIMLI_CAST
    std::vector < T > tmp(a.size(), 0.0) ;
    for (Index i = 0; i < a.size(); i ++) tmp[i] = a[i];
    std::sort(tmp.begin(), tmp.end());

    Vector < T > ret(tmp);
    return ret;
#endif // fixme .. implement me without std::vector
//     Vector < T > t(a);
//     std::sort(t.begin(), t.end());
    return Vector < T > (0);
}

/*! Returning a copy of the vector and replacing all consecutive occurrences
    of a value by a single instance of that value. e.g. [0 1 1 2 1 1] -> [0 1 2 1].
    To remove all double values from the vector use an additionally sorting.
    e.g. unique(sort(v)) gets you [0 1 2]. */
template < class T > Vector< T > unique(const Vector < T > & a){
#ifndef PYGIMLI_CAST
    std::vector < T > tmp(a.size()), u;
    for (Index i = 0; i < a.size(); i ++) tmp[i] = a[i];
    std::unique_copy(tmp.begin(), tmp.end(), back_inserter(u));

    Vector < T > ret(u);
    return ret;
    #endif // fixme .. implement me without std::vector
    return Vector < T >(0);
}


#ifndef PYGIMLI_CAST
//** Beware! this is not thread safe.
template< class ValueType > struct indexCmp {
    indexCmp(const Vector < ValueType > & arr) : arr_(arr) {}
    bool operator()(const Index a, const Index b) const {
        return arr_[a] < arr_[b];
    }
    const Vector < ValueType > & arr_;
};
#endif

template < class ValueType >
void sort(const Vector < ValueType > & unsorted,
          Vector < ValueType > & sorted,
          IndexArray & indexMap){
#ifndef PYGIMLI_CAST
    indexMap.resize(unsorted.size());

    for (Index i=0; i < unsorted.size(); i++) indexMap[i] = i;

    // Sort the index map, using unsorted for comparison

    std::vector < Index > tmp(indexMap.size(), 0.0) ;
    for (Index i = 0; i < indexMap.size(); i ++) tmp[i] = indexMap[i];
    std::sort(tmp.begin(), tmp.end(), indexCmp< ValueType >(unsorted));
    for (Index i = 0; i < indexMap.size(); i ++) indexMap[i] = tmp[i];

    sorted = unsorted(indexMap);
#endif
}

template < class ValueType >
IndexArray sortIdx(const Vector < ValueType > & unsorted){
    IndexArray indexMap;
    Vector < ValueType > sorted;
    sort(unsorted, sorted, indexMap);
    return indexMap;
}


// template < template < class T > class Vec, class T >
// std::ostream & operator << (std::ostream & str, const Vec < T > & vec){
//     for (Index i = 0; i < vec.size(); i ++) str << vec[i] << " ";
//     return str;
// }

template < class T >
std::ostream & operator << (std::ostream & str, const std::vector < T > & vec){
    for (Index i = 0; i < vec.size(); i ++) str << vec[i] << " ";
    return str;
}

template < class T >
std::ostream & operator << (std::ostream & str, const Vector < T > & vec){
    for (Index i = 0; i < vec.size(); i ++) str << vec[i] << " ";
    return str;
}

template < class ValueType, class A >
std::ostream & operator << (std::ostream & str,
                            const __VectorExpr< ValueType, A > & vec){
    for (Index i = 0; i < vec.size(); i ++) str << vec[i] << " ";
    return str;
}

template < class ValueType >
std::istream & operator >> (std::istream & str,
                            Vector< ValueType > & vec){
    THROW_TO_IMPL
    return str;
}

/*!
Return a RVector with increasing values of size(n+1) filled with : [0, i*a + (i-1)*x, .. ,last]
*/
template < class ValueType >
Vector< ValueType > increasingRange2(const ValueType & a,
                                     const ValueType & last, Index n){
    if (abs(a) < 1e-12){
        throwError("Can't create increasing range for start value of: " +
        str(a) );
    }

    if (a == last){
        throwError("start and end need to be different: " +
        str(a) + " "  + str(last));
    }

    if (sign(a) != sign(last)){
        throwError("Can't create increasing range from [0 " + str(a) + " to " + str(last) + "]");
    }

    if (n < 3){
        throwError("need at least n > 2" + str(a) + " n(" + str(n) +") "+ str(last));
    }


    ValueType x = (last- (n-1) * a) / (n-2);

    Placeholder x__;
    RVector y(n); y.fill(x__ * a);
    for (Index i = 2; i < n; i ++ ) y[i] += x*(i-1);

    return y;
}

/*!
Return a RVector with increasing values of size(n+1) filled with : 0, first, ... ,last
*/
template < class ValueType >
Vector< ValueType > increasingRange(const ValueType & first,
                                    const ValueType & last, Index n){
    if (sign(first) != sign(last)){
        throwError("cant increase range from [0 " + str(first) + " to " + str(last) + "]");
    }

    Placeholder x__;
    RVector y(n + 1); y.fill(x__);

    ValueType dy = (last - first * n) / (sum(y) - ValueType(n));

    if (dy < 0.0){
//         std::cout << "decreasing number of layers: " << n << " " << dy << std::endl;
        return increasingRange(first, last, n-1);
    }

    ValueType yval = 0.0;
    for (Index i = 0; i < n; i ++){
        yval = yval + first + dy * i;
        y[i + 1] = yval;
    }
    //y = fliplr(y) * -1.0 //#+ g.max(y)
    return y;
    //return fliplr(y) * -1.0;
}

template < class ValueType >
Vector < std::complex < ValueType > > toComplex(const Vector < ValueType > & re,
                                                const Vector < ValueType > & im){
    Vector < std::complex < ValueType > > cv(re.size());
    for (Index i = 0; i < cv.size(); i ++)
        cv[i] = std::complex < ValueType >(re[i], im[i]);
    return cv;
}

// template < class ValueType >
// Vector < std::complex < ValueType > > toComplex(ValueType re, const Vector < ValueType > & im){
//     return toComplex(Vector < ValueType > (im.size(), re), im);
// }

// template < class ValueType >
// Vector < std::complex < ValueType > > toComplex(const Vector < ValueType > & re, ValueType im=0){
//     return toComplex(re, Vector < ValueType > (re.size(), im));
// }

inline CVector toComplex(const RVector & re, double im=0.){
     return toComplex(re, RVector(re.size(), im));
}
inline CVector toComplex(double re, const RVector & im){
     return toComplex(RVector(im.size(), re), im);
}

/*! Convert absolute and phase (default in mrad) values to complex values.
* To get the vice versa use abs(cvector) and phase(cvector). */
inline CVector polarToComplex(const RVector & mag, const RVector & phi,
                              bool mRad=false){
    log(Warning, "polarToComplex .. Do not use me" );
    ASSERT_EQUAL(mag.size(), phi.size())
    if (mRad){
        return polarToComplex(mag, phi / 1000.0, false);
    } else {
        return toComplex(RVector(mag * cos(phi)), RVector(-mag * sin(phi)));
    }
}

template < class ValueType >
Vector < std::complex < ValueType > > operator * (const Vector < std::complex< ValueType > > & cv,
                                                   const Vector < ValueType > & v){
    return cv * toComplex(v);
}

template < class ValueType >
Vector < std::complex < ValueType > > operator * (const Vector < ValueType > & v,
                                                   const Vector < std::complex< ValueType > > & cv){
    return cv * toComplex(v);
}

template < class ValueType >
Vector < std::complex < ValueType > > operator / (const std::complex< ValueType > & v,
                                                  const Vector < std::complex< ValueType > > & cv){
    Vector < std::complex< ValueType > > ret(cv.size());
    for (Index i = 0; i < ret.size(); i ++ ) {
        ret[i] = v / cv[i];
    }
    return ret;
}

template < class ValueType >
Vector < std::complex < ValueType > > operator / (const ValueType & v,
                                                  const Vector < std::complex< ValueType > > & cv){
    return std::complex< ValueType >(v) / cv;
}

template < class ValueType, class A >
Vector < ValueType > real(const __VectorExpr< std::complex< ValueType >, A > & a){
    return real(Vector < std::complex< ValueType > >(a));
}

template < class ValueType >
Vector < ValueType > real(const Vector < std::complex< ValueType > > & cv){
    Vector < ValueType > v(cv.size());
    for (Index i = 0; i < cv.size(); i ++) v[i] = cv[i].real();
    return v;
}

template < class ValueType, class A >
Vector < ValueType > imag(const __VectorExpr< std::complex< ValueType >, A > & a){
    return imag(Vector < std::complex< ValueType > >(a));
}

template < class ValueType >
Vector < ValueType > imag(const Vector < std::complex< ValueType > > & cv){
    Vector < ValueType > v(cv.size());
    for (Index i = 0; i < cv.size(); i ++) v[i] = cv[i].imag();
    return v;
}

template < class ValueType, class A >
Vector < ValueType > angle(const __VectorExpr< std::complex< ValueType >, A > & a){
    return angle(Vector < std::complex< ValueType > >(a));
}

template < class ValueType >
Vector < ValueType > angle(const Vector < std::complex< ValueType > > & z){
    Vector < ValueType > v(z.size());
    for (Index i = 0; i < z.size(); i ++) v[i] = std::atan2(imag(z[i]), real(z[i]));
    return v;
}

inline RVector angle(const RVector & b, const RVector & a){
    ASSERT_EQUAL_SIZE(b, a)
    RVector v(b.size());
    for (Index i = 0; i < b.size(); i ++) v[i] = std::atan2(b[i], a[i]);
    return v;
}

template < class ValueType >
Vector < ValueType > phase(const Vector < std::complex< ValueType > > & z){
    Vector < ValueType > v(z.size());
    for (Index i = 0; i < z.size(); i ++) v[i] = std::arg(z[i]);
    return v;
}


template < class ValueType, class A >
Vector < ValueType > abs(const __VectorExpr< std::complex< ValueType >, A > & a){
    return abs(Vector < std::complex< ValueType > >(a));
}

template < class ValueType >
Vector < ValueType > abs(const Vector < std::complex< ValueType > > & cv){
    return sqrt(real(cv * conj(cv)));
}

template < class ValueType, class A >
Vector < std::complex< ValueType > > conj(const __VectorExpr< std::complex< ValueType >, A > & a){
    return conj(Vector < std::complex< ValueType > >(a));
}

template < class ValueType >
Vector < std::complex< ValueType > > conj(const Vector < std::complex< ValueType > > & cv){
    Vector < std::complex< ValueType > > v(cv.size());
    for (Index i = 0; i < cv.size(); i ++) v[i] = conj(cv[i]);
    return v;
}

inline RVector TmpToRealHACK(const RVector & v){ return v; }
inline RVector TmpToRealHACK(const CVector & v){ __M return real(v); }

#define DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR(OP) \
template <class T, class U > \
inline std::complex< T > operator OP (const std::complex< T > & lhs, const U & rhs) { \
    std::complex< T > ret; \
    return ret OP##= rhs; \
} \

DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR(*)
DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR(/)
DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR(-)
DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR(+)

#undef DEFINE_SCALAR_COMPLEX_BINARY_OPERATOR

inline IVector toIVector(const RVector & v){
    IVector ret(v.size());
    for (Index i = 0; i < ret.size(); i ++) ret[i] = int(v[i]);
    return ret;
}


template < class ValueType >
std::string save(const Vector< ValueType > & a,
                 const std::string & filename,
                 IOFormat format=Ascii){
    return a.save(filename, format);
}

template < class ValueType >
bool load(Vector< ValueType > & a, const std::string & filename,
          IOFormat format = Ascii,
          bool verbose=true){
    return a.load(filename, format);
}

} // namespace GIMLI

#ifndef PYGIMLI_CAST
namespace std {
    template<> struct hash< std::complex < double > > {
        GIMLI::Index operator()(const std::complex < double > & p) const noexcept {
            return GIMLI::hash(p.real(), p.imag());
        }
    };
    template<> struct hash< GIMLI::RVector > {
        GIMLI::Index operator()(const GIMLI::RVector & p) const noexcept {
            return p.hash();
        }
    };
    template<> struct hash< GIMLI::IndexArray > {
        GIMLI::Index operator()(const GIMLI::IndexArray & p) const noexcept {
            return p.hash();
        }
    };
    template<> struct hash< GIMLI::IVector > {
        GIMLI::Index operator()(const GIMLI::IVector & p) const noexcept {
            return p.hash();
        }
    };
    template<> struct hash< std::map< std::string, GIMLI::RVector > > {
        GIMLI::Index operator()(const std::map< std::string, GIMLI::RVector > & p) const noexcept {
            GIMLI::Index seed = 0;
            for (auto & x: p) {
                hashCombine(seed, x.first, x.second);
            }
            return seed;
        }
    };
}
#endif // PYGIMLI_CAST

