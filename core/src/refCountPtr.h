/******************************************************************************
 *   Copyright (C) 2014-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_REFCOUNTPTR__H
#define _GIMLI_REFCOUNTPTR__H

#include "gimli.h"

namespace GIMLI{

// Base class for reference counted objects
// Scott Meyers: More Effective C++ Item 29 Source Code
class RCObject {
public:
    void addReference(){
        __M
        ++refCount_;
    }

    void removeReference(){
        __M
        if (--refCount_ == 0) delete this;
    }

    void markUnshareable(){
        shareable_ = false;
    }

    bool isShareable() const{
        return shareable_;
    }

    bool isShared() const {
        return refCount_ > 1;
    }

protected:
    RCObject()
        : refCount_(0), shareable_(true) {}

    RCObject(const RCObject& rhs)
        : refCount_(0), shareable_(true) {}

    RCObject& operator=(const RCObject& rhs){
        return *this;
    }

    virtual ~RCObject() = 0;

private:
    Index refCount_;
    bool shareable_;
};

 //*! Template class for smart pointers-to-T objects;
 /*! Template class for smart pointers-to-T objects;
  * must support the RCObject interface
  * */
template<class T > class RefCountIPtr {
public:

    RefCountIPtr(T * realPtr = 0)
    : counter_(new CountHolder){
        counter_->pointee_ = realPtr;
        init_();
    }

    RefCountIPtr(const RefCountIPtr & rhs)
    : counter_(rhs.counter_){
        init_();
    }

    RefCountIPtr(RefCountIPtr & rhs)
    : counter_(rhs.counter_){
        THROW_TO_IMPL
    }

    ~RefCountIPtr(){
        counter_->removeReference();
    }

    RefCountIPtr & operator = (const RefCountIPtr & rhs){
        __M
        if (counter_ != rhs.counter_) {
            counter_->removeReference();
            counter_ = rhs.counter_;
            init_();
        }
        return *this;
    }

    RefCountIPtr & operator = (RefCountIPtr & rhs){
        __M
        if (counter_ != rhs.counter_) {
            counter_->removeReference();
            counter_ = rhs.counter_;
            init_();
        }
        THROW_TO_IMPL
        return *this;
    }

    T * operator->() {
        __M
        THROW_TO_IMPL
        return counter_->pointee_;
    }

    T * operator->() const{
        __M
        return counter_->pointee_;
    }

    T & operator*(){
        __M
        THROW_TO_IMPL
        return *counter_->pointee_;
    }

    T & operator*() const {
        __M
        return *counter_->pointee_;
    }

private:
    struct CountHolder : public RCObject {
        ~CountHolder() { delete pointee_; }
        T * pointee_;
    };

    void init_(){
        if (counter_->isShareable() == false) {
            T *oldValue = counter_->pointee_;
            counter_ = new CountHolder;
            counter_->pointee_ = oldValue ? new T(*oldValue) : 0;
        }
        counter_->addReference();
    }

    CountHolder *counter_;

};

} // namespace GIMLI

#endif //_GIMLI_REFCOUNTPTR__H
