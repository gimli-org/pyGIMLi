/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *   Thomas Günther thomas@resistivity.net                                    *
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

#ifndef _GIMLI_INVERSIONBASE__H
#define _GIMLI_INVERSIONBASE__H

#include "gimli.h"

namespace GIMLI{

//! Inversion base template
/*! Inversion base template, provide inversion interface. */
template < class ModelValType > class InversionBase {
public:

    typedef Vector < ModelValType > ModelVector;

    /*! Default constructor*/
    InversionBase(){
    }

    /*! Default destructor */
    virtual ~InversionBase(){
    }

    /*! Abstract method for running the whole inversion. */
    virtual const ModelVector & run() = 0;

    virtual void setModel(const ModelVector & model) = 0;

    virtual void setReferenceModel(const ModelVector & model) = 0;

    /*! Abstract method for returning a const reference to the model vector. */
    virtual const ModelVector & model() const = 0;

    virtual void setData(const ModelVector & data)   = 0;

    virtual void setError(const ModelVector & err, bool isRelative=true)   = 0;

    virtual void setTransData(Trans< ModelVector > & t)   = 0;

    virtual void setTransModel(Trans< ModelVector > & t)   = 0;

    virtual void setLambda(double l) = 0;

    virtual void setMaxIter(int maxiter) = 0;

    virtual ModellingBase * forwardOperator() = 0;

    virtual void setForwardOperator(ModellingBase & fop) = 0;

    virtual const ModelVector & cWeight() const = 0;

    virtual void abort() = 0;

    virtual uint iter() const = 0;

    virtual double chi2() const = 0;

    virtual bool isRunning() const = 0;


protected:

};

} // namespace GIMLI

#endif // _GIMLI_INVERSIONBASE__H
