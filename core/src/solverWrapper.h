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

#ifndef _GIMLI_SOLVERWRAPPER__H
#define _GIMLI_SOLVERWRAPPER__H

#include "gimli.h"

namespace GIMLI{

class DLLEXPORT SolverWrapper{
public:
    
    SolverWrapper(bool verbose=false);

    virtual ~SolverWrapper();

    virtual void setMatrix(const RSparseMatrix & S){THROW_TO_IMPL;};

    virtual void setMatrix(const CSparseMatrix & S){THROW_TO_IMPL;};

    virtual void solve(const RVector & rhs, RVector & solution) = 0;

    virtual void solve(const CVector & rhs, CVector & solution){ THROW_TO_IMPL;}

    std::string name() const { return name_; }
    
protected:

    std::string name_;
    bool dummy_;
    bool verbose_;
    bool isComplex_;
    uint dim_;
    long nVals_;
    double dropTol_;
    double tolerance_;
    double maxiter_;

};

} //namespace GIMLI;

#endif // _GIMLI_SOLVERWRAPPER__H
