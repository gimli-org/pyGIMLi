/******************************************************************************
 *   Copyright (C) 2005-2024 by the GIMLi development team                    *
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

#ifndef _GIMLI_CHOLMODWRAPPER__H
#define _GIMLI_CHOLMODWRAPPER__H

#include "gimli.h"
#include "solverWrapper.h"

namespace GIMLI{

class DLLEXPORT CHOLMODWrapper : public SolverWrapper {
public:
    CHOLMODWrapper(RSparseMatrix & S, bool verbose=false, int stype=-2, bool forceUmfpack=false);
    CHOLMODWrapper(CSparseMatrix & S, bool verbose=false, int stype=-2, bool forceUmfpack=false);

    virtual ~CHOLMODWrapper();

    static bool valid();

    virtual void setMatrix(RSparseMatrix & S);

    virtual void setMatrix(CSparseMatrix & S);

    virtual void solve(const RVector & rhs, RVector & solution);

    virtual void solve(const CVector & rhs, CVector & solution);

protected:
    void init();

    void free();

    int initializeMatrix_(RSparseMatrix & S);

    int initializeMatrix_(CSparseMatrix & S);

    template < class ValueType >
    void init_(SparseMatrix < ValueType > & S, int stype);

    template < class ValueType >
    int initMatrixChol_(SparseMatrix < ValueType > & S, int xType);

    template < class ValueType >
    void solveCHOL_(const Vector < ValueType > & rhs, Vector < ValueType > & solution);

    template < class ValueType >
    void solveUmf_(const Vector < ValueType > & rhs, Vector < ValueType > & solution);

    int factorise_();

    int stype_;

    void *c_;
    void *A_;
    void *L_;

    bool useUmfpack_;
    bool forceUmfpack_;
    void *Numeric_;
    void *NumericD_;
    int * Ap_;
    int * Ai_;
    int * ApR_;
    int * AiR_;

    RVector *AxV_;
    RVector *AzV_;
};

} //namespace GIMLI;

#endif // _GIMLI_CHOLMODWRAPPER__H
