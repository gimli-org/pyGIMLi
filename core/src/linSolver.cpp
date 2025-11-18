/******************************************************************************
 *   Copyright (C) 2007-2024 by the GIMLi development team                    *
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

#include "linSolver.h"
#include "sparsematrix.h"
#include "ldlWrapper.h"
#include "cholmodWrapper.h"

namespace GIMLI{

RVector LinSolver::operator()(const RVector & rhs) {
    return this->solve(rhs);
}

CVector LinSolver::operator()(const CVector & rhs) {
    return this->solve(rhs);
}

LinSolver::LinSolver(bool verbose)
    : SolverWrapper(verbose){
     
    init_();
    setSolverType(AUTOMATIC);
}

LinSolver::LinSolver(RSparseMatrix & S, bool verbose)
    : SolverWrapper(verbose) {
    init_();
    setSolverType(AUTOMATIC);
    setMatrix(S);
}

LinSolver::LinSolver(RSparseMapMatrix & S, bool verbose)
    : SolverWrapper(verbose) {
    init_();
    setSolverType(AUTOMATIC);
    cacheMatrix_ = new RSparseMatrix(S);

    setMatrix(dynamic_cast < RSparseMatrix &>( *cacheMatrix_));
}

LinSolver::LinSolver(RSparseMatrix & S, SolverType solverType, bool verbose)
    : SolverWrapper(verbose) {
    init_();
    setSolverType(solverType);
    setMatrix(S);
}

LinSolver::LinSolver(CSparseMatrix & S, bool verbose)
    : SolverWrapper(verbose) {
    init_();
    setSolverType(AUTOMATIC);
    setMatrix(S);
}

LinSolver::LinSolver(CSparseMatrix & S, SolverType solverType, bool verbose)
    : SolverWrapper(verbose) {
    init_();
    setSolverType(solverType);
    setMatrix(S);
}

void LinSolver::init_(){
    rows_ = 0;
    cols_ = 0;
    solver_ = 0;
    cacheMatrix_ = 0;
}

LinSolver::~LinSolver(){
    if (cacheMatrix_) {
        delete cacheMatrix_;
        cacheMatrix_ = 0;
    }
    if (solver_) {
        delete solver_;
        solver_ = 0;
    }
}

void LinSolver::setSolverType(SolverType solverType){
   solverType_ = solverType;
   if (solverType_ == AUTOMATIC){
        solverType_ = UNKNOWN;

        if (LDLWrapper::valid()){
            solverType_ = LDL;
        }
        if (CHOLMODWrapper::valid()){
            solverType_ = CHOLMOD;
        }
    }
}

void LinSolver::setMatrix(RSparseMatrix & S, int stype){
    initialize_(S, stype);
}

void LinSolver::setMatrix(CSparseMatrix & S, int stype){
    initialize_(S, stype);
}

void LinSolver::solve(const RVector & rhs, RVector & solution){
    ASSERT_VEC_SIZE(rhs, cols_)
    solution.resize(rows_);
    if (rhs.size() != cols_){
        std::cerr << WHERE_AM_I << " rhs size mismatch: " << cols_ << "  " << rhs.size() << std::endl;
    }
    if (solver_) solver_->solve(rhs, solution);
}

void LinSolver::solve(const CVector & rhs, CVector & solution){
    ASSERT_VEC_SIZE(rhs, cols_)
    solution.resize(rows_);
    if (rhs.size() != cols_){
        std::cerr << WHERE_AM_I << " rhs size mismatch: " << cols_ << "  " << rhs.size() << std::endl;
    }
    if (solver_) solver_->solve(rhs, solution);
}

RVector LinSolver::solve(const RVector & rhs){
    ASSERT_VEC_SIZE(rhs, cols_)
    RVector solution(rhs.size());
    if (solver_) solver_->solve(rhs, solution);
    return solution;
}

CVector LinSolver::solve(const CVector & rhs){
    ASSERT_VEC_SIZE(rhs, cols_)
    CVector solution(rhs.size());
    if (solver_) solver_->solve(rhs, solution);
    return solution;
}

void LinSolver::initialize_(RSparseMatrix & S, int stype){
    rows_ = S.rows();
    cols_ = S.cols();
    setSolverType(solverType_);

    switch(solverType_){
        case LDL:     solver_ = new LDLWrapper(S, verbose_); break;
        case CHOLMOD: solver_ = new CHOLMODWrapper(S, verbose_, stype); break;
        case UMFPACK: solver_ = new CHOLMODWrapper(S, verbose_, stype, true); break;
        case UNKNOWN:
    default:
        std::cerr << WHERE_AM_I << " no valid solver found"  << std::endl;
    }
}

void LinSolver::initialize_(CSparseMatrix & S, int stype){
    rows_ = S.rows();
    cols_ = S.cols();
    setSolverType(solverType_);

    switch(solverType_){
        case LDL:     solver_ = new LDLWrapper(S, verbose_); break;
        case CHOLMOD: solver_ = new CHOLMODWrapper(S, verbose_, stype); break;
        case UMFPACK: solver_ = new CHOLMODWrapper(S, verbose_, stype, true); break;
        case UNKNOWN:
    default:
        std::cerr << WHERE_AM_I << " no valid solver found"  << std::endl;
    }
}

std::string LinSolver::solverName() const {
    if (solver_) return solver_->name();
    __MS("no solver initialized")

    switch(solverType_){
        case LDL:     return "LDL"; break;
        case CHOLMOD: return "CHOLMOD"; break;
        case UMFPACK: return "UMFPACK"; break;
        case UNKNOWN:
        default: return " no valid solver installed";
    }
}

} // namespace GIMLI

