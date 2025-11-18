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

#include "cholmodWrapper.h"
#include "vector.h"
#include "sparsematrix.h"

#if CHOLMOD_FOUND
    #define CHOLMOD_MAXMETHODS 9
    #define UF_long long

    #include <cholmod.h>
    #define USE_CHOLMOD 1

    #if defined(UMFPACK_FOUND)
        #if (UMFPACK_FOUND == TRUE)
            #include <umfpack.h>
            #define USE_UMFPACK 1
        #endif
    #endif
#endif

namespace GIMLI{

#if USE_CHOLMOD
    bool CHOLMODWrapper::valid() { return true; }
#else
    bool CHOLMODWrapper::valid() { return false; }
#endif

CHOLMODWrapper::CHOLMODWrapper(RSparseMatrix & S, bool verbose, int stype,
                               bool forceUmfpack)
    : SolverWrapper(verbose), stype_(stype), forceUmfpack_(forceUmfpack){
    
    Numeric_ = nullptr;
    NumericD_= nullptr;
    c_ = NULL;
    A_ = NULL;
    L_ = NULL;
    AxV_ = NULL;
    AzV_ = NULL;
    Ap_ = NULL;
    Ai_ = NULL;
    ApR_ = NULL;
    AiR_ = NULL;
    setMatrix(S);
}
CHOLMODWrapper::CHOLMODWrapper(CSparseMatrix & S, bool verbose, int stype,
                               bool forceUmfpack)
    : SolverWrapper(verbose), stype_(stype), forceUmfpack_(forceUmfpack){
    Numeric_ = nullptr;
    NumericD_= nullptr;
    c_ = NULL;
    A_ = NULL;
    L_ = NULL;
    AxV_ = NULL;
    AzV_ = NULL;
    Ap_ = NULL;
    Ai_ = NULL;
    ApR_ = NULL;
    AiR_ = NULL;
    setMatrix(S);
}

CHOLMODWrapper::~CHOLMODWrapper(){
    this->free();
}

void CHOLMODWrapper::free(){
#if USE_CHOLMOD
    if (L_) cholmod_free_factor((cholmod_factor**)(&L_), (cholmod_common*)c_);
//     ** We did not allocate the matrix so we dont need to free it
//      cholmod_free_sparse(&A_, c_);

    cholmod_finish((cholmod_common*)c_);

    if (A_) delete (cholmod_sparse*)A_;
    A_ = nullptr;
    if (c_) delete (cholmod_common*)c_;
    c_ = nullptr;

#if USE_UMFPACK
    if (Numeric_) umfpack_zi_free_numeric (&Numeric_);
    Numeric_ = nullptr;
    if (NumericD_) umfpack_di_free_numeric (&NumericD_);
    NumericD_ = nullptr;
#endif
    // try allways
    if (AxV_) delete AxV_;
    AxV_ = nullptr;
    if (AzV_) delete AzV_;
    AzV_ = nullptr;

    if (ApR_) delete [] ApR_;
    ApR_ = nullptr;

    if (AiR_) delete [] AiR_;
    AiR_ = nullptr;

#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

void CHOLMODWrapper::setMatrix(RSparseMatrix & S){
    this->free();
    init_(S, stype_);
}
void CHOLMODWrapper::setMatrix(CSparseMatrix & S){
    this->free();
    init_(S, stype_);
}

template < class ValueType >
void CHOLMODWrapper::init_(SparseMatrix < ValueType > & S, int stype){

    dim_ = S.size();
    nVals_ = S.nVals();

    useUmfpack_ = false;
    Numeric_ = 0;
    NumericD_ = 0;
    c_ = NULL;
    A_ = NULL;
    L_ = NULL;
    AxV_ = NULL;
    AzV_ = NULL;
    Ap_ = NULL;
    Ai_ = NULL;
    ApR_ = NULL;
    AiR_ = NULL;

    if (stype == -2){
        stype_ = S.stype();
    } else {
        stype_ = stype;
    }

#if USE_CHOLMOD
    int ret = 0;
    c_ = new cholmod_common;
    ret =  cholmod_start((cholmod_common*)c_);
    if (ret) dummy_ = false;
    initializeMatrix_(S);
#elif USE_UMFPACK
    initializeMatrix_(S);
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
}

int CHOLMODWrapper::initializeMatrix_(CSparseMatrix & S){

    if (!dummy_){
        // check for non-hermetian
        if (S.stype() == 0){ //  matrix is full
            for (Index i = 0; i < S.size(); i++){
                for (int j = S.vecColPtr()[i]; j < S.vecColPtr()[i + 1]; j ++){
                    if (std::imag(S.vecVals()[j]) != 0.0){
                        if (S.vecVals()[j] == S.getVal(S.vecRowIdx()[j], i)){
                            // non-hermetian symmetric
                            // if (verbose_) std::cout << "non-hermetian symmetric "
                            // " matrix found .. switching to umfpack." << std::endl;
                            useUmfpack_ = true;
                            i = S.size();
                            break;
                        }
                    }
                }
            }
        }

        if (forceUmfpack_) useUmfpack_ = true;

        if (useUmfpack_){
#if USE_UMFPACK
// __MS("new C")
            Ap_ = (int*)S.colPtr();
            Ai_ = (int*)S.rowIdx();
            AxV_ = new RVector(real(S.vecVals()));
            AzV_ = new RVector(imag(S.vecVals()));
            double * Ax_ = &((*AxV_)[0]);
            double * Az_ = &((*AzV_)[0]);

            double *null = (double *) NULL;
            void *Symbolic;

            // if (verbose_) std::cout << "Using umfpack .. " << std::endl;

            umfpack_zi_symbolic (S.nRows(), S.nRows(), Ap_, Ai_, Ax_, Az_, &Symbolic, null, null) ;
            umfpack_zi_numeric (Ap_, Ai_, Ax_, Az_, Symbolic, &Numeric_, null, null) ;
            umfpack_zi_free_symbolic (&Symbolic);
            name_ = "Umfpack";
            return 1;
#else
            std::cerr << WHERE_AM_I << " umfpack not installed" << std::endl;
#endif
        } else {
#if USE_CHOLMOD
            name_ = "Cholmod";
            return initMatrixChol_(S, CHOLMOD_COMPLEX);
#endif
			return 0;
        }
    } // if ! dummy
    return 0;
}

int CHOLMODWrapper::initializeMatrix_(RSparseMatrix & S){
    if (!dummy_){
        // check symmetry
        if (S.stype() == 0){ //  matrix is full
            for (Index i = 0; i < S.size(); i++){
                for (int j = S.vecColPtr()[i]; j < S.vecColPtr()[i + 1]; j ++){
//                     __MS(i << " " << j << " " << S.vecColPtr()[i] << " " << S.vecColPtr()[i + 1] << " " << S.vecRowIdx()[j])
                    if (S.vecVals()[j] != 0.0){
//                 __MS(S.vecVals()[j] << " "  << S.getVal(S.vecRowIdx()[j], i))
                        if (::fabs(S.vecVals()[j] - S.getVal(S.vecRowIdx()[j], i, false)) > 1e-12){
                            // non-symmetric
                            if (verbose_) {
                                log(Info,
                                    "non-symmetric matrix found (", i,
                                    S.vecRowIdx()[j],
                                    S.vecVals()[j] , " ",
                                    S.getVal(S.vecRowIdx()[j], i), "tol:",
                                    ::fabs(S.vecVals()[j] - S.getVal(S.vecRowIdx()[j], i, false)),
                                    ").. switching to umfpack." );
                            }
                            useUmfpack_ = true;
                            i = S.size();
                            break;
                        }
                        //i = S.size();
                        //break;
                    }
                }
            }
        }
        
        if (forceUmfpack_ || getEnvironment("BERTUSEUMFPACK", 0) == 1) useUmfpack_ = true;

        if (useUmfpack_){
#if USE_UMFPACK

            log(Info, "Using umfpack.");
            int * ApT = (int*)S.colPtr();
            int * AiT = (int*)S.rowIdx();
            double * AxT = &S.vecVals()[0];

            AxV_ = new RVector(S.vecVals());

            double * null = (double *) NULL;

            ApR_ = new int[S.vecColPtr().size()];
            AiR_ = new int[S.vecRowIdx().size()];
            double *Ax_ =  new double[S.vecVals().size()]; // delete me!!
            //our crs format need to be transposed first

            int *P=0, *Q=0;
            (void) umfpack_di_transpose(S.nRows(), S.nCols(), ApT, AiT, AxT, P, Q, ApR_, AiR_, Ax_);

            for (uint i = 0; i < S.vecVals().size(); i++) (*AxV_)[i] = Ax_[i];

            void *Symbolic;

            // if (verbose_) std::cout << "Using umfpack .. " << std::endl;
            // beware transposed matrix here
            Stopwatch sw;
            // __MS(sw.duration(true))
            (void) umfpack_di_symbolic(S.nCols(), S.nRows(), ApR_, AiR_, Ax_, &Symbolic, null, null) ;
            // __MS(sw.duration(true)) // test maybe just 1% ??
            (void) umfpack_di_numeric(ApR_, AiR_, Ax_, Symbolic, &NumericD_, null, null) ;
            // __MS(sw.duration(true))
            umfpack_di_free_symbolic(&Symbolic);
            // __MS(sw.duration(true))
            name_ = "Umfpack";

            return 1;
#else
            std::cerr << WHERE_AM_I << " umfpack not installed" << std::endl;
#endif
        } else {
#if USE_CHOLMOD
            name_ = "Cholmod";
            return initMatrixChol_(S, CHOLMOD_REAL);
#endif
        }
        return 0;
    } // if ! dummy
    return 0;
}

template < class ValueType >
int CHOLMODWrapper::initMatrixChol_(SparseMatrix < ValueType > & S, int xType){
    if (!dummy_){
#if USE_CHOLMOD

        A_ = new cholmod_sparse;
        ((cholmod_sparse*)A_)->nrow  = S.nRows();           /* number of rows */
        ((cholmod_sparse*)A_)->ncol  = S.nCols();           /* number of columns */
        ((cholmod_sparse*)A_)->nzmax = S.nVals();           /* maximum number of entries */
        ((cholmod_sparse*)A_)->p     = (void*)S.colPtr();   /* column pointers (size n+1) or col indices (size nzmax) */
        ((cholmod_sparse*)A_)->i     = (void*)S.rowIdx();   /* row indices, size nzmax */

        ((cholmod_sparse*)A_)->x     = S.vals();     /* numerical values, size nzmax */
        ((cholmod_sparse*)A_)->stype = stype_;

        ((cholmod_sparse*)A_)->itype = CHOLMOD_INT;
        ((cholmod_sparse*)A_)->xtype = xType;  // data type for the pattern (Real, complex, zcomplex)
        ((cholmod_sparse*)A_)->dtype = CHOLMOD_DOUBLE; // data type for complex or real (float/double)
        ((cholmod_sparse*)A_)->packed = true;
        ((cholmod_sparse*)A_)->sorted = true; // testen, scheint schneller, aber hab ich das immer?

        factorise_();
#else
        std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
    } // if ! dummy
    return 0;
}

int CHOLMODWrapper::factorise_(){
    if (!dummy_){
#if USE_CHOLMOD
         if (useUmfpack_){
             // allready factorized in matrix init;
             //std::cerr << WHERE_AM_I << " factorize for umfpack called" << std::endl;
         } else {
            if (verbose_) cholmod_print_sparse((cholmod_sparse *)A_, "A", (cholmod_common*)c_);
            // cholmod_print_common("common parameter", (cholmod_common*)c_);


            // CHOLMOD_SIMPLICIA == 0
            // CHOLMOD_AUTO == 1
            // CHOLMOD_SUPERNODAL == 2
            // __MS(((cholmod_common *)c_)->supernodal)
            // ((cholmod_common *)c_)->supernodal=1;

            // __MS(((cholmod_common *)c_)->nmethods)
            // ((cholmod_common *)c_)->nmethods=0;
            // ((cholmod_common *)c_)->current=3;


            L_ = cholmod_analyze((cholmod_sparse*)A_,
                                 (cholmod_common*)c_);		    /* analyze */
            // __MS(((cholmod_factor *)L_)->is_super)
            cholmod_factorize((cholmod_sparse*)A_,
                              (cholmod_factor*)L_,
                              (cholmod_common*)c_);		    /* factorize */

        if (verbose_) std::cout << "CHOLMOD analyzed preordering: "
                                << ((cholmod_factor *)(L_))->ordering << std::endl;

        if (verbose_) cholmod_print_factor((cholmod_factor *)L_, "L", (cholmod_common*)c_);
    }
    return 1;
#else
        std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
    }
    return 0;
}

template < class ValueType >
    void CHOLMODWrapper::solveCHOL_(const Vector < ValueType > & rhs,
                                   Vector < ValueType > & solution){
    ASSERT_VEC_SIZE(rhs, this->dim_)
    ASSERT_VEC_SIZE(solution, this->dim_)
    if (!dummy_){
#if USE_CHOLMOD
        cholmod_dense * b = cholmod_zeros(((cholmod_sparse*)A_)->nrow,
                                         1,
                                         ((cholmod_sparse*)A_)->xtype,
                                         (cholmod_common*)c_);
        cholmod_dense * r = cholmod_zeros(((cholmod_sparse*)A_)->nrow,
                                          1,
                                          ((cholmod_sparse*)A_)->xtype,
                                          (cholmod_common*)c_);

        ValueType * bx = (ValueType*)b->x;
        for (uint i = 0; i < dim_; i++) bx[i] = rhs[i];

        cholmod_dense * x = cholmod_solve(CHOLMOD_A,
                                          (cholmod_factor *)L_,
                                          b,
                                          (cholmod_common *)c_);       /* solve Ax=b */

        if (((cholmod_sparse*)A_)->stype == 0){
            double al[2] = {0,0}, be[2] = {1,0};       /* basic scalars */
            cholmod_sdmult((cholmod_sparse*)A_, 0, be, al, x, r, (cholmod_common*)c_);
            bx = (ValueType *)r->x; /* ret = Ax */
            for (uint i = 0; i < dim_; i++) solution[i] = conj(bx[i]);

            //conj here .. check crs->ccs format or transpose before use

        } else {
            bx = (ValueType *)x->x; /* ret = x */
            for (uint i = 0; i < dim_; i++) solution[i] = bx[i];
        }
        cholmod_free_dense(&x, (cholmod_common*)c_);
        cholmod_free_dense(&r, (cholmod_common*)c_);
        cholmod_free_dense(&b, (cholmod_common*)c_);
#else
    std::cerr << WHERE_AM_I << " cholmod not installed" << std::endl;
#endif
    }
}

void CHOLMODWrapper::solve(const RVector & rhs, RVector & solution){
    ASSERT_VEC_SIZE(rhs, this->dim_)
    ASSERT_VEC_SIZE(solution, this->dim_)
    if (!dummy_){

        if (useUmfpack_){
#if USE_UMFPACK
            double * Ax_ = &((*AxV_)[0]);

            double *null = (double *) NULL ;

            (void) umfpack_di_solve(UMFPACK_A, ApR_, AiR_, Ax_,
                                    &solution[0], &rhs[0],
                                    NumericD_, null, null) ;


// int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
// int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
// double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
// // double b [ ] = {8., 45., -3., 3., 19.} ;
//
// void *Symbolic, *Numeric ;
// (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
// (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
// umfpack_di_free_symbolic (&Symbolic) ;
//
// (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, &solution[0], &rhs[0], Numeric, null, null) ;
// umfpack_di_free_numeric (&Numeric) ;
//

#endif
        } else {
            return solveCHOL_(rhs, solution);
        }
    }
}

void CHOLMODWrapper::solve(const CVector & rhs, CVector & solution){
    ASSERT_VEC_SIZE(rhs, this->dim_)
    ASSERT_VEC_SIZE(solution, this->dim_)
    if (!dummy_){

        if (useUmfpack_){
#if USE_UMFPACK
            double * Ax_ = &((*AxV_)[0]);
            double * Az_ = &((*AzV_)[0]);

            double *null = (double *) NULL ;

            RVector xx(rhs.size());
            RVector xz(rhs.size());
            RVector bx(real(rhs));
            RVector bz(imag(rhs));

            umfpack_zi_solve(UMFPACK_A, Ap_, Ai_, Ax_, Az_,
                             &xx[0], &xz[0], &bx[0], &bz[0],
                             Numeric_, null, null) ;

            // fix crs->ccs sparseformat and cleanup/unify umfpackwrapper for this

            // set booth to NULL or they will be wrongly deleted ..

            solution = toComplex(xx, xz);
#endif
        } else {
            return solveCHOL_(rhs, solution);
        }
    }
}


} //namespace GIMLI
