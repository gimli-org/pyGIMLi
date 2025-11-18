/******************************************************************************
 *   Copyright (C) 2007-2024 by the GIMLi development team                    *
 *   Carsten Ruecker carsten@resistivity.net                                  *
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

#ifndef _GIMLI_LINSOLVER__H
#define _GIMLI_LINSOLVER__H

#include "gimli.h"
#include "solverWrapper.h"
#include <cmath>
#include <map>

namespace GIMLI{

class SolverWrapper;

enum SolverType{AUTOMATIC,LDL,CHOLMOD,UMFPACK,UNKNOWN};


class DLLEXPORT LinSolver : public SolverWrapper{
public:
    LinSolver(bool verbose=false);
    
    LinSolver(RSparseMatrix & S, bool verbose=false);

    LinSolver(RSparseMapMatrix & S, bool verbose=false);

    LinSolver(CSparseMatrix & S, bool verbose=false);

    LinSolver(RSparseMatrix & S, SolverType solverType, bool verbose=false);

    LinSolver(CSparseMatrix & S, SolverType solverType, bool verbose=false);

    ~LinSolver();

    RVector operator()(const RVector & rhs);
    CVector operator()(const CVector & rhs);

    void solve(const RVector & rhs, RVector & solution);
    void solve(const CVector & rhs, CVector & solution);

    RVector solve(const RVector & rhs);
    CVector solve(const CVector & rhs);

    void setSolverType(SolverType solverType=AUTOMATIC);
    
    // void setSolver(const std::string & name);

    /*! Forwarded to the wrapper to overwrite settings within S. stype =-2 -> use S.stype()*/
    void setMatrix(RSparseMatrix & S, int stype=-2);

    /*! Verbose level = -1, use Linsolver.verbose(). */
    void setMatrix(CSparseMatrix & S, int stype=-2);

    SolverType solverType() const { return solverType_; }

    std::string solverName() const;

    std::string name() const { return solverName(); }

protected:
    void init_();

    void initialize_(RSparseMatrix & S, int stype);
    void initialize_(CSparseMatrix & S, int stype);

    MatrixBase * cacheMatrix_;
    SolverType      solverType_;
    SolverWrapper * solver_;

    uint rows_;
    uint cols_;
};

template < class Mat, class Vec > int solveLU(const Mat & A, Vec & x, const Vec & b){

	//** from TETGEN

    Mat lu(A);

    int N = 0;
    uint n = b.size();
    int ps[n];

    double scales[n];
    double pivot, biggest, mult, tempf;
    uint pivotindex = 0, tmpIdx = 0;

    for (uint i = N; i < n + N; i++) {
        // Find the largest element in each row for row equilibration
        biggest = 0.0;
        for (uint j = N; j < n + N; j++) {
            if (biggest < (tempf = std::fabs(lu[i][j]))) biggest  = tempf;
        }

        if (biggest != 0.0) scales[i] = 1.0 / biggest;
        else {
            scales[i] = 0.0;
            std::cerr << WHERE_AM_I << " Zero row: singular matrix" << std::endl;
            return false;                      // Zero row: singular matrix.
        }
        ps[i] = i;                                 // Initialize pivot sequence.
    }

    for (uint k = N; k < n + N - 1; k++) {      // For each column.
        // Find the largest element in each column to pivot around.
        biggest = 0.0;
        for (uint i = k; i < n + N; i++) {
            if (biggest < (tempf = std::fabs(lu[ps[i]][k]) * scales[ps[i]])) {
                biggest = tempf;
                pivotindex = i;
            }
        }
        if (biggest == 0.0) {
            std::cerr << WHERE_AM_I << " Zero column: singular matrix" << std::endl;
            return false;                         // Zero column: singular matrix.
        }

        if (pivotindex != k) {                         // Update pivot sequence.
            tmpIdx = ps[k];
            ps[k] = ps[pivotindex];
            ps[pivotindex] = tmpIdx;
            //      *d = -(*d);               // ...and change the parity of d.
        }

        // Pivot, eliminating an extra variable  each time
        pivot = lu[ps[k]][k];

        for (uint i = k + 1; i < n + N; i++) {
            lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
            if (mult != 0.0) {
                for (uint j = k + 1; j < n + N; j++) lu[ps[i]][j] -= mult * lu[ps[k]][j];
            }
        }
    }

    double dot = 0.0;

    // Vector reduction using U triangular matrix.
    for (uint i = N; i < n + N; i++) {
        dot = 0.0;
        for (uint j = N; j < i + N; j++) dot += lu[ps[i]][j] * x[j];
        x[i] = b[ps[i]] - dot;
    }

    // Back substitution, in L triangular matrix.
    for (int i = n + N - 1; i >= N; i--) {
        dot = 0.0;
        for (uint j = i + 1; j < n + N; j++) dot += lu[ps[i]][j] * x[j];
        x[i] = (x[i] - dot) / lu[ps[i]][i];
    }

    if (rms(Vec(A * x - b)) > 1e-9){
        std::cerr << "rms(A * x -b) " << rms((const Vec)(A * x - b)) << std::endl;
    //     std::cout << Vector(A * x - b) << std::endl;
    //     std::cout << A << b << x << std::endl;
        return -1;
    }
    return 1;
}


} // namespace GIMLI{

#endif //_GIMLI_LINSOLVER__H
