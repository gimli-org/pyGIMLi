/***************************************************************************
 *   Copyright (C) 2006-2014 by the GIMLi development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _GIMLI_LDLWRAPPER__H
#define _GIMLI_LDLWRAPPER__H

#include "gimli.h"
#include "solverWrapper.h"

namespace GIMLI{

class DLLEXPORT LDLWrapper : public SolverWrapper {
public:
    LDLWrapper(RSparseMatrix & S, bool verbose=false);
    
    LDLWrapper(CSparseMatrix & S, bool verbose=false);
    
    virtual ~LDLWrapper();
 
    static bool valid();

    int factorise();
    
    virtual void setMatrix(RSparseMatrix & S);

    virtual void solve(const RVector & rhs, RVector & solution);
    
protected:
  
    int initialize_(RSparseMatrix & S);
    void free_();

    int * colPtr_;
    int * rowIdx_;
    double * vals_;

    int * Li_;
    int * Lp_;
    double * Lx_;
    double * D_;

    bool preordering_;
    int * P_;
    int * Pinv_;
};

} //namespace GIMLI;

#endif // _GIMLI_LDLWRAPPER__H
