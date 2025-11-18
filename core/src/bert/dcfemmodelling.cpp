/******************************************************************************
 *   Copyright (C) 2006-2024 by the resistivity.net development team          *
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

#include "dcfemmodelling.h"
#include "bertJacobian.h"
#include "bertMisc.h"
#include "bertDataContainer.h"
#include "datamap.h"
#include "electrode.h"

#include <datacontainer.h>
#include <elementmatrix.h>
#include <expressions.h>

#include <interpolate.h>
#include <linSolver.h>
#include <matrix.h>
#include <memwatch.h>
#include <mesh.h>
#include <numericbase.h>

#include <regionManager.h>
#include <shape.h>
#include <sparsematrix.h>
#include <stopwatch.h>
#include <vectortemplates.h>

#undef HAVE_LIBBOOST_THREAD

#ifdef HAVE_LIBBOOST_THREAD
#include <boost/thread.hpp>
#endif

namespace GIMLI{

void setComplexResistivities(Mesh & mesh,
                             const std::map < float, Complex > & aMap){
    std::map< float, Complex >::const_iterator itm;

    RVector am(mesh.cellCount());
    RVector ph(mesh.cellCount());

    if (aMap.size() != 0){
        for (Index i = 0, imax = mesh.cellCount(); i < imax; i++){
            itm = aMap.find(float(mesh.cell(i).marker()));
            if (itm != aMap.end()) {
                am[mesh.cell(i).id()] = std::real((*itm).second);
                ph[mesh.cell(i).id()] = std::imag((*itm).second);
            }
        }
    }
    // Assuming aMap is Ohm , phase(rad)
    setComplexResistivities(mesh, am, ph);
}

void setComplexResistivities(Mesh & mesh,
                             const RVector & am,
                             const RVector & ph){
    setComplexResistivities(mesh, polarToComplex(am, ph, true));
}

void setComplexResistivities(Mesh & mesh, const CVector & z){
    mesh.addData("AttributeReal", real(z));
    mesh.addData("AttributeImag", imag(z));
}

CVector getComplexResistivities(const Mesh & mesh){
    if (!mesh.haveData("AttributeReal") || !mesh.haveData("AttributeImag")){
        throwError(WHERE_AM_I +
                        " complex resistivity values expected but non found");
    }
    RVector re(mesh.data("AttributeReal"));
    RVector im(mesh.data("AttributeImag"));
    return toComplex(re, im);
}

void setComplexData(DataContainer & data,
                    const RVector & re,
                    const RVector & im){
    __MS("setComplexData")
    setComplexData(data, toComplex(re, -im));
}

void setComplexData(DataContainer & data, const CVector & z){
    __MS("setComplexData")
    data.set("u", abs(z));
    data.set("ip", -angle(z) * 1000);
}

CVector getComplexData(const DataContainer & data){
    if (!data.allNonZero("rhoa") || !data.exists("ip")){
        throwError(WHERE_AM_I  + " We need rhoa and ip to get complex data.");
    }
    RVector am(data("rhoa"));
    RVector ph(data("ip"));
    return polarToComplex(am, ph, true);
}

template < class Vec > bool checkIfMapFileExistAndLoadToVector(const std::string & filename, Vec & v){
    bool fromOne = true;

    if (fileExist(filename)){
        std::map < float, float > iMap(loadFloatMap(filename));
        if ((uint)rint(iMap.begin()->first) == 0){
            fromOne = false;
        }
        for (std::map <float, float>::iterator it = iMap.begin(); it != iMap.end(); it ++){

            if (it->first==-1){
                v[v.size() - 1] = it->second;
            } else {
                    uint idx = (uint)rint(it->first);
                    //std::cout << "idx: " << idx << " " << it->second <<" " << fromOne <<std::endl;
                    if (fromOne){
                        if (idx <= v.size() && idx > 0) v[idx - 1] = it->second;
                    } else {
                        if (idx < v.size() && idx >= 0) v[idx] = it->second;
                    }
            }
        }
        return true;
    }
    return false;
}

template < class ValueType >
void assembleStiffnessMatrixHomogenDirichletBC(SparseMatrix < ValueType > & S,
                                               const IndexArray & nodeID,
                                               std::vector < Vector < ValueType > > & rhs){

    for (Index i = 0; i < nodeID.size(); i ++){
        S.cleanRow(nodeID[i]);
        S.cleanCol(nodeID[i]);
        S.setVal(nodeID[i], nodeID[i], 1.0);
        if (rhs.size() == S.rows()){
            for (Index j = 0; j < rhs.size(); j ++) rhs[j][nodeID[i]] = ValueType(0.0);
        }
    }
}

template < class ValueType >
void assembleStiffnessMatrixHomogenDirichletBC(SparseMatrix < ValueType > & S,
                                               const IndexArray & nodeID){
    std::vector < Vector < ValueType > > rhs(0);
    assembleStiffnessMatrixHomogenDirichletBC(S, nodeID, rhs);
}

template < class ValueType >
void dcfemDomainAssembleStiffnessMatrix(SparseMatrix < ValueType > & S, const Mesh & mesh,
                                        const Vector < ValueType > & atts,
                                        double k, bool fix){
    S.clean();
    uint countRho0 = 0, countforcedHomDirichlet = 0;

    if (!S.valid()) S.buildSparsityPattern(mesh);

    ElementMatrix < double > Se, Stmp;

    if (atts.size() != mesh.cellCount()){
       throwLengthError(WHERE_AM_I + " attribute size missmatch" + str(atts.size())
                       + " != " + str(mesh.cellCount()));
    }
    ValueType rho = 0.0;
    Stopwatch swatch(true);
    //unsigned long sCount = 0;

    for (uint i = 0; i < mesh.cellCount(); i++){
        rho = atts[mesh.cell(i).id()];
        //** rho == 0.0 may happen while secondary field assemblation
        if (GIMLI::abs(rho) > TOLERANCE){
            if (k > 0.0){
//             Stmp = Se.u2(mesh.cell(i));
//             Stmp *= k * k;
//             Stmp += Se.ux2uy2uz2(mesh.cell(i));

                Stopwatch s(true);
                Se.u2(mesh.cell(i));
                //sCount += s.cycleCounter().toc();
                Se *= k * k;
                Se += Stmp.ux2uy2uz2(mesh.cell(i));
            } else {
                Se.ux2uy2uz2(mesh.cell(i));
            }
	S.add(Se, 1./rho);
//             Se *= 1.0 / rho;
//             S += Se;
        } else {
            //std::cout << WHERE_AM_I << " " << rho << std::endl;
        }
        if (rho < ValueType(0.0) && fix) countRho0++;
    }

    //std::cout << "assemble time: " << swatch.cycleCounter().toc() << " " << sCount << "  " << swatch.duration()  << std::endl;
    if (fix){
        IndexArray fixSingNodesID;
        for (uint i = 0; i < S.size(); i ++){
            if (::fabs(S.getVal(i, i) < TOLERANCE)) {
                fixSingNodesID.push_back(i);
                countforcedHomDirichlet++;
            }
        }
        assembleStiffnessMatrixHomogenDirichletBC(S, fixSingNodesID);
    }

    if (countRho0){
        std::cout << WHERE_AM_I << " WARNING! " << countRho0
                << " cells with rho <= 0.0 found." << std::endl;
    }
    if (countforcedHomDirichlet++){
        std::cout << WHERE_AM_I << " WARNING! " << countforcedHomDirichlet
                << " nodes forced to homogen Dirichlet to fix singularity of stiffness matrix" << std::endl;
    }
}

void dcfemDomainAssembleStiffnessMatrix(RSparseMatrix & S, const Mesh & mesh,
                                        double k, bool fix){
    dcfemDomainAssembleStiffnessMatrix(S, mesh, mesh.cellAttributes(),
                                       k, fix);
}
void dcfemDomainAssembleStiffnessMatrix(CSparseMatrix & S, const Mesh & mesh,
                                        double k, bool fix){
    CVector res(getComplexResistivities(mesh));
    // if (min(imag(res)) < TOLERANCE){
    //     log(Error, "Check *******************************************");
    // }
    dcfemDomainAssembleStiffnessMatrix(S, mesh, res, k, fix);
}
template < class ValueType >
void dcfemBoundaryAssembleStiffnessMatrix(SparseMatrix < ValueType > & S,
                                          const Mesh & mesh,
                                          const Vector < ValueType > & atts,
                                          const RVector3 & source,
                                          double k){
    ElementMatrix < double > Se;
    std::set < Node * > homDirNodes;
    for (Index i = 0, imax = mesh.boundaryCount(); i < imax; i++){
        int marker = mesh.boundary(i).marker();
        if (marker < 0){
            switch (marker){
            case MARKER_BOUND_HOMOGEN_NEUMANN: break;
            case MARKER_BOUND_MIXED:{

                ValueType rho(0.0);
                Cell * cell = mesh.boundary(i).leftCell();
                if (!cell) cell = mesh.boundary(i).rightCell();
                if (!cell) cell = findCommonCell(mesh.boundary(i).shape().nodes());
                if (cell) {
//                     rho = cell->attribute();
                    rho = atts[cell->id()];
                    //rho = 1.0;
                } else {
                    mesh.exportVTK("FailBC");
                    mesh.save("FailBC");
                    throwError(" no cell found for boundary. can't determine mixed boundary conditions. See FailBC exports." + str(i));
                }

                if (GIMLI::abs(rho) < TOLERANCE){
                    std::cerr << WHERE_AM_I << " parameter rho == 0.0 found " << rho << std::endl;
                }
                Se.u2(mesh.boundary(i));
                S.add(Se, (mixedBoundaryCondition(mesh.boundary(i), source, k) / rho));
                //Se *= (mixedBoundaryCondition(mesh.boundary(i), source, k) / rho);
                // S += Se;
            } break;
            case MARKER_BOUND_HOMOGEN_DIRICHLET:
                for (Index n = 0; n < mesh.boundary(i).nodeCount(); n ++){
                    homDirNodes.insert(&mesh.boundary(i).node(n));
                }
                break;
            case MARKER_BOUND_DIRICHLET: THROW_TO_IMPL; break;
            default:
	       //	std::cerr << WHERE_AM_I << " boundary condition for marker " << marker << " not
                 //defined." << std::endl;
                break;
            }
        }
    }

    IndexArray vecHomDirNodes;
    for (std::set< Node * >::iterator it = homDirNodes.begin(); it != homDirNodes.end(); it ++){
        vecHomDirNodes.push_back((*it)->id());
    }
    assembleStiffnessMatrixHomogenDirichletBC(S, vecHomDirNodes);
}

void dcfemBoundaryAssembleStiffnessMatrix(RSparseMatrix & S, const Mesh & mesh,
                                          const RVector3 & source,
                                          double k){
    dcfemBoundaryAssembleStiffnessMatrix(S, mesh, mesh.cellAttributes(),
                                         source, k);
}

void dcfemBoundaryAssembleStiffnessMatrix(CSparseMatrix & S, const Mesh & mesh,
                                          const RVector3 & source,
                                          double k){
    dcfemBoundaryAssembleStiffnessMatrix(S, mesh, getComplexResistivities(mesh),
                                         source, k);
}

void assembleCompleteElectrodeModel_(RSparseMatrix & S,
                                    const std::vector < ElectrodeShape * > & elecs,
                                    uint oldMatSize, bool lastIsReferenz,
                                    const RVector & contactImpedances){
    RSparseMapMatrix mapS(S);
    ElementMatrix < double > Se;

    uint nElectrodes = elecs.size();
    mapS.setRows(oldMatSize + nElectrodes);
    mapS.setCols(oldMatSize + nElectrodes);

    // RVector vContactImpedance( nElectrodes, 1.0); // Ohm * m^2

    // bool hasImp = checkIfMapFileExistAndLoadToVector("contactImpedance.map",  vContactImpedance);

    bool hasImp = true;
    RVector vContactResistance(nElectrodes, 1.0); // Ohm
    bool hasRes = checkIfMapFileExistAndLoadToVector("contactResistance.map", vContactResistance);

    for (uint elecID = 0; elecID < nElectrodes; elecID ++){

        //** some scale value, can used for contact impedance
        double sumArea = elecs[elecID]->domainSize();
        uint mat_ID = oldMatSize + elecID;
        // __MS(elecID)
        // __MS(sumArea)
        // __MS(elecs[elecID])
        elecs[elecID]->setMID(mat_ID);

        double contactResistance = vContactResistance[elecID];
        double contactImpedance  = contactImpedances[elecID];

        std::vector < MeshEntity * > electrodeEnts(elecs[elecID]->entities());
        if (hasImp || hasRes){
            if (sumArea < TOLERANCE){ //** point electrode
                contactImpedance = 1.0;
                sumArea = 1.0;
            } else {
                if (hasRes) contactImpedance = contactResistance * sumArea;
                else if (hasImp) contactResistance = contactImpedance / sumArea;
            }
            if (sumArea != 1.0){
                std::cout << "Electrode " << elecs[elecID]->id()
                    << " Contact- resistance: "<< contactResistance << " Ohm"
                    << " - impedance: " << contactImpedance  << " Ohm m^2"
                    << " - area: " << sumArea << " m^2" << std::endl;
            }
        }

        //std::cout << "electrode facet contact impedance: " << contactImpedance << std::endl;
        for (uint j = 0; j < electrodeEnts.size(); j ++){

            Se.u(*electrodeEnts[j]);

            if (lastIsReferenz && elecID == nElectrodes - 1){
                Se /= contactImpedance;
            } else {
                Se /= -contactImpedance;
            }

            //** C
            mapS.addToCol(mat_ID, Se);
                //** C'
            mapS.addToRow(mat_ID, Se);
//            std::cout << Se<< std::endl;
            //if (::fabs(contactImpedance - 1.0) > TOLERANCE){
                //** B +=
                Se.u2(*electrodeEnts[j]);
                Se /= contactImpedance;
                mapS += Se;
//             } else {
//                 std::cout << " cem: without contactImpedance " << std::endl;
//             }
        } // for each: electrode entity

        //** G
        if (lastIsReferenz){
            throwError("CEM with lastIsReferenz is currently not supported. Please add Reference Electrode");
            //**!! this leads to nonpositive definite S .. pls check
            if (elecID != nElectrodes- 1){
                std::cout << " cem: last is reference" << std::endl;
                uint refID = oldMatSize + nElectrodes -1;

                mapS[refID][mat_ID] = sumArea / contactImpedance;
                mapS[mat_ID][refID] = sumArea / contactImpedance;

                mapS[refID][refID]   = 2.0 * sumArea / contactImpedance;
                mapS[mat_ID][mat_ID] = 2.0 * sumArea / contactImpedance;
            }
        } else {
            if (::fabs(sumArea) < TOLERANCE){
                //** asuming node electrode
                mapS[mat_ID][mat_ID] = 1.0;
            } else {
                mapS[mat_ID][mat_ID] = sumArea / contactImpedance;
            }
        }
    } // for each: electrode
    S = mapS;
}

void assembleCompleteElectrodeModel(RSparseMatrix & S,
                                    const std::vector < ElectrodeShape * > & elecs,
                                    uint oldMatSize, bool lastIsReferenz,
                                    const RVector & contactImpedances){
    assembleCompleteElectrodeModel_(S, elecs, oldMatSize, lastIsReferenz, contactImpedances);
}

void assembleCompleteElectrodeModel(CSparseMatrix & S,
                                    const std::vector < ElectrodeShape * > & elecs,
                                    uint oldMatSize, bool lastIsReferenz,
                                    const RVector & contactImpedances){
    THROW_TO_IMPL
}

double mixedBoundaryCondition(const Boundary & boundary, const RVector3 & source, double k){
    if (!source.valid()){
        std::cerr << WHERE_AM_I << " no valid source found " << std::endl;
        return 0.0;
    }
    double mirrorPlaneZ = 0.0;
    RVector3 sourceMir(source);
    uint dim = 3;
    if (k > 0) dim = 2;
    sourceMir[dim - 1] = 2.0 * mirrorPlaneZ - source[dim - 1];

    RVector3 facetPos(boundary.center());
    RVector3 norm(boundary.norm());
    RVector3 r(source - facetPos);
    RVector3 rMir(sourceMir - facetPos);
    double rAbs = r.abs(), rMirAbs = rMir.abs();

//   std::cout << " S: " << source << " S': " << sourceMir
// 	    << " F: " << facetPos << " B: " << boundary.node(0) << " " << boundary.node(1) << " n: " << norm
// 	    << " r: "  << r<< " r'" << rMir << std::endl;

    double result = 0.0;
    enum spaceConfig{HALFSPACE,FULLSPACE,MIRRORSOURCE} config = MIRRORSOURCE;

    if (k == 0){ // 3D
        result = ((rMirAbs * rMirAbs) * std::fabs(r.dot(norm))    / rAbs      +
                   (rAbs * rAbs)       * std::fabs(rMir.dot(norm)) / rMirAbs) /
                   (rMirAbs * rAbs * (rAbs + rMirAbs));

                //(|r'|² * |r * n| / |r| + |r|² * |r' * n| / |r'|) / (|r'|*|r|* (|r|+|r'|))
//     switch(config){
//     case HALFSPACE:
//       result =::fabs(r.scalar(norm)) / (rAbs * rAbs);
//       break;
//     case FULLSPACE: TO_IMPL break;
//     case MIRRORSOURCE:
    // nach Bing & Greenhalgh

    //    cout << alpha << std::endl;
//       break;
//     default:
//       std::cerr << WHERE_AM_I << " Warning SpaceConfigEnum = " << config << std::endl;
//       break;
//     }
    } else { // 2.5D
        switch(config){
        case HALFSPACE:
            if (std::fabs(besselK0(rAbs * k)) < 1e-40) return 0.0;
            result = k * std::fabs(r.dot(norm)) / rAbs * besselK1(rAbs * k) / besselK0(rAbs * k);
        break;
        case FULLSPACE:
            if (std::fabs(besselK0(rAbs * k)) < 1e-40) return 0.0;
            //ca: is there missing factor 2??????????
            result = k * std::fabs(r.dot(norm)) / rAbs * besselK1(rAbs * k) / besselK0(rAbs * k);
        break;
        case MIRRORSOURCE:
            if ((::fabs(besselK0(rAbs * k)) < TOLERANCE) ||
                    (::fabs(besselK0(rMirAbs * k)) < TOLERANCE)) return 0.0;

            result = k * (::fabs(r.dot(norm)) / rAbs * besselK1(rAbs * k) +
		      ::fabs(rMir.dot(norm)) / rMirAbs  * besselK1(rMirAbs * k)) /
	       (besselK0(rAbs * k) + besselK0(rMirAbs * k));
        break;
        }
    }

    if (std::isnan(result) || std::isinf(result) || std::fabs(result) < TOLERANCE){
            std::cerr << WHERE_AM_I << " Warning " << result << std::endl;
            std::cerr << "Source: " << source << std::endl;
            std::cerr << "n: " << norm << std::endl;
            std::cerr << "r: " << r << " rMir " << rMir << std::endl;
            std::cerr << "besselK1(rAbs * k) " << besselK1(rAbs * k) << " k " << k << std::endl;
            std::cerr << "rMirAbs " << rMirAbs << " rAbs " << rAbs << std::endl;
    }

    return result;
}

DCMultiElectrodeModelling::DCMultiElectrodeModelling(bool verbose)
    : ModellingBase(verbose) {
    init_();
}

DCMultiElectrodeModelling::DCMultiElectrodeModelling(Mesh & mesh, bool verbose)
    : ModellingBase(verbose) {
    init_();
    setMesh(mesh);
}

DCMultiElectrodeModelling::DCMultiElectrodeModelling(DataContainerERT & dataContainer, bool verbose)
    : ModellingBase(dataContainer, verbose){
    init_();
}

DCMultiElectrodeModelling::DCMultiElectrodeModelling(Mesh & mesh, DataContainerERT & dataContainer, bool verbose)
    : ModellingBase(dataContainer, verbose){
    init_();
    setMesh(mesh);
}

DCMultiElectrodeModelling::~DCMultiElectrodeModelling(){
    if (subSolutions_ && subpotOwner_) {
        delete subSolutions_;
    }

    if (electrodeRef_ && electrodeRef_ != electrodes_.back()){
        delete electrodeRef_;
    }

    if (primDataMap_) delete primDataMap_;

    for_each(electrodes_.begin(), electrodes_.end(), deletePtr());
}

void DCMultiElectrodeModelling::init_(){
    analytical_          = false;
    topography_          = false;
    neumannDomain_       = true;
    lastIsReferenz_      = false;
    complex_             = false;
    setSingValue_        = true;

    subpotOwner_         = false;
    subSolutions_        = NULL;

    electrodeRef_        = NULL;
    JIsRMatrix_          = true;
    JIsCMatrix_          = false;

    solver_              = nullptr;
    buildCompleteElectrodeModel_    = false;
    dipoleCurrentPattern_           = false;

    primDataMap_ = new DataMap();

    byPassFile_ = "bypass.map";

    Index nThreads = getEnvironment("BERTTHREADS", 0, verbose_);
    nThreads = getEnvironment("BERT_NUM_THREADS", 0, verbose_);
    if (nThreads > 0) setThreadCount(nThreads);
    //log(Info, "DCMultiElectrodeModelling init with thread count:", nThreads);
}

void DCMultiElectrodeModelling::setComplex(bool c) {
    if (complex_ != c){

        if (subSolutions_ && subpotOwner_) {
            delete subSolutions_;
            subSolutions_ = 0;
        }
        complex_=c;
    }
}

DataContainerERT & DCMultiElectrodeModelling::dataContainer() const{
    return dynamic_cast < DataContainerERT & >(*dataContainer_);
}

void DCMultiElectrodeModelling::deleteMeshDependency_(){
    for_each(electrodes_.begin(), electrodes_.end(), deletePtr()); electrodes_.clear();
    electrodeRef_        = NULL;
}

void DCMultiElectrodeModelling::assembleStiffnessMatrixDCFEMByPass(RSparseMatrix & S){
    assembleStiffnessMatrixDCFEMByPass_(S);
}

void DCMultiElectrodeModelling::assembleStiffnessMatrixDCFEMByPass(CSparseMatrix & S){
    //assembleStiffnessMatrixDCFEMByPass_(S);
}

template < class ValueType >
void DCMultiElectrodeModelling::assembleStiffnessMatrixDCFEMByPass_(SparseMatrix < ValueType > & _S){

    std::vector < std::pair< Index, Index> > byPassPair;
    std::vector < std::pair< Index, Index> > byPassNodesPair;
    std::vector < double > resistance;
    std::vector < double > resistanceNode;
    std::vector < std::string > row;

    if (fileExist(byPassFile_)){
        if (verbose_) std::cout << byPassFile_ << " found in working path. Applying them." << std::endl;
        std::fstream file; openInFile(byPassFile_, &file, true);

        // Check bypass File
        while(!file.eof()) {
            row  = getNonEmptyRow(file);
            if (row.size() == 3){
                Index eK1idx = toInt(row[0]);
                Index eK2idx = toInt(row[1]);
                double resis = toFloat(row[2]);

                Index a1 = 0, a2 = 0;

                if (eK1idx < 1){
                    std::cerr << WHERE_AM_I << " bypass electrode unknown: " << eK1idx << " please choose electrode indices from 1. Ignoring." << std::endl;
                    continue;
                }
                if (eK1idx < electrodes_.size()) {
                    a1 = electrodes_[eK1idx - 1]->mID();
                } else {
                    std::cerr << WHERE_AM_I << " bypass electrode unknown " << eK1idx << " e-size = "
                            << electrodes_.size() << " Ignoring."<< std::endl;
                    continue;
                }
                if (eK1idx < electrodes_.size()) {
                    a2 = electrodes_[eK2idx - 1]->mID();
                } else {
                    std::cerr << WHERE_AM_I << " bypass electrode unknown " << eK1idx << " e-size = "
                            << electrodes_.size() << " Ignoring."<< std::endl;
                    continue;
                }

                byPassNodesPair.push_back(std::pair< Index, Index>(a1, a2));
                resistanceNode.push_back(resis);
            } else if (row.size() == 2){
                int nodeMarker = toInt(row[0]);
                double resis = toFloat(row[1]);

                IndexArray nodesIDX(mesh_->findNodesIdxByMarker(nodeMarker));
                if (nodesIDX.size()){
                    for (Index i=0; i < nodesIDX.size()-1; i ++ ){
                        byPassNodesPair.push_back(std::pair< Index, Index>(nodesIDX[i], nodesIDX[i+1]));
                        resistanceNode.push_back(resis);
                    }
                } else {
                    std::cerr << WHERE_AM_I << " Warning! cannot requested node marker ("+str(nodeMarker)+") for bypass.map.\n" <<
                    "Expect either: \n int(ElectrodeID) int(ElectrodeID) double(resistance)\n or \n"<<
                    "int(NodeMarker) double(resistance)" << row.size() << std::endl;
                }

            } else if (row.size() > 0){
                std::cerr << WHERE_AM_I << " Warning! Wrong format for bypass.map.\n" <<
                    "Expect either: \n int(ElectrodeID) int(ElectrodeID) double(resistance)\n or \n"<<
                    "int(NodeMarker) double(resistance)" << row.size() << std::endl;
            }
        } // while file
        file.close();
    } // if file can be read

    //## looking for CEM nodes
    if (bypassNodeIdx_.size()){
        std::map < float, float > rMap;
        if (fileExist("electrodeBypassResistances.map")){
            rMap = loadFloatMap("electrodeBypassResistances.map");
        }

        for (Index j = 0; j < bypassNodeIdx_.size(); j ++){
            int marker = bypassNodeIdx_[j];

            IndexArray nodesIDX(mesh_->findNodesIdxByMarker(marker));

            for (Index i = 0; i < nodesIDX.size()-1; i ++){
                byPassNodesPair.push_back(std::pair< Index, Index>(nodesIDX[i],
                                                                   nodesIDX[i+1]));

                if (rMap.count(float(marker))){
                    resistanceNode.push_back(rMap[float(marker)]);
                } else {
                    resistanceNode.push_back(1e-6);
                }
            }
        }
    }

    RSparseMapMatrix S(_S);
    for (Index i = 0; i < byPassNodesPair.size(); i ++){
        Index a1 = byPassNodesPair[i].first;
        Index a2 = byPassNodesPair[i].second;
        if (verbose_) std::cout << "Bypass nodes: " << a1 << "-" << a2 << " with resistance: "
                          << resistanceNode[i] << std::endl;
        double val = 1.0 / resistanceNode[i];

        S[a1][a1] += val;
        S[a2][a2] += val;
        S[a1][a2] -= val;
        S[a2][a1] -= val;
    }
    _S = S;
}

void DCMultiElectrodeModelling::updateMeshDependency_(){
    if (subSolutions_) subSolutions_->clear();

    for_each(electrodes_.begin(), electrodes_.end(), deletePtr());
    electrodes_.clear();

    electrodeRef_        = NULL;
    //** Try to analyse the geometrie, check for topography and looking
    // for surface Z-Koordinate
    topography_     = false;
    neumannDomain_  = true;
    surfaceZ_       = -MAX_DOUBLE;
    bool init       = false;

    for (Index i = 0; i < mesh_->boundaryCount(); i++){
        if (mesh_->boundary(i).marker() == MARKER_BOUND_MIXED ||
            mesh_->boundary(i).marker() == MARKER_BOUND_HOMOGEN_DIRICHLET ||
            mesh_->boundary(i).marker() == MARKER_BOUND_DIRICHLET){
            neumannDomain_ = false;
        }

        if (mesh_->boundary(i).marker() == MARKER_BOUND_HOMOGEN_NEUMANN &&
            !topography_){
            if (init == false) {
                surfaceZ_ = mesh_->boundary(i).center()[mesh_->dim() -1];
                init = true;
            } else {
                if (mesh_->boundary(i).center()[mesh_->dim() -1] != surfaceZ_){
                    if (verbose_) {
                        std::cout << "Found topography for surface="
                                  << surfaceZ_ << " : "
                                  << mesh_->boundary(i).center()[mesh_->dim() -1] << std::endl;
                    }
                    topography_ = true;
                }
            }
        }
    }

    if (neumannDomain_) {
        if (verbose_) {
            std::cout << "Found Neumann domain. Setting topography=1." << std::endl;
        }
        topography_ = true;
    }

    //** when we calculate 2,5D we never have neumannDomain
    if (mesh_->dim() == 2){
        if (neumannDomain_){
            neumannDomain_ = false;
            if (verbose_) {
                std::cout << "Found Neumann domain. but 2.5D -> neumann: false" << std::endl;
            }
        }
    }
//    mesh_->exportBoundaryVTU("meshBound");
    if (mesh_->haveData("AttributeReal") && mesh_->haveData("AttributeImag")){
        this->setComplex(true);
    }

    //## new mesh but old data .. so we need search electrodes again
    searchElectrodes_();
}

void DCMultiElectrodeModelling::updateDataDependency_(){

    if (subSolutions_) subSolutions_->clear();

    for_each(electrodes_.begin(), electrodes_.end(), deletePtr());
    electrodes_.clear();

    electrodeRef_        = NULL;
    if (mesh_) searchElectrodes_();
}

void DCMultiElectrodeModelling::setContactImpedances(const RVector & zi){
    vContactImpedance_ = zi;
}

void DCMultiElectrodeModelling::searchElectrodes_(){
    // __MS(*(this->mesh_))
    if (!mesh_){
        throwError("DCMultiElectrodeModelling::searchElectrodes_() have no mesh defined");
    }
    if (electrodes_.size() > 0) return;

    passiveCEM_.clear();

    //** step 1 search the mesh for signs of electrodes
    std::vector < Index > sourceIdx = mesh_->findNodesIdxByMarker(MARKER_NODE_ELECTRODE);
    IndexArray refSourceIdx  = mesh_->findNodesIdxByMarker(MARKER_NODE_REFERENCEELECTRODE);
    calibrationSourceIdx_    = mesh_->findNodesIdxByMarker(MARKER_NODE_CALIBRATION);

    //** looking for CEM-boundaries
    std::map< int, std::vector< MeshEntity * > > electrodeFaces;

    for (Index i = 0, imax = mesh_->boundaryCount(); i < imax; i++){
        int marker = mesh_->boundary(i).marker();
        if (marker <= MARKER_BOUND_ELECTRODE + 1){ // +1 since -9999 is passive body
            electrodeFaces[MARKER_BOUND_ELECTRODE - marker].push_back(&mesh_->boundary(i));
        }
    }

    //** looking for CEM-cells
    for (Index i = 0, imax = mesh_->cellCount(); i < imax; i++){
        int marker = mesh_->cell(i).marker();

        if (marker <= MARKER_BOUND_ELECTRODE and marker > MARKER_FIXEDVALUE_REGION){
            electrodeFaces[MARKER_BOUND_ELECTRODE - marker].push_back(&mesh_->cell(i));
        }
    }

    //** looking for CEM-nodes (bypass)
    std::map< int, std::vector< Node * > > bypassNodes;
    for (Index i = 0, imax = mesh_->nodeCount(); i < imax; i++){
        int marker = mesh_->node(i).marker();
        // no nails. either CEM faces OR CEM bypass nodes
        if (marker <= MARKER_BOUND_ELECTRODE &&
            !electrodeFaces.count(MARKER_BOUND_ELECTRODE - marker)){

            if (!bypassNodes.count(MARKER_BOUND_ELECTRODE - marker)) {
                bypassNodeIdx_.push_back(marker);
            }
            bypassNodes[MARKER_BOUND_ELECTRODE - marker].push_back(&mesh_->node(i));
        }
    }

    std::list < ElectrodeShape * > cemElectrodes;
    std::list < ElectrodeShape * > passiveCEMBodies;
    for (std::map< int, std::vector< MeshEntity * > >::iterator it = electrodeFaces.begin();
         it != electrodeFaces.end(); it ++){
        if (it->first >= 0){
            cemElectrodes.push_back(new ElectrodeShapeDomain(it->second));
            cemElectrodes.back()->setId(it->first);
        } else {
            passiveCEMBodies.push_back(new ElectrodeShapeDomain(it->second));
        }
                //** transform to start with electrode-id == 0;
    }
    Index nodeECounter = 0, nodeBCounter = 0, freeECounter = 0;
    Index cemECounter = 0, passiveCEMbodyCounter = 0;

    if (dataContainer_){
        //!!** match the known electrode to list of electrodes from the datacontainer
        R3Vector ePos(dataContainer_->sensorPositions());

        //!!** For 2d-problems we have to check if xy or xz coordinates are given
        //!!** only xy is valid so it is necessary to swap the coordinates
        if (mesh_->dim() == 2){
            if ((zVari(ePos) || max(abs(z(ePos))) > 0) &&
                (!yVari(ePos) && max(abs(y(ePos))) < 1e-8)){

                if (verbose_) std::cout << "Warning! swap YZ coordinates for sensor positions to meet mesh dimensions." << std::endl;

                swapYZ(ePos);
            }
        }

        for (uint i = 0; i < ePos.size(); i ++){
            bool match = false;
            //** match the known CEM-electrodes
            for (std::list< ElectrodeShape * >::iterator it = cemElectrodes.begin();
                 it != cemElectrodes.end(); it ++){
                    std::cout << "req.pos:" << ePos[i] << " epos("
                                << (*it)->id() << "): " << (*it)->pos()
                                << " dist: " <<  ePos[i].distance((*it)->pos())
                                << " e-size: " << std::sqrt((*it)->domainSize())*2.0
                                << std::endl;

                if (ePos[i].distance((*it)->pos()) < std::sqrt((*it)->domainSize()) * 2.0 ||
                    (uint)(*it)->id() == i) {
                    electrodes_.push_back((*it));
                    electrodes_.back()->setId(i);
                    electrodes_.back()->setPos(ePos[i]);
                    cemElectrodes.erase(it);
                    cemECounter++;
                    match = true;
                    break;
                }
            }
            //** match the known bypass-electrodes
            if (!match){
                for (std::map< int, std::vector< Node * > >::iterator it = bypassNodes.begin();
                     it != bypassNodes.end(); it ++){
// //                     std::cout << ePos[i] << " " << mesh_->node(*it).pos() <<
// //                             " " << ePos[i].dist(mesh_->node(*it).pos()) << std::endl;
                     if (ePos[i].dist(it->second[0]->pos()) < 0.01){ //CR 1cm?? really??
                        electrodes_.push_back(new ElectrodeShapeNodesWithBypass(it->second));
                        electrodes_.back()->setId(i);
                        bypassNodes.erase(it);
                        nodeBCounter++;
                        match = true;
                        break;
                    }
                }
            }

            //** match the known node-electrodes
            if (!match){
                for (std::vector< Index >::iterator it = sourceIdx.begin(); it != sourceIdx.end(); it ++){
                    //  std::cout << ePos[i] << " " << mesh_->node(*it).pos() <<
                    //                          " " << ePos[i].dist(mesh_->node(*it).pos()) << std::endl;
                    if (ePos[i].dist(mesh_->node(*it).pos()) < 0.01){ //CR 1cm?? really??
                        electrodes_.push_back(new ElectrodeShapeNode(mesh_->node(*it)));
                        electrodes_.back()->setId(i);
                        sourceIdx.erase(it);
                        nodeECounter++;
                        match = true;
                        break;
                    }
                }
            }

            //** fill the missing with node independent electrodes
            if (!match){
                Cell * cell = mesh_->findCell(ePos[i]);
                if (cell){
                    electrodes_.push_back(new ElectrodeShapeEntity(*cell, ePos[i]));
                } else{
                    electrodes_.push_back(new ElectrodeShape(ePos[i]));
                    std::cerr << WHERE_AM_I << " " << ePos[i] << " mesh " << mesh_->boundingBox() << std::endl;

                    for (auto &n : sourceIdx){
                        std::cout << ePos[i] << " " << mesh_->node(n).pos() <<
                                                " " << ePos[i].dist(mesh_->node(n).pos()) << std::endl;
                    }
                    throwError("There is a requested electrode that does not match the given mesh. ");
                }

//                     std::cout << ePos[i] << std::endl;
//                     std::cout << *cell << std::endl;

                electrodes_.back()->setId(i);
                freeECounter++;
            }
        } // for all in ePos

    } else { //** no dataContainer_

        //** add all remaining CEM-electrodes
        for (std::list< ElectrodeShape * >::iterator it = cemElectrodes.begin();
            it != cemElectrodes.end(); it ++){
            electrodes_.push_back((*it));
                    //** transform to start with electrode-id == 0;
            electrodes_.back()->setId(cemECounter);
            cemECounter++;
        }

        //** add all known bypass-electrodes
        for (std::map< int, std::vector< Node * > >::iterator it = bypassNodes.begin();
            it != bypassNodes.end(); it ++){
// //                     std::cout << ePos[i] << " " << mesh_->node(*it).pos() <<
// //                             " " << ePos[i].dist(mesh_->node(*it).pos()) << std::endl;
            electrodes_.push_back(new ElectrodeShapeNodesWithBypass(it->second));
            electrodes_.back()->setId(cemECounter + nodeBCounter);
            nodeBCounter++;
        }

        //** add all remaining electrodes to the calculation
        for (std::vector < Index >::iterator it = sourceIdx.begin(); it != sourceIdx.end(); it ++){
            electrodes_.push_back(new ElectrodeShapeNode(mesh_->node(*it)));
            electrodes_.back()->setId(nodeECounter + nodeBCounter + cemECounter);
            nodeECounter++;
        }
    }

    //** add passive cem bodies
    for (std::list< ElectrodeShape * >::iterator it = passiveCEMBodies.begin();
         it != passiveCEMBodies.end(); it ++){
        passiveCEM_.push_back((*it));
        passiveCEM_.back()->setId(-1);
        passiveCEMbodyCounter ++;
    }

    if (cemECounter > 0 || passiveCEMbodyCounter > 0) buildCompleteElectrodeModel_ = true;

    if (verbose_) {
        if (dataContainer_){
            std::cout << "Found datafile: " << dataContainer_->sensorCount()
                    << " electrodes" << std::endl;
        }
        if (cemECounter) {
            std::cout << "Found: " << cemECounter << " cem-electrodes" << std::endl;
        }
        if (nodeBCounter) {
            std::cout << "Found: " << nodeBCounter << " bypass-electrodes" << std::endl;
        }
        if (nodeECounter) {
            std::cout << "Found: " << nodeECounter << " node-electrodes" << std::endl;
        }
        if (freeECounter) {
            std::cout << "Found: " << freeECounter << " free-electrodes" << std::endl;
        }
        if (passiveCEMbodyCounter) {
            std::cout << "Found: " << passiveCEMbodyCounter << " passive cem bodies" << std::endl;
        }
    }

    //** Two special cases:
    //** Just one reference electrode node is possible,
    if (refSourceIdx.size()) {
        if (verbose_) std::cout   << "Found: " << refSourceIdx.size()
                                    << " reference electrode node." << std::endl;
        electrodeRef_ = new ElectrodeShapeNode(mesh_->node(refSourceIdx[0]));
    }

    //** though several calibration point nodes are.
    if (calibrationSourceIdx_.size()) {
        if (verbose_) std::cout << "Found: " << calibrationSourceIdx_.size()
                                    << " calibration node." << std::endl;
    }

    if (electrodes_.size() > 0){

        //** looking for source center Position;
        RVector3 sumPos(0.0, 0.0, 0.0);
        for (uint i = 0; i < electrodes_.size(); i ++){
            if (electrodes_[i]->valid()) sumPos += electrodes_[i]->pos();
        }
        sourceCenterPos_ = sumPos / (double)electrodes_.size();

        if (kValues_.empty() || weights_.empty()){
            if (dataContainer_){
                initKWaveList(*mesh_, kValues_, weights_, dataContainer_->sensorPositions(), verbose_);
            } else {
                initKWaveList(*mesh_, kValues_, weights_, verbose_);
            }
        }
    } else {
        //throwError(WHERE_AM_I+ " Warning ! Found neighter electrode nodes nor electrode facets, don't know what to do. ");
        std::cout << "Warning! Found neigther electrode nodes nor electrode facets, don't know what to do. " << std::endl;
    }

    if (neumannDomain_){
        if (calibrationSourceIdx_.size() == 0) {
            std::cout << "Warning! Neumann domain without calibration (potential reference) point. "
                        << "This may lead to a non-positive definite matrix. LDL can solve instead of"
                            " CHOLMOD, but better you add VIP with calibration marker "
                        << MARKER_NODE_CALIBRATION << std::endl;
            std::cout << "Choose first node as calibration node " << std::endl;
            calibrationSourceIdx_.push_back(0);
        }

        if (!electrodeRef_ && !dipoleCurrentPattern_) {
            if (verbose_) {
                std::cout << "Found Neumann domain without reference electrode. " << std::endl
                            << "Choose last electrode as reference. " << std::endl;
            }
            if (electrodes_.size() == 0){
                log(Warning, "No electrodes defined.");
                return;
            }
            electrodeRef_ = electrodes_.back();
            lastIsReferenz_ = true;
        }
    } else { // no neumann
        if (verbose_) std::cout << "Found non-Neumann domain" << std::endl;
        if (calibrationSourceIdx_.size()){
            std::cout << "Non-Neumann domain conflicts with given calibration point. Ignoring calibration." << std::endl;
            calibrationSourceIdx_.clear();
        }
    }
}

RVector DCMultiElectrodeModelling::createDefaultStartModel(){
    RVector vec(this->regionManager().parameterCount(), 0.0);
    if (dataContainer_ != NULL){
         vec.fill(median(dataContainer_->get("rhoa")));
    } else {
        std::cerr << WHERE_AM_I << " No data container given. " << std::endl;
    }
    return vec;
}

RVector DCMultiElectrodeModelling::response(const RVector & model,
                                            double background){

    if (min(abs(dataContainer_->get("k"))) < TOLERANCE){
        if (!(this->topography() || buildCompleteElectrodeModel_)){
            dataContainer_->set("k",
                              this->calcGeometricFactor(this->dataContainer()));
            log(Warning, " data contains no K-factors but we calculate them "
                         " analytically for the response call");

        } else {
            throwError(WHERE_AM_I + " data contains no K-factors ");
        }
    }

    if (!this->mesh_){
        log(Critical, "Found no mesh, so cannot calculate a response.");
    }
    if (complex()){

        if (min(abs(model) < TOLERANCE)){
            model.save("modelFail.vector");
            log(Critical, " complex response for abs model with negative or zero resistivity is not defined.");
        }

        DataMap dMap(response_(toComplex(model(0, model.size()/2),
                                         model(model.size()/2, model.size())),
                               Complex(background, -9e99)));

        RVector respRe(dMap.data(this->dataContainer(), false, false));
        RVector respIm(dMap.data(this->dataContainer(), false, true));

        CVector resp(toComplex(respRe, respIm) * dataContainer_->get("k"));
        return cat(real(resp), imag(resp));

    } // if complex

    //** to following can lead to problematic situations https://gitlab.com/resistivity-net/bert/issues/41
    // but costs one forward calculation so its probably better to
    // remove it (temporary comment)
    // if (::fabs(max(model) - min(model)) < TOLERANCE){
    //     if (!this->topography_ ){
    //         return RVector(dataContainer_->size(), min(model));
    //     } else {
    //         std::cout << "Calcuating numerical response for homogeneous model due to topography";
    //     }
    // }

    if (min(model) < TOLERANCE){
        model.save("modelFail.vector");
        log(Critical, " response for model with negative or zero resistivity is not defined.:",
                     min(model), max(model));
    }

    DataMap dMap(response_(model, background));
    RVector resp(round(dMap.data(this->dataContainer()), 1e-10));
    RVector respRez(round(dMap.data(this->dataContainer(), true), 1e-10));

    if (resp.size() != dataContainer_->size() || respRez.size() != dataContainer_->size()){
        throwError(WHERE_AM_I + " size wrong: " + str(dataContainer_->size())
        + " " + str(resp.size()) + " " + str(respRez.size()));
    }

    resp    *= dataContainer_->get("k");
    respRez *= dataContainer_->get("k");

    RVector modelReciprocity((resp - respRez) / (resp + respRez) * 2.0);

    if (verbose_){
        if (min(resp) < 0 && 1){
            std::cout << "Found neg. resp (saving)." << std::endl;
                for (uint i = 0; i < resp.size(); i ++){
                    if (resp[i ] < 0) {
                        int a = (*dataContainer_)("a")[i];
                        int b = (*dataContainer_)("b")[i];
                        int m = (*dataContainer_)("m")[i];
                        int n = (*dataContainer_)("n")[i];

                        RVector ab(mesh_->nodeCount(), 0.0), mn(mesh_->nodeCount(), 0.0);
                        if (a != -1) ab = solutions_[a];
                        if (b != -1) ab -= solutions_[b];
                        if (m != -1) mn = solutions_[m];
                        if (n != -1) mn -= solutions_[n];
                        std::cout << i << " " << resp[i] << " " << respRez[i]<< std::endl;
                        std::cout << a << " " << b << " " << m << " " << n << std::endl;

                        mesh_->addData("ab-pot", prepExportPotentialData(ab));
                        mesh_->addData("mn-pot", prepExportPotentialData(mn));
                        //mesh_->addData("sens-mn-pot", prepExportSensitivityData(jacobian));
                        mesh_->exportVTK("negResp");

                        break;
                        //std::cout << (*dataContainer_)[i] << std::endl;
                    }
                }
//std::cout << find(resp < 0) << std::endl;
                mesh_->save("negResp");
                save(mesh_->cellAttributes(), "negResp-Atts");

                save(resp, "resp.vec");
                save(respRez, "respRez.vec");
            } // if found neg. Responses
            std::cout << "Response: min = " << min(resp)
                        << " max = " << max(resp) << " mean = " << mean(resp) << std::endl;
            std::cout << "Reciprocity rms(modelReciprocity) "
                        << rms(modelReciprocity) * 100.0 << "%, "
                        << "max: " << max(modelReciprocity) * 100.0 << "%" << std::endl;
//        std::cout << "Reciprocity sum = " << sum(modelReciprocity) << std::endl;

 //       std::cout << 13 << " " << resp[13] << " " << respRez[13] << std::endl;
    }
    return sqrt(abs(resp * respRez));
}

void DCMultiElectrodeModelling::mapERTModel(const CVector & model, Complex background){
    if (model.size() == this->mesh_->cellCount()){
        setComplexResistivities(*mesh_, model);
    } else {
        RVector re(createMappedModel(real(model), background.real()));
        // RVector im(re*0.0);
        // log(Warning, "imag part forced to zero");
        RVector im(createMappedModel(imag(model), -9e99));
        setComplexResistivities(*mesh_, toComplex(re, im));
    }
}

void DCMultiElectrodeModelling::mapERTModel(const RVector & model,
                                            double background){
    if (model.size() == this->mesh_->cellCount()){
        this->mesh_->setCellAttributes(model);
    } else {
        mapModel(model, background);
    }
}

template < class ValueType >
DataMap DCMultiElectrodeModelling::response_(const Vector < ValueType > & model,
                                             ValueType background){
    if (verbose_) std::cout << "Calculating response for model: min = " << min(model)
                            << " max = " << max(model) << std::endl;
    DataMap dMap;

    this->mapERTModel(model, background);

    if (dataContainer_ != NULL){

        if (dipoleCurrentPattern_){
            THROW_TO_IMPL
            calculate(this->dataContainer(), false);
//             resp    = dataContainer_->get("u");
            //*** No reciprocity for dipole-current pattern.
        } else {
            calculate(dMap);
        }
    } else {
        throwError(WHERE_AM_I + " no response without data container");
    }
    return dMap;
}

template < class ValueType >
Matrix < ValueType > * DCMultiElectrodeModelling::prepareJacobianT_(const Vector< ValueType > & model){
    this->searchElectrodes_();
    if (dataContainer_){
        if (!subSolutions_){
            if (verbose_) {
                std::cout << "Creating new subpotentials for createJacobian."  << std::endl;
            }
            subpotOwner_ = true;
            subSolutions_ = new Matrix< ValueType >;
        }  else {
            if (verbose_) {
                std::cout << "Using existing subpotentials for createJacobian."  << std::endl;
            }
        }
        Matrix< ValueType > *u = NULL;
        u = dynamic_cast< Matrix< ValueType > * >(subSolutions_);
        if (u->rows() == 0){
            if (verbose_) {
                std::cout << "Subpotentials matrix is empty."  << std::endl;
            }
// //            std::cout << WHERE_AM_I << " " << mean(model) << " " << model.size() << std::endl;
//__MS(toc__)
            this->mapERTModel(model, ValueType(-9e99)); // very slow
//__MS(toc__)
            bool oldAna = this->analytical();

            this->setAnalytical(!(this->topography() ||
                                  buildCompleteElectrodeModel_ ||
                                  stdDev(model) > TOLERANCE*1e5));
            if (verbose_) {
                std::cout << "Calculating subpotentials analytical for createJacobian: "
                          << this->analytical() << " ("
                          << "top: " << this->topography() << "|"
                          << "cem: " << buildCompleteElectrodeModel_ << "|"
                          << "het: " << bool(stdDev(model) > TOLERANCE*1e5) << ")" << std::endl;
            }
//__MS(toc__)
            //** first geometric factors .. since this will overwrite u
            if (!dataContainer_->allNonZero("k")){
                dataContainer_->set("k",
                              this->calcGeometricFactor(this->dataContainer(),
                                                        model.size()));
            }

            DataContainerERT tmp(this->dataContainer());

            this->calculate(tmp);

            /*! We have to scale subSolutions_ for the analytical solution to match the model */
            if (this->analytical()){
                if (verbose_) std::cout << "Scale subpotentials with " << model[0] << std::endl;

                for (uint i = 0, imax = u->rows(); i < imax; i ++) {
                    (*u)[i] *= model[0];
                }
            }
            this->setAnalytical(oldAna);
        } // if u.rows()
        return u;
    } else {
        throwError(WHERE_AM_I + " no data structure given");
    }
    return 0;
}

RMatrix * DCMultiElectrodeModelling::prepareJacobian_(const RVector & model){
    return prepareJacobianT_(model);
}
CMatrix * DCMultiElectrodeModelling::prepareJacobian_(const CVector & model){
    return prepareJacobianT_(model);
}

void DCMultiElectrodeModelling::createJacobian_(const RVector & model,
                                                const RMatrix & u, RMatrix * J){

    std::vector < std::pair < Index, Index > > matrixClusterIds;

// MEMINFO
//         save(*u, "pots.bmat");

    createSensitivityCol(*J, *mesh_, this->dataContainer(), u,
                         weights_, kValues_,
                         matrixClusterIds, nThreads_, verbose_);

// MEMINFO
    double sensMatDropTol = getEnvironment("BERT_SENSMATDROPTOL", 0.0, verbose_);

    RSparseMapMatrix * Jsparse = 0;

    if (matrixClusterIds.size() > 0){
        // just test clustering here
        Index nData = matrixClusterIds[0].first;
        Index nModel = matrixClusterIds[0].second;

        if (sensMatDropTol > 0.0){
            // breaks possible blockmatrix sizes
            if (complex_) {THROW_TO_IMPL
            }

            delete jacobian_;
// MEMINFO
            jacobian_ = new RSparseMapMatrix(nData, nModel);
            Jsparse = dynamic_cast< RSparseMapMatrix  * >(jacobian_);
        } else {
            J->resize(nData, nModel);
        }


        for (uint c = 1; c < matrixClusterIds.size(); c ++){
// MEMINFO
            Index start = matrixClusterIds[c].first;
            Index end = matrixClusterIds[c].second;

            if (Jsparse){
                Jsparse->importCol("sensPart_" + str(start) + "-"
                                   + str(end) + ".bmat",
                                   sensMatDropTol, start);
            } else { // no drop tol
                RMatrix Jcluster("sensPart_" + str(start) + "-" + str(end));

                for (Index i = 0; i < J->rows(); i ++){
                    (*J)[i].setVal(Jcluster[i], start, end);
                }
            }
// MEMINFO
        } // for each clustering
    } // if clustering
// MEMINFO

    if (!Jsparse){

        if (model.size() == J->cols()){
            RVector m2(model*model);
            if (model.size() == J->cols()){
                for (uint i = 0; i < J->rows(); i ++) {
                    (*J)[i] /= (m2 / dataContainer_->get("k")[i]);
                }
            }
        }
        if (verbose_){
            RVector sumsens(J->rows());
            for (Index i = 0, imax = J->rows(); i < imax; i ++){
                sumsens[i] = sum((*J)[i]);
            }

            std::cout << "sens sum: median = " << median(sumsens)
                          << " min = " << min(sumsens)
                          << " max = " << max(sumsens) << std::endl;
        }
    } else {
        for (RSparseMapMatrix::iterator it = Jsparse->begin();
             it != Jsparse->end(); it ++){

            Index row = (*it).first.first;
            Index col = (*it).first.second;
            (*it).second /= (model[col] * model[col]) / dataContainer_->get("k")[row];
        }
        if (verbose_){
            std::cout << "S sparsity: " << ((double)Jsparse->nVals() /
            (Jsparse->rows() * Jsparse->cols())) * 100.0 << "% " << std::endl;
        }
    }
}

void DCMultiElectrodeModelling::createJacobian_(const CVector & model,
                                                const CMatrix & u, CMatrix * J){

    std::vector < std::pair < Index, Index > > matrixClusterIds;

    createSensitivityCol(*J,
                         *this->mesh_, this->dataContainer(),
                         u,
                         this->weights_, this->kValues_,
                         matrixClusterIds, this->nThreads_, this->verbose_);

    if (model.size() == J->cols()){
        // __MS("check")
        CVector m2(model*model);
        //CVector m2(model*conj(model));
        if (model.size() == J->cols()){
            for (Index i = 0; i < J->rows(); i ++) {
                (*J)[i] /= (m2 / dataContainer_->get("k")[i]);
            }
        }
    } else  {
        __M
        log(Error, "size mismatch");
    }
    if (verbose_){
        CVector sumsens(J->rows());
        for (Index i = 0, imax = J->rows(); i < imax; i ++){
            sumsens[i] = sum((*J)[i]);
        }

        // std::cout << "sens sum: median = " << median(sumsens)
        //           << " min = " << min(sumsens)
        //           << " max = " << max(sumsens) << std::endl;
    }
}

void DCMultiElectrodeModelling::createJacobian(const RVector & model){
    if (complex_){
        CVector cMod(toComplex(model(0, model.size()/2),
                               model(model.size()/2, model.size())));

        CMatrix * u = this->prepareJacobianT_(cMod);

        if (!JIsCMatrix_){
            // log(Warning, "delete non complex Jacobian and create a new CMatrix");
            delete jacobian_;
            jacobian_ = new CMatrix();
            JIsCMatrix_ = true;
            JIsRMatrix_ = false;

        }
        CMatrix * J = dynamic_cast< CMatrix * >(jacobian_);
        this->createJacobian_(cMod, *u, J);
    } else {
        RMatrix * u = this->prepareJacobianT_(model);
        if (!JIsRMatrix_){
            log(Warning, "delete non real Jacobian and create a new RMatrix");
            delete jacobian_;
            jacobian_ = new RMatrix();
            JIsRMatrix_ = true;
            JIsCMatrix_ = false;
        }

        RMatrix * J = dynamic_cast< RMatrix * >(jacobian_);
        this->createJacobian_(model, *u, J);
    }
}

void DCMultiElectrodeModelling::createConstraints(){
    ModellingBase::createConstraints();
}

void DCMultiElectrodeModelling::createCurrentPattern(std::vector < ElectrodeShape * > & eA,
                                                     std::vector < ElectrodeShape * > & eB,
                                                     bool reciprocity){
    this->searchElectrodes_();
    int nElecs = (int)electrodes_.size();

    if (dipoleCurrentPattern_){

        //** this is only useful for the forward calculation since the reciprocity potentials are needed for sensitivity calculation.
        if (dataContainer_){
            //** reciprocity disabled
            std::set < Index > inject(this->dataContainer().currentPattern(false));
            if (verbose_) std::cout << "Found " << inject.size()
                                        << " dipole-current pattern" << std::endl;
            eA.resize(inject.size(), NULL);
            eB.resize(inject.size(), NULL);
            CurrentPattern cp;
            uint i = 0;

            for (std::set < Index >::iterator it = inject.begin(); it != inject.end(); it ++, i++){
                currentPatternIdxMap_[(*it)] = i;
                cp = this->dataContainer().currentPatternToElectrode((*it));
                if (cp.first < nElecs && cp.second < nElecs){
                    //std::cout << cp.first << " " << cp.second << std::endl;
                    if (cp.first > -1)  eA[i] = electrodes_[cp.first];
                    if (cp.second > -1) eB[i] = electrodes_[cp.second];
                } else {
                    std::cerr << WHERE_AM_I << " requested electrode a = " << cp.first
                                             << " or b = " << cp.second << " do not exist." << std::endl;
                }
            }
        } else {
            std::cerr << WHERE_AM_I << " no data structure given" << std::endl;
        }
    } else { // simple pol or dipole with reference node

        for (int i = 0; i < nElecs; i ++){
            if (electrodes_[i] != electrodeRef_ && electrodes_[i]->id() > -1){
                eA.push_back(electrodes_[i]);
                eB.push_back(electrodeRef_);
            }
        }
    }
}

RVector DCMultiElectrodeModelling::calcGeometricFactor(const DataContainerERT & data,
                                                       Index nModel){
    if (verbose_) std::cout << "Obtaining geometric factors";
    if (!this->topography() && !buildCompleteElectrodeModel_) {
        if (verbose_) std::cout << " (analytical)" << std::endl;
        return geometricFactors(data, mesh_->dimension(), false);
    }

    if (electrodes_.size() == 0){
        this->searchElectrodes_();
    }

    if (primDataMap_->electrodes().size() != electrodes_.size()){
        if (verbose_) std::cout << " (numerical)" << std::endl;
        RVector atts(mesh_->cellAttributes());
        if (nModel > 0) {
            this->mapERTModel(RVector(nModel, 1.0), -9e99);
        } else {
            mesh_->setCellAttributes(RVector(mesh_->cellCount(), 1.0));
        }

        this->calculate(*primDataMap_);

        mesh_->setCellAttributes(atts);
    } else {
        if (verbose_) std::cout << " (recover)" << std::endl;
        THROW_TO_IMPL
    }

    return 1.0 / (primDataMap_->data(data) + TOLERANCE);
}

void DCMultiElectrodeModelling::calculate(DataContainerERT & data, bool reciprocity){

    if (dipoleCurrentPattern_){
        if (complex_) THROW_TO_IMPL

        std::vector < ElectrodeShape * > eA, eB;
        //** no reciprocity while using dipole current pattern
        createCurrentPattern(eA, eB, false);
        calculate(eA, eB);

        if (buildCompleteElectrodeModel_ && potentialsCEM_.rows()){
            std::cout << "Save cemMatrix.matrix for debugging purposes" << std::endl;
            saveMatrixRow(potentialsCEM_, "cemMatrix.matrix");
        }

        RVector u(data.size());
        uint currentIdx = 0;

        for (uint i = 0; i < data.size(); i ++){
            long abPattern = data.electrodeToCurrentPattern(data("a")[i],
                                                            data("b")[i]);

            if (currentPatternIdxMap_.count(abPattern)){
                currentIdx = currentPatternIdxMap_[abPattern];
            } else {
                std::cerr << WHERE_AM_I
                          << " cannot find pattern index. " << std::endl;
            }

            if (data("m")[i] > -1) {
                if (buildCompleteElectrodeModel_){
                    u[i] = potentialsCEM_[currentIdx][data("m")[i]];
                } else {
                    u[i] = electrodes_[data("m")[i]]->pot(solutions_[currentIdx]);
                }
            }
            if (data("n")[i] > -1) {
                if (buildCompleteElectrodeModel_){
                    u[i] -= potentialsCEM_[currentIdx][data("n")[i]];
                } else {
                    u[i] -= electrodes_[data("n")[i]]->pot(solutions_[currentIdx]);
                }
            }

//             if (reciprocity){
//                 if (currentPatternIdxMap_.count(mnPattern)){
//                     currentIdx = currentPatternIdxMap_[mnPattern];
//                 } else {
//                    std::cerr << WHERE_AM_I << " cannot find pattern index. " << std::endl;
//                 }
//
//                 uAB_MN = 0.0;
//                 if (data(i).a() > -1) uAB_MN  = potentialsCEM_[currentIdx][data(i).a()];
//                 if (data(i).b() > -1) uAB_MN -= potentialsCEM_[currentIdx][data(i).b()];
//                 ur[i] = uAB_MN;
//             }
        } //** for each in data

        data.set("u", u);
//         if (reciprocity) data.add("urez", ur);
    } else {
        DataMap dMap;
        this->calculate(dMap);
        if (complex_) {
            log(Error, " DCMultiElectrodeModelling::calculate, don't use this.");
            setComplexData(data, dMap.data(data), dMap.data(data, false, true));
        } else {
            data.set("u", dMap.data(data));
        }
    }
}

void DCMultiElectrodeModelling::calculate(DataMap & dMap){
    //! create current pattern;

    if (dipoleCurrentPattern_){
        throwError(WHERE_AM_I + " Unable to calculate(datamap) using dipoleCurrentPattern. Use calculate(DataContainer) instead ");
    }
    std::vector < ElectrodeShape * > eA, eB;

    createCurrentPattern(eA, eB, true);
    calculate(eA, eB);

    if (buildCompleteElectrodeModel_ && potentialsCEM_.rows()) {
        if (verbose_) std::cout << "Building collectmatrix from CEM matrix appendix." << std::endl;
        dMap.collect(electrodes_, potentialsCEM_, buildCompleteElectrodeModel_);
    } else {
        dMap.collect(electrodes_, solutions_);
    }
}

class CalculateMT{
    public:
    CalculateMT(DCMultiElectrodeModelling * fop, const std::vector < ElectrodeShape * > & eA,
                 const std::vector < ElectrodeShape * > & eB,
                 uint kIdx, RMatrix & mat)
    : fop_(fop), eA_(&eA), eB_(&eB), kIdx_(kIdx), mat_(& mat) {
    }
    void operator()(){
        fop_->calculateK(*eA_, *eB_, *mat_, kIdx_);
    }

    DCMultiElectrodeModelling * fop_;

    const std::vector < ElectrodeShape * > * eA_;
    const std::vector < ElectrodeShape * > * eB_;
    uint kIdx_;
    RMatrix * mat_;
};

void DCMultiElectrodeModelling::calculate(const std::vector < ElectrodeShape * > & eA,
                                          const std::vector < ElectrodeShape * > & eB){

    if (!subSolutions_) {
        subpotOwner_ = true;
        if (complex_){
            subSolutions_ = new CMatrix(0, 0);
        } else {
            subSolutions_ = new RMatrix(0, 0);
        }
    }
    uint nCurrentPattern = eA.size();
    subSolutions_->resize(nCurrentPattern * kValues_.size(),
                          mesh_->nodeCount());

    solutions_.clear();
    if (complex_){
        solutions_.resize(2 * nCurrentPattern, mesh_->nodeCount());
    } else {
        solutions_.resize(nCurrentPattern, mesh_->nodeCount());
    }
    // MEMINFO

    Stopwatch swatch(true);

    // create or find primary potentials
    preCalculate(eA, eB);

    for (Index kIdx = 0; kIdx < kValues_.size(); kIdx ++){
        //if (verbose_ && kValues_.size() > 1) std::cout << "\r" << kIdx + 1 << "/" << kValues_.size();

        if (complex_){
            calculateK(eA, eB, dynamic_cast< CMatrix & > (*subSolutions_), kIdx);
        } else {
            calculateK(eA, eB, dynamic_cast< RMatrix & > (*subSolutions_), kIdx);
        }
    }

    for (Index kIdx = 0; kIdx < kValues_.size(); kIdx ++){
        for (Index i = 0; i < nCurrentPattern; i ++) {
            if (kIdx == 0) {
                if (complex_) {
                    RVector re(real((dynamic_cast< CMatrix & > (*subSolutions_))[i + kIdx * nCurrentPattern]));
                    RVector im(imag((dynamic_cast< CMatrix & > (*subSolutions_))[i + kIdx * nCurrentPattern]));
                    solutions_[i] = re * weights_[kIdx];
                    solutions_[i+ nCurrentPattern] = im * weights_[kIdx];
                } else {
                    solutions_[i] = (dynamic_cast< RMatrix & > (*subSolutions_))[i + kIdx * nCurrentPattern] * weights_[kIdx];
                }
            } else {
                if (complex_){
                    RVector re(real((dynamic_cast< CMatrix & > (*subSolutions_))[i + kIdx * nCurrentPattern]));
                    RVector im(imag((dynamic_cast< CMatrix & > (*subSolutions_))[i + kIdx * nCurrentPattern]));
                    solutions_[i] += re * weights_[kIdx];
                    solutions_[i + nCurrentPattern] += im * weights_[kIdx];
                } else {
                    solutions_[i] += (dynamic_cast< RMatrix & > (*subSolutions_))[i + kIdx * nCurrentPattern] * weights_[kIdx];
                }
            }
        }
    }

    if (verbose_) std::cout << "Forward: ";
    swatch.stop(verbose_);
    // MEMINFO
}

template < class ValueType >
void DCMultiElectrodeModelling::calculateKAnalyt(const std::vector < ElectrodeShape * > & eA,
                                                 const std::vector < ElectrodeShape * > & eB,
                                                 Matrix < ValueType > & solutionK,
                                                 double k, int kIdx) const {

    uint nCurrentPattern = eA.size();
    if (solutionK.rows() < (kIdx + 1) * nCurrentPattern) {
        throwLengthError(WHERE_AM_I + " workspace size insufficient" + str(solutionK.rows())
            + " " + str((kIdx+1)*nCurrentPattern));
    }

    for (uint i = 0; i < nCurrentPattern; i ++) {
        solutionK[i + kIdx * nCurrentPattern] *= ValueType(0.0);
        if (eA[i]) solutionK[i + kIdx * nCurrentPattern] = exactDCSolution(*mesh_, eA[i], k, surfaceZ_, setSingValue_);
        if (eB[i]) solutionK[i + kIdx * nCurrentPattern] -= exactDCSolution(*mesh_, eB[i], k, surfaceZ_, setSingValue_);
    }
}

template < class ValueType >
void DCMultiElectrodeModelling::calculateK_(const std::vector < ElectrodeShape * > & eA,
                                            const std::vector < ElectrodeShape * > & eB,
                                            Matrix < ValueType > & solutionK, int kIdx){
    bool debug = false;
    Stopwatch swatch(true);

    if (debug) std::cout << "Starting calculateK ... " << std::endl;

    uint nCurrentPattern = eA.size();
    double k = kValues_[kIdx];

    if (solutionK.rows() < (kIdx+1) * nCurrentPattern) {
        throwLengthError(WHERE_AM_I + " workspace size insufficient" + str(solutionK.rows())
            + " " + str((kIdx+1) * nCurrentPattern));
    }

    if (analytical_) {
        return calculateKAnalyt(eA, eB, solutionK, k, kIdx);
    }

    if (debug) std::cout << "Building sparsity pattern ... " ;

    SparseMatrix < ValueType > S_;
    S_.buildSparsityPattern(*mesh_);

// MEMINFO

//** START  assemble matrix
    if (debug) std::cout << "Assembling system matrix ... " ;
    dcfemDomainAssembleStiffnessMatrix(S_, *mesh_, k);
    dcfemBoundaryAssembleStiffnessMatrix(S_, *mesh_, sourceCenterPos_, k);

    uint oldMatSize = mesh_->nodeCount();

    if (buildCompleteElectrodeModel_){
        int lastValidElectrode = electrodes_.size();
//         //** check for passive bodies that lay behind valid electrodes
//         for (uint i = 0; i < electrodes_.size(); i ++){
//             if (electrodes_[i]->id() < 0) {
//                 lastValidElectrode = i;
//                 if (verbose_) std::cout << "active electrode count: " << lastValidElectrode << std::endl;
//                 break;
//             }
//         }

        std::vector < ElectrodeShape * > elecs;
        for (Index i = 0; i < electrodes_.size(); i ++) elecs.push_back(electrodes_[i]);

        if (electrodeRef_ && electrodeRef_ != electrodes_[lastValidElectrode]) {
            electrodeRef_->setId(electrodes_.size());
//             std::cout << WHERE_AM_I << "//** FIXME this fails with passive bodies " << std::endl;
//             //** FIXME this fails with passive bodies
            elecs.push_back(electrodeRef_);
        }
        if (debug) std::cout << "Assembling complete electrode model ... " << std::endl;

        for (Index i = 0; i < passiveCEM_.size(); i ++) elecs.push_back(passiveCEM_[i]);

        if (vContactImpedance_.size() == 0){
            vContactImpedance_.resize(elecs.size(), 1.0); // Ohm
            // RVector vContactResistance(nElectrodes, 1.0); // Ohm
            // RVector vContactImpedance( nElectrodes, 1.0); // Ohm * m^2

            bool hasImp = checkIfMapFileExistAndLoadToVector("contactImpedance.map",
                                                             vContactImpedance_);
            if (hasImp){
                if (verbose_) std::cout << "Loaded: contactImpedance.map." << std::endl;
            }
        }

        assembleCompleteElectrodeModel(S_, elecs, oldMatSize, lastIsReferenz_,
                                           vContactImpedance_);

        potentialsCEM_.resize(nCurrentPattern, lastValidElectrode);

    } // end CEM

    this->assembleStiffnessMatrixDCFEMByPass(S_);

    assembleStiffnessMatrixHomogenDirichletBC(S_, calibrationSourceIdx_);

// MEMINFO
    //** END assemble matrix

    //** START solving
    SolverWrapper * solver;
    bool solverNeedsDelete = false;
    if (solver_ != nullptr){
        solver = solver_;
        solver->setMatrix(S_);
    } else {
        solver = new LinSolver(debug);
        solverNeedsDelete = true;
        dynamic_cast< LinSolver * >(solver)->setMatrix(S_, 1);
        if (debug) std::cout << "Factorize (" << solver->name() << ") matrix ... " << swatch.duration() << std::endl;
    }

// MEMINFO

    Vector < ValueType > sol(S_.cols());

    for (uint i = 0; i < nCurrentPattern; i ++){
        if (verbose_ && k == 0){
            std::cout << "\r " << i << " (" << swatch.duration(true) << "s)";
        }

        RVector rTmp(S_.rows(), 0.0);
        if (eA[i]) eA[i]->assembleRHS(rTmp,  1.0, oldMatSize);
        if (eB[i]) eB[i]->assembleRHS(rTmp, -1.0, oldMatSize);
        Vector < ValueType >rhs(rTmp);

        solver->solve(rhs, sol);

        /*if (i==4){
            __MS(*mesh_)
            __MS(min(rhs) << " " << max(rhs))
    	    __MS(min(sol) << " " << max(sol))

	    S_.save("S-gimlii_2.matrix");
             rhs.save("rhs_2.vec");
             save(sol, "sol_2.vec");
            __MS(eA[i])
             __MS(eB[i])
             mesh_->addData("sol" + str((kIdx) * nCurrentPattern + i), sol);
             mesh_->addData("rhs", rhs);
             mesh_->addData("soll", log(abs(sol)));
             mesh_->exportVTK("sol_2");
             exit(0);
         }
	 */
//         mesh_->addData("sol" + str((kIdx) * nCurrentPattern + i), sol);
//         S_.save("S-gimli.matrix");
//         rhs.save("rhs.vec");

//          save(rhs[i],  "rhs.vec");
//          save(sol, "sol.vec");

//         __MS(complex_)
//         __MS(min(sol)<< " " << max(sol))
//         __MS(i<< " " << k)
        if (norml2(S_ * sol - rhs) / norml2(rhs) > 1e-6){
//                 S_.save("S.matrix");
//                 save(rhs[i], "rhs.vec");
//                 save(sol, "sol.vec");
            std::cout   << " Ooops: Warning!!!! Solver: " << solver->name()
                            << " fails with rms(A *x -b)/rms(b) > tol: "
                            << norml2(S_ * sol - rhs)<< std::endl;

//             S_.save("S-gimli.matrix");
//             rhs.save("rhs.vec");
//             exit(0);
        }
        solutionK[i + kIdx * nCurrentPattern].setVal(sol, 0, oldMatSize);

        if (buildCompleteElectrodeModel_){
            //__MS(sol(oldMatSize, sol.size()- passiveCEM_.size()));
            potentialsCEM_[i] = TmpToRealHACK(sol(oldMatSize, sol.size() - passiveCEM_.size()));
        }

        // no need for setSingValue here .. numerical primpotentials have
        // some "proper" non singular value
//         if (setSingValue_){
//             if (eA[i]) eA[i]->setSingValue(solutionK[i], mesh_->cellAttributes(),  1.0, k);
//             if (eB[i]) eB[i]->setSingValue(solutionK[i], mesh_->cellAttributes(), -1.0, k);
//         }
    }
// MEMINFO
    // we dont need reserve the memory
    S_.clean();
    if (solverNeedsDelete == true){
        delete solver;
    }
}


void DCMultiElectrodeModelling::calculateK(const std::vector < ElectrodeShape * > & eA,
                                           const std::vector < ElectrodeShape * > & eB,
                                           RMatrix & solutionK, int kIdx){
    calculateK_(eA, eB, solutionK, kIdx);
}


void DCMultiElectrodeModelling::calculateK(const std::vector < ElectrodeShape * > & eA,
                                           const std::vector < ElectrodeShape * > & eB,
                                           CMatrix & solutionK, int kIdx){
    calculateK_(eA, eB, solutionK, kIdx);
}

void DCSRMultiElectrodeModelling::updateMeshDependency_(){
    DCMultiElectrodeModelling::updateMeshDependency_();
    if (primMeshOwner_ && primMesh_) {
        delete primMesh_;
        primMesh_ = 0;
    }
    if (primPot_) {
        if (verbose_) std::cout<< " updateMeshDependency:: cleaning primpot" << std::endl;

        primPot_->clear();
        if (primPotOwner_) {
            delete primPot_;
            primPot_ = NULL;
        }
    }
}

void DCSRMultiElectrodeModelling::updateDataDependency_(){
    DCMultiElectrodeModelling::updateDataDependency_();
    if (primPot_) {
        if (verbose_) std::cout<< " updateDataDependency:: cleaning primpot" << std::endl;
        primPot_->clear();
        if (primPotOwner_) {
            delete primPot_;
            primPot_ = NULL;
        }
    }
}

void DCSRMultiElectrodeModelling::setPrimaryMesh(const std::string & meshname){
    if (primPotFileBody_.find(NOT_DEFINED) == std::string::npos){
        primMesh_ = new Mesh;
        try {
            primMesh_->load(meshname);
        } catch (std::exception & e) {
            std::cerr << "Cannot load mesh: " << e.what() << std::endl;
            delete primMesh_;
        }
        primMeshOwner_ = true;
    }
}

void DCSRMultiElectrodeModelling::checkPrimpotentials_(const std::vector < ElectrodeShape * > & eA,
                                                       const std::vector < ElectrodeShape * > & eB){
    uint nCurrentPattern = eA.size();
    Stopwatch swatch(true);

    if (!primPot_) {

        //! First check if primPot can be recovered by loading binary matrix
        if (verbose_) std::cout << "Allocating memory for primary potential..." ;

        primPot_ = new RMatrix(nCurrentPattern * kValues_.size(), mesh_->nodeCount());
        primPot_->rowFlag().fill(0);
        primPotOwner_ = true;
        if (verbose_) std::cout << "... " << swatch.duration(true) << std::endl;

        if (primPotFileBody_.rfind(".bmat") != std::string::npos){
            std::cout << std::endl << "No primary potential for secondary field. Recovering " + primPotFileBody_ << std::endl;
            loadMatrixSingleBin(*primPot_, primPotFileBody_);
            std::cout << std::endl << " ... done " << std::endl;
// MEMINFO
            //! check for sizes here!!!
        } else {
            //! no binary matrix so calculate if topography present and no alternativ filename given
            if (topography() && primPotFileBody_.find(NOT_DEFINED) != std::string::npos){
                if (verbose_) {
                    std::cout << std::endl
                            << "No primary potential for secondary field calculation with topography."
                            << std::endl;
                }

                if (!primMesh_){
                    primMesh_ = new Mesh;
// dangerous can lead to extremly large meshes
                    *primMesh_ = mesh_->createP2();
//                    *primMesh_ = *mesh_;
                    primMeshOwner_ = true;

                    if (verbose_) {
                        std::cout<< "Creating P2-Primmesh:\t"; primMesh_->showInfos();
                    }
                }

                primMesh_->setCellAttributes(1.0);
                DCMultiElectrodeModelling f(this->dataContainer(), verbose_);
                f.setMesh(*primMesh_);
                RMatrix primPotentials;

                f.collectSubPotentials(primPotentials);
                f.calculate(*primDataMap_);

                if (verbose_){
                    std::cout << "Interpolating to secondary mesh" << std::endl;
                }

//                 for (Index i = 0; i < primPotentials.rows(); i ++ ){
//                     primMesh_->addData("p"+str(i), log(abs(primPotentials[i])));
//                 }

                interpolate(*primMesh_, primPotentials, mesh_->positions(), *primPot_, verbose_);
                primPot_->rowFlag().fill(1);

//                 for (Index i = 0; i < primPot_->rows(); i ++ ){
//                     mesh_->addData("s"+str(i), log(abs((*primPot_)[i])));
//                 }
//
//                 primMesh_->exportVTK("prim");
//                 mesh_->exportVTK("sec");
                delete primMesh_;
                primMesh_ = 0;
//                   save(*primPot_, "primPot");
            }
        }
    } else {// if (!primPot_)
        // we assume that given primPotentials are ok
        if (verbose_){
            std::cout << "Using existing primary potentials." << std::endl;
        }

        primPot_->rowFlag().fill(1);
    }

    //! First check ready!

    //! now check if all necessary primary potentials loaded or calculated (topography).
    //! On demand load from single file or calculate analytically
    if (primPot_->rows() != kValues_.size() * nCurrentPattern ||
         primPot_->cols() != mesh_->nodeCount()){

        std::cout << "Warning! primary potential matrix size is invalid for secondary field calculation." << std::endl
                  << "\tMatrix size is " << primPot_->rows() << " x " <<  primPot_->cols()
                  << "instead of " <<  kValues_.size() * nCurrentPattern << " x " << mesh_->nodeCount() << std::endl;

        primPot_->resize(kValues_.size() * nCurrentPattern, mesh_->nodeCount());
        primPot_->rowFlag().fill(0);
    }
    bool initVerbose = verbose_;
    //RVector prim(mesh_->nodeCount());

    for (uint kIdx = 0; kIdx < kValues_.size(); kIdx ++){
        double k = kValues_[kIdx];

        for (uint i = 0; i < nCurrentPattern; i ++){
            uint potID = (i + kIdx * nCurrentPattern);
            if (primPot_->rowFlag()[potID] == 0) {

            //std::cout << "not enough primPot entries " << primPot_->rows() << " " << nCurrentPattern << std::endl;
                //!** primary potential vector is unknown

                if (primPotFileBody_.find(NOT_DEFINED) != std::string::npos){
                //!** primary potential file body is NOT_DEFINED so we try to determine ourself
                    if (initVerbose){
                        std::cout << std::endl << "No primary potential for secondary field calculation. "
                                                  "Calculating analytically..." << std::endl;
                        initVerbose = false;
                    }
                    // PLS CHECK some redundancy here see DCMultiElectrodeModelling::calculateKAnalyt
                    if (eA[i]) (*primPot_)[potID]  = exactDCSolution(*mesh_, eA[i], k, surfaceZ_, setSingValue_);
                    if (eB[i]) (*primPot_)[potID] -= exactDCSolution(*mesh_, eB[i], k, surfaceZ_, setSingValue_);
                    //print("\r", k , "/", kValues_.size(), "/", nCurrentPattern);

                } else {
                    if (initVerbose){
                        std::cout << std::endl << "No primary potential for secondary field calculation. "
                                                  "Loading potentials." << std::endl;
                        initVerbose = false;
                    }
                //!** primary potential file body is given so we load it
                    if (k == 0.0){
                        //!** load 3D potential
                        if (initVerbose) std::cout << std::endl << "Loading primary potential: "
                                        << primPotFileBody_ + "." + str(i) + ".pot" << std::endl;
                        load((*primPot_)[potID], primPotFileBody_ + "." + str(i) + ".pot", Binary);
                    } else {
                        //!** else load 2D potential
                        //!** first try new style "name_Nr.s.pot"
                        if (!load((*primPot_)[potID], primPotFileBody_ + "." +
                                    str(kIdx * nCurrentPattern + i) + ".s.pot", Binary, false)){

                            if (!load((*primPot_)[potID], primPotFileBody_ + "." +
                                str(i) + "_" + str(kIdx) + ".pot", Binary)){
                                throwError(WHERE_AM_I + " neither new-style potential ("
                                    + primPotFileBody_ + ".XX.s.pot nor old-style ("
                                    + primPotFileBody_ + ".XX_k.pot) found");
                            }
                        }
                    } //! else load 2d pot
                } //! else load pot
                //** current primary potential is loaded or created, set flag to 1
                primPot_->rowFlag()[potID] = 1;
            } //! if primPot[potID] == 0
        } //! for each currentPattern
    } //** for each k

//     std::cout << swatch.duration() << std::endl;
//     exit(0);
}

void DCSRMultiElectrodeModelling::preCalculate(const std::vector < ElectrodeShape * > & eA,
                                                const std::vector < ElectrodeShape * > & eB){
    //! check for valid primary potentials, calculate analytical or numerical if nessecary
    checkPrimpotentials_(eA, eB);
    mesh1_ = *mesh_;
    mesh1_.setCellAttributes(1.0);
// MEMINFO
}

void DCSRMultiElectrodeModelling::calculateK(const std::vector < ElectrodeShape * > & eA,
                                             const std::vector < ElectrodeShape * > & eB,
                                             RMatrix & solutionK, int kIdx) {
    bool debug = false;
    if (complex_){
        THROW_TO_IMPL
    }
    Stopwatch swatch(true);
    double k = kValues_[kIdx];

    uint nCurrentPattern = eA.size();
    if (solutionK.rows() < (kIdx + 1) * nCurrentPattern) {
        throwLengthError(WHERE_AM_I + " workspace size insufficient" + str(solutionK.rows())
            + " " + str((kIdx + 1)*nCurrentPattern));
    }

    if (analytical_){
        calculateKAnalyt(eA, eB, solutionK, k, kIdx);
        return ;
    }
// MEMINFO

    RSparseMatrix S_;
    S_.buildSparsityPattern(*mesh_);
    bool singleVerbose = verbose_;

// MEMINFO

    dcfemDomainAssembleStiffnessMatrix(       S_, *mesh_, k);
    dcfemBoundaryAssembleStiffnessMatrix(     S_, *mesh_, sourceCenterPos_, k);
    assembleStiffnessMatrixHomogenDirichletBC(S_, calibrationSourceIdx_);
//     S_.save("S.mat");
//     exit(1);

    RSparseMatrix S1(S_);

//     RVector tmpRho(mesh_->cellAttributes());
//     mesh_->setCellAttributes(1.0);
    dcfemDomainAssembleStiffnessMatrix(  S1, mesh1_, k);
    dcfemBoundaryAssembleStiffnessMatrix(S1, mesh1_, sourceCenterPos_, k);
    assembleStiffnessMatrixHomogenDirichletBC(S1, calibrationSourceIdx_);

    // if (verbose_) std::cout << "Assembling system matrix (SR) ... " << std::endl;
//     if (verbose_) std::cout << "Assembling: " << swatch.duration() << std::endl;
    //mesh_->setCellAttributes(tmpRho);

// MEMINFO
    SolverWrapper *solver;
    bool solverNeedsDelete = false;

    if (solver_ != nullptr){
        solver = solver_;
        solver->setMatrix(S_);
    } else {
        solver = new LinSolver(debug);
        solverNeedsDelete = true;
        dynamic_cast< LinSolver * >(solver)->setMatrix(S_, 1);
        if (debug) std::cout << "Factorize (" << solver->name() << ") matrix ... " << swatch.duration() << std::endl;
    }


// MEMINFO

    //std::cout << WHERE_AM_I << std::endl;
    RVector rhs(S_.rows()), prim(rhs.size());

    for (uint i = 0; i < nCurrentPattern; i ++){

        if (primPot_->rows() <= (i + kIdx * nCurrentPattern)){
            throwError(WHERE_AM_I + " this should not happen, pls check primpots ");
            //** we do not have the primPot;
        } else {
            prim = (*primPot_)[i + kIdx * nCurrentPattern];
        }

        //** determine resistivity at the source location
        double rhoSourceA = 0.0, rhoSourceB = 0.0, rhoSource = 0.0;
        uint count = 0;
        if (eA[i]) {
            rhoSourceA = eA[i]->geomMeanCellAttributes();
            if (rhoSourceA > TOLERANCE){
                rhoSource += rhoSourceA; count++;
            } else {
                std::cout << eA[i]->id() << " "<< eA[i]->pos() << " "
                << eA[i]->geomMeanCellAttributes() << std::endl;
                std::cerr << WHERE_AM_I << " WARNING! rhoSourceA < TOLERANCE: " << std::endl;
            }
        }
        if (eB[i]) {
            rhoSourceB = eB[i]->geomMeanCellAttributes();
            if (rhoSourceB > TOLERANCE){
                rhoSource += rhoSourceB; count++;
            } else {
                std::cout << eB[i]->id() << " "<< eB[i]->pos() << " "
                    << eB[i]->geomMeanCellAttributes() << std::endl;
                std::cerr << WHERE_AM_I << " WARNING! rhoSourceB < TOLERANCE: " << std::endl;
            }
        }

//S_sig*u_s = (sig_0 * S1 - S_sig) u_p = sig_0 * S1 * u_p - S_sig * u_p
        rhoSource = rhoSource / count;
        prim *= rhoSource;
        rhs = S1 * prim / rhoSource - S_ * prim;

//         bool newWay = true;
//         if (newWay){
//         rhs = S1 * prim / rhoSource - S_ * prim;
//             //rhs = (1.0 / (rhoSource)) * S1 * prim - S_ * prim;
//         } else {

// //             RSparseMatrix Stmp(S_);
// //             RVector tmpRho2(mesh_->cellAttributes());
// //             for (uint t = 0; t < mesh_->cellCount(); t ++) {
// //                 if (std::fabs(mesh_->cell(t).attribute() - rhoSource) < 1e-10) {
// //                     //std::cout << "mesh_->cell(t).setAttribute(0.0); " << std::endl;
// //                     mesh_->cell(t).setAttribute(0.0);
// //                 } else {
// //                     mesh_->cell(t).setAttribute(1.0 /
// //                                 ( 1.0 / rhoSource - 1.0 / mesh_->cell(t).attribute())) ;
// //                 }
// //             }
// //             dcfemDomainAssembleStiffnessMatrix(Stmp, *mesh_, k, false);
// //             dcfemBoundaryAssembleStiffnessMatrix(Stmp, *mesh_, sourceCenterPos_, k);
// //             mesh_->setCellAttributes(tmpRho2);
// //
// //             rhs = Stmp * prim;
//         }

        //** fill calibration points
        for (uint j = 0; j < calibrationSourceIdx_.size(); j ++) {
            rhs[calibrationSourceIdx_[j]] = 0.0;
        }
        solutionK[i + (kIdx * nCurrentPattern)] *= 0.0;
        solver->solve(rhs, solutionK[i + (kIdx * nCurrentPattern)]);

        solutionK[i + (kIdx * nCurrentPattern)] += prim;

        if (singleVerbose) singleVerbose = false;
    }

    if (solverNeedsDelete == true){
        delete solver;
    }

}

} // namespace GIMLI{
