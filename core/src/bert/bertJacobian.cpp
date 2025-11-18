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

#include "bertJacobian.h"

#include <calculateMultiThread.h>
#include <elementmatrix.h>
#include <memwatch.h>
#include <meshentities.h>
#include <shape.h>
#include <stopwatch.h>

#include <shape.h>

namespace GIMLI{

#if USE_BOOST_THREAD
    #include <boost/thread.hpp>
    boost::mutex eraseMutex__;
#else
    #include <thread>
    #include <mutex>
    std::mutex eraseMutex__;
#endif

template < class ValueType > class CreateSensitivityColMT : public GIMLI::BaseCalcMT{
public:
  CreateSensitivityColMT(Matrix < ValueType >          & S,
                         const std::vector < Cell * >  & para,
                         const DataContainerERT        & data,
                         const Matrix < ValueType >    & pots,
                         const std::map< long, uint >  & currPatternIdx,
                         const RVector                 & weights,
                         const RVector                 & k,
                         bool calc1,
                         bool verbose)
    : BaseCalcMT(verbose), S_(&S), para_(&para), //cellMapIndex_ (&cellMapIndex),
        data_(&data), pots_(&pots), currPatternIdx_(&currPatternIdx),
        weights_(&weights), k_(&k), calc1_(calc1){
            nData_ = data.size();
            nElecs_ = data.sensorCount();
    }

    virtual ~CreateSensitivityColMT(){}

    virtual void calc(){

        if (calc1_){
            calc1();
        }else {
            calc2();
        }
    }

    virtual void calc2(){
        ElementMatrix < double > S_i;
        ElementMatrix < double > S1_i;

        Cell * cell = NULL;
        int modelIdx = 0;

        const Vector < ValueType > *va;
        const Vector < ValueType > *vb;
        const Vector < ValueType > *vm;
        const Vector < ValueType > *vn;
        Vector < ValueType > dummy((*pots_)[0].size(), ValueType(0));

        const RVector *da = &(*data_)("a");
        const RVector *db = &(*data_)("b");
        const RVector *dm = &(*data_)("m");
        const RVector *dn = &(*data_)("n");

        for (Index cellID = start_; cellID < end_; cellID ++) {

            cell    = (*para_)[cellID];
            modelIdx = cell->marker();

            if (modelIdx < 0) continue;

            S1_i.ux2uy2uz2(*cell);

            for (Index kIdx = 0; kIdx < weights_->size(); kIdx ++){
                S_i.u2(*cell);
                S_i *= (*k_)[kIdx] * (*k_)[kIdx];
                S_i += S1_i;
                int a = 0, b = 0, m = 0, n = 0;
                for (Index dataIdx = 0; dataIdx < nData_; dataIdx ++ ){

                    a = (int)(*da)[dataIdx];
                    b = (int)(*db)[dataIdx];
                    m = (int)(*dm)[dataIdx];
                    n = (int)(*dn)[dataIdx];

                    if (a > -1) va = &(*pots_)[a + nElecs_ * kIdx]; else va = &dummy;
                    if (b > -1) vb = &(*pots_)[b + nElecs_ * kIdx]; else vb = &dummy;
                    if (m > -1) vm = &(*pots_)[m + nElecs_ * kIdx]; else vm = &dummy;
                    if (n > -1) vn = &(*pots_)[n + nElecs_ * kIdx]; else vn = &dummy;

                    (*S_)[dataIdx][modelIdx] +=
                       S_i.mult((*va), (*vb), (*vm), (*vn)) * (*weights_)[kIdx];
                }
            }
        }
    }

    virtual void calc1(){
        // log(Debug, "Thread #" + str(_threadNumber) + ": on CPU " + str(schedGetCPU()) +
        //            " slice " + str(start_) + ":" + str(end_));
        bool haveCurrentPatterns = false;

        if (currPatternIdx_->size() * weights_->size() == pots_->rows()) {
        //** we have current and measurements pattern instead of pol pol potentials;
            haveCurrentPatterns = true;
        }

        ElementMatrix < double > S_i;
        Cell * cell = NULL;
        int modelIdx = 0;
        // Index si, sj;

        const Vector < ValueType > *va;
        const Vector < ValueType > *vb;
        const Vector < ValueType > *vm;
        const Vector < ValueType > *vn;

        const RVector *da = &(*data_)("a");
        const RVector *db = &(*data_)("b");
        const RVector *dm = &(*data_)("m");
        const RVector *dn = &(*data_)("n");

        Vector < ValueType > dummy((*pots_)[0].size(), ValueType(0));

        for (Index cellID = start_; cellID < end_; cellID ++) {

            cell    = (*para_)[cellID];
            modelIdx = cell->marker();

            if (modelIdx < 0) continue;

            if (verbose_) {
                // 	cout << "\r";
                // 	for (int i = 0; i < tNr; i ++) cout << "\t\t\t";
                // 	cout <<	cellID << "/" << para_->size() - 1;
            }

            S_i.ux2uy2uz2(*cell);
            Index cellNodeCount = cell->nodeCount();

//             ValueType tmpPotA = ValueType(0);
//             ValueType tmpPotM = ValueType(0);
            ValueType sum = ValueType(0);

            int a = 0, b = 0, m = 0, n = 0;

            double weightsFactor = 1.0;
            //** if weights_->size() > 1, assuming 2.5D so we need to double the weights
            //** we integrate from 0 to \infty but need we need from -\infty to \infty
            if (weights_->size() > 1) weightsFactor = 2.0;

            for (Index dataIdx = 0; dataIdx < nData_; dataIdx ++ ){

                if (haveCurrentPatterns){
                    a = currPatternIdx_->find(data_->electrodeToCurrentPattern(a, b))->second;
                    m = currPatternIdx_->find(data_->electrodeToCurrentPattern(m, n))->second;
                    b = -1;
                    n = -1;
                } else {
                    a = (int)(*da)[dataIdx];
                    b = (int)(*db)[dataIdx];
                    m = (int)(*dm)[dataIdx];
                    n = (int)(*dn)[dataIdx];
                }

                for (Index kIdx = 0; kIdx < weights_->size(); kIdx ++){
                    sum = ValueType(0);

                    if (a > -1) va = &(*pots_)[a + nElecs_ * kIdx]; else va = &dummy;
                    if (b > -1) vb = &(*pots_)[b + nElecs_ * kIdx]; else vb = &dummy;
                    if (m > -1) vm = &(*pots_)[m + nElecs_ * kIdx]; else vm = &dummy;
                    if (n > -1) vn = &(*pots_)[n + nElecs_ * kIdx]; else vn = &dummy;

                    (*S_)[dataIdx][modelIdx] += S_i.mult((*va), (*vb), (*vm), (*vn)) * (weightsFactor * (*weights_)[kIdx]);
                    continue;
//                     std::cout << cell->id() << std::endl;
                    for (Index i = 0; i < cellNodeCount; i ++){
//                         si = S_i.idx(i);
//                         std::cout << cell->node(i).id() << std::endl;
                        for (Index j = 0; j < cellNodeCount; j ++){
//                             sj = S_i.idx(j);
                            // very most time criticle section here
//                             tmpPotA  = (*va)[si] - (*vb)[si];
//                             tmpPotM  = (*vm)[sj] - (*vn)[sj];
//
//                             sum += S_i.getVal(i, j) * tmpPotA * tmpPotM;

//                             sum += S_i.getVal(i, j) *
//                                     ((*va)[si] - (*vb)[si]) *
//                                     ((*vm)[sj] - (*vn)[sj]);
                            sum += S_i.getVal(i, j) *
                                    ((*va)[S_i.idx(i)] - (*vb)[S_i.idx(i)]) *
                                    ((*vm)[S_i.idx(j)] - (*vn)[S_i.idx(j)]);

                      /*  std::cout << i<<" "<<j<<" "<< S_i.getVal(i, j)<<" "<< (*va)[S_i.idx(i)]<<" "<<
                        (*vb)[S_i.idx(i)]<<" "<< (*vm)[S_i.idx(j)]<<" "<< (*vn)[S_i.idx(j)] << std::endl;
                      */
                        }
                    }
                    if (isInfNaN(sum)){
                        std::cerr << WHERE_AM_I << std::endl;
                    }
// 	  if ((*cellMapIndex_)[cellID] - 2 < 0 ||
// 	       (*cellMapIndex_)[cellID] - 2 > nModel-1){
// 	    std::cerr << WHERE_AM_I << " index out of bounds: [0 -- " << nModel-1 << "]" << cellMapIndex[cellID] - 2 << std::endl;
// 	  }

//	  (*S_)[dataIdx][(*cellMapIndex_)[cellID] - 2] += sum * (data_->k(dataIdx) * 2.0 * (*weights_)[kIdx]);

                    /*! Dangerous!!! mt writing into matrix*/
	       //(*S_)[dataIdx][(*cellMapIndex_)[cellID]] += sum * (2.0 * (*weights_)[kIdx]);
                    {
//                   #ifdef HAVE_LIBBOOST_THREAD
//                   boost::mutex::scoped_lock lock(eraseMutex__); // slows down alot
//                   #endif

                        (*S_)[dataIdx][modelIdx] += sum * (weightsFactor * (*weights_)[kIdx]);
//                         std::cout << "b: " << dataIdx<<" "<<modelIdx<<" "<<(*S_)[dataIdx][modelIdx]<< " "
//                         << sum * (weightsFactor * (*weights_)[kIdx]) << std::endl;
                    }
//                     exit(1);
                } // for each k
//                      exit(1);
            } // for each data
        } // for each cellID
    }

protected:
    Matrix < ValueType >            * S_;
    const std::vector < Cell * >    * para_;
    const DataContainerERT          * data_;
    const Matrix < ValueType >      * pots_;
    const std::map< long, uint >    * currPatternIdx_;
    const RVector                   * weights_;
    const RVector                   * k_;
    uint                            nData_;
    uint                            nElecs_;
    bool                            calc1_;

};

bool lessCellMarker(const Cell * c1, const Cell * c2) { return c1->marker() < c2->marker(); }

template < class ValueType >
void createSensitivityCol_(Matrix < ValueType > & S,
                          const Mesh & mesh,
                          const DataContainerERT & data,
                          const Matrix < ValueType > & pots,
                          const RVector & weights,
                          const RVector & k,
                          std::vector < std::pair < Index, Index > > & matrixClusterIds,
                          uint nThreads, bool verbose){

MEMINFO

    Index nData  = data.size();
    Index nModel = max(mesh.cellMarkers()) + 1;
    Index maxRows = weights.size() * data.sensorCount();

    if (pots.rows() >= maxRows){
//         if (pots[0].size() < mesh.nodeCount()){
//             std::stringstream str; str << WHERE_AM_I << " potential matrix colsize to small. "
//                                        << pots[0].size()  << "< " << mesh.nodeCount() << std::endl;
//         }
    } else {
        std::stringstream str1; str1 << WHERE_AM_I << " potential matrix rowsize to small."
                                   << pots.rows() << " < " << maxRows << std::endl;
        throwLengthError(str1.str());
    }

    Stopwatch swatch(true);
    std::map< long, uint > currPatternIdx;
    //std::cout << "CreateSensitivityColMT " << nThreads << std::endl;

    std::vector< Cell * > cells(mesh.findCellByMarker(0, -1));
    std::sort(cells.begin(), cells.end(), lessCellMarker);

    double maxMemSize = max(0.0, getEnvironment("SENSMATMAXMEM", 0.0, verbose));
    double maxSizeNeeded = MByte(nData * nModel * sizeof(double));

    if (maxMemSize > 0 && verbose){
        std::cout << "Size of S: " << maxSizeNeeded << " MB" << std::endl;
    }

    //** avoid MT problems
    for (std::vector< Cell * >::iterator it = cells.begin();
         it != cells.end(); it ++){
        (*it)->pShape()->invJacobian();
    }

//     ShapeFunctionCache::instance().shapeFunctions(cells[0]->shape());
//     ShapeFunctionCache::instance().deriveShapeFunctions(cells[0]->shape(), 0);
//     ShapeFunctionCache::instance().deriveShapeFunctions(cells[0]->shape(), 1);
//     ShapeFunctionCache::instance().deriveShapeFunctions(cells[0]->shape(), 2);
//     ShapeFunctionCache::instance().shapeFunctions(*cells[0]);
//     ShapeFunctionCache::instance().deriveShapeFunctions(*cells[0], 0);
//     ShapeFunctionCache::instance().deriveShapeFunctions(*cells[0], 1);
//     ShapeFunctionCache::instance().deriveShapeFunctions(*cells[0], 2);
//     //cells[0]->createShapefunctionts();

    if (maxMemSize > 0 && maxMemSize < maxSizeNeeded){

        uint modelCluster = std::floor((double)nModel / (maxSizeNeeded / maxMemSize));

        if (modelCluster < 1) {
            throwError(WHERE_AM_I + " sorry, size of single sensitivity-row exceeds memory limitations.");
        }

        std::cout << "Size of S cluster: " << MByte(nData * modelCluster * sizeof(ValueType)) << " MB" << std::endl;
        std::cout << "Using model cluster " << nModel << " x " << modelCluster << std::endl;

        S.resize(nData, modelCluster);

        matrixClusterIds.clear();
        matrixClusterIds.push_back(std::pair < Index, Index >(nData, nModel));

        bool calc1 = getEnvironment("SENSMAT1", false, true);
        for (uint i = 0; i < nModel; i += modelCluster ){
MEMINFO
            Index start = i;
            Index end   = min(start + modelCluster, nModel);
            std::cout << " " << start << " " << end<< std::endl;

            S.resize(nData, end - start);

            std::vector< Cell * > cellsCluster(mesh.findCellByMarker(start, end));
            std::sort(cellsCluster.begin(), cellsCluster.end(), lessCellMarker);

MEMINFO
            // subtract marker start index
            for (std::vector< Cell * >::iterator it = cellsCluster.begin(); it != cellsCluster.end(); it ++){
                (*it)->setMarker((*it)->marker() - start);
            }

            S *= ValueType(0);
MEMINFO

            distributeCalc(CreateSensitivityColMT< ValueType >(S, cellsCluster,
                                                               data, pots,
                                                               currPatternIdx,
                                                               weights, k,
                                                               calc1,
                                                               verbose),
                           cellsCluster.size(), nThreads, verbose);

MEMINFO

            //** fight against the Lorenz butterfly
            //** 1e-8 is to coarse, need adaptive tolerance
            //S.round(1e-8);
MEMINFO

            S.save("sensPart_" + str(start) + "-" + str(end));

            matrixClusterIds.push_back(std::pair < Index, Index >(start, end));

            // add marker start index
            for (std::vector< Cell * >::iterator it = cellsCluster.begin(); it != cellsCluster.end(); it ++){
                (*it)->setMarker((*it)->marker() + start);
            }
MEMINFO
        }

        S.clear();
    } else {

        // __MS(nData << " "<< S.rows() << " "<<  nModel << " "<<  S.cols())
        if (S.rows() != nData || S.cols() != nModel) S.resize(nData, nModel);
        S *= ValueType(0);
MEMINFO
        if (verbose){
            std::cout << "S(" << numberOfCPU() << "/" << nThreads; //**check!!!
            #if USE_BOOST_THREAD
            std::cout << "-boost::mt";
            #else
            std::cout << "-std::mt";
            #endif
            std::cout << "): " << swatch.duration() << ":";
//swatch.stop(verbose);
        }
        bool calc1 = getEnvironment("SENSMAT1", false, true);
        if (useOMP()){

        } else {
            distributeCalc(CreateSensitivityColMT< ValueType >(S, cells, data,
                                                           pots, currPatternIdx,
                                                           weights, k, calc1, verbose),
                            cells.size(), nThreads, verbose);
        }
         if (verbose){
             swatch.stop(verbose);
         }
MEMINFO
        //** fight against the Lorenz butterfly
        //** 1e-8 is to coarse, need adaptive tolerance
        //S.round(1e-8);
    }
}

void createSensitivityCol(RMatrix & S,
                          const Mesh & mesh,
                          const DataContainerERT & data,
                          const RMatrix & pots,
                          const RVector & weights,
                          const RVector & k,
                          std::vector < std::pair < Index, Index > > & matrixClusterIds,
                          uint nThreads, bool verbose){
    createSensitivityCol_(S, mesh, data, pots, weights, k, matrixClusterIds, nThreads, verbose);
}

void createSensitivityCol(CMatrix & S,
                          const Mesh & mesh,
                          const DataContainerERT & data,
                          const CMatrix & pots,
                          const RVector & weights,
                          const RVector & k,
                          std::vector < std::pair < Index, Index > > & matrixClusterIds,
                          uint nThreads, bool verbose){
    createSensitivityCol_(S, mesh, data, pots, weights, k, matrixClusterIds, nThreads, verbose);
}


void sensitivityDCFEMSingle(const std::vector < Cell * > & para, const RVector & p1, const RVector & p2,
		       RVector & sens, bool verbose){
    uint nCells = para.size();
    if (sens.size() != nCells) sens.resize(nCells);

    ElementMatrix < double > S_i;
    double sum = 0.0, a_jk = 0.0, tmppot = 0.0;
    //  cout << nCells << std::endl;

    for (uint i = 0; i < nCells; i ++){
    //cout << "Nr. " << i << std::endl;
        S_i.ux2uy2uz2(*para[i]);

        sum = 0.0;
        for (int j = 0, jmax = para[i]->nodeCount(); j < jmax; j ++){
            for (int k = 0, kmax = para[i]->nodeCount(); k < kmax; k ++){
	       a_jk = S_i.getVal(j, k);
	       tmppot = p1[S_i.idx(j)] * p2[S_i.idx(k)];
	       sum += a_jk * tmppot;
	//	cout << "\tS_mn: " << a_jk << "\tp1*p2: " << tmppot << "\tp*S_mn: " << a_jk * tmppot << "\tsum: " << sum << std::endl;
            }
        }
        sens[i] = sum;
    }
}

RVector prepExportSensitivityData(const Mesh & mesh, const RVector & data, double logdrop){
    Index nModel = unique(sort(mesh.cellMarkers())).size();

    ASSERT_EQUAL(nModel, data.size())

    //data have always the right length since it comes from S directly
    RVector modelSizes(nModel, 0.0);
    for (Index i = 0; i < mesh.cellCount(); i ++ ){
        modelSizes[mesh.cell(i).marker()] += mesh.cell(i).size();
    }

    return logTransDropTol(data/modelSizes, logdrop, true)(mesh.cellMarkers());

    // //RVector tmp(data/mesh.cellSizes());
    // if ((uint)data.size() != (uint)mesh.cellCount()){

    //     throwLengthError(WHERE_AM_I + " Datasize missmatch: " + str(mesh.cellCount())+
    //                         " " + str(data.size()));
    // } else {
    //     //for (uint i = 0; i < tmp.size(); i ++) tmp[i] = tmp[i] / mesh.cell(i).shape().domainSize();
    // }
    // RVector tmp(data/mesh.cellSizes());

    // RVector s(sign(tmp));

    // double tmpMax = max(abs(tmp));
    // tmp /= tmpMax;

    // for (uint i = 0; i < tmp.size(); i ++) {
    //     tmp[i] = std::fabs(tmp[i] / logdrop);
    //     if (tmp[i] < 1.0) tmp[i] = 1.0;
    // }

    // tmp = log10(tmp);
    // tmp /= max(tmp) * s;
    // return tmp;
}

void exportSensitivityVTK(const std::string & fileName,
                          const Mesh & mesh, const RVector & data,
                          double logdrop){
    std::map< std::string, RVector > res;
    res.insert(std::make_pair("Sensitivity" ,
                              prepExportSensitivityData(mesh, data, logdrop)));
    mesh.exportVTK(fileName, res);
}

// void exportSensMatrixDC(const std::string & filename, const Mesh & mesh, const RMatrix & S) {
//     exportSensMatrixDC(filename, mesh, S
// }

void exportSensMatrixDC(const std::string & filename, const Mesh & mesh,
                        const RMatrix & S, const IVector & idx,
                        double logdrop) {
    std::map< std::string, RVector > res;

    for (std::map < std::string, RVector >::const_iterator
            it = mesh.dataMap().begin();
            it != mesh.dataMap().end(); it ++){
        res.insert(std::make_pair(it->first, it->second));
    }

    std::string add;
//     RVector tmp(mesh.cellCount());

    for (size_t i = 0; i < S.rows(); i ++) {
        if (i < 100000) add = "0";
        if (i < 10000) add = "00";
        if (i < 1000) add = "000";
        if (i < 100) add = "0000";
        if (i < 10) add = "00000";

//         for (uint j = 0; j < tmp.size(); j ++) tmp[j] = S[i][mesh.cell(j).marker()];

//#res.insert(std::make_pair("sens-" + add + str(i), log10(RVector(abs(tmp)))));

            res.insert(std::make_pair("sens-" + add + str(i),
                        prepExportSensitivityData(mesh, S[i], logdrop)));

   }
   mesh.exportVTK(filename, res);
}

RVector coverageDC(const RMatrix & sensMatrix) {
    RVector cov;
    if (sensMatrix.rows() > 0) {
        cov.resize(sensMatrix.cols(), 0.0);

        for (size_t i = 0; i < sensMatrix.rows(); i ++) {
            cov += abs(sensMatrix[i]);
        }
    } else {
        std::cout << "Sensmatrix invalid" << std::endl;
    }
    return cov;
}

RVector coverageDCtrans(const MatrixBase & S, const RVector & dd, const RVector & mm) {
    RVector cov;

    if (S.rows() > 0) {
        cov.resize(S.cols(), 0.0);
    } else {
        std::cout << "Sensmatrix invalid" << std::endl;
    }

    if (S.rtti() == GIMLI_MATRIX_RTTI){
        const RMatrix *Sl = dynamic_cast < const RMatrix * >(&S);

        for (size_t i = 0; i < S.rows(); i ++) {
            cov += abs((*Sl)[i] * dd[i]);
        }
    } else if (S.rtti() == GIMLI_SPARSE_MAP_MATRIX_RTTI){

        const RSparseMapMatrix * Sl = dynamic_cast< const RSparseMapMatrix * >(&S);

        for (RSparseMapMatrix::const_iterator it = Sl->begin(); it != Sl->end(); it ++){
            Index row = (*it).first.first;
            Index col = (*it).first.second;
            cov[col] += (*it).second * dd[row];
        }
    } else {
        CERR_TO_IMPL
    }

    return cov / abs(mm);
}

RVector createCoverage(const MatrixBase & S, const Mesh & mesh){
    return createCoverage(S, mesh, RVector(S.rows(), 1.0), RVector(S.cols(), 1.0));
}

RVector createCoverage(const MatrixBase & S, const Mesh & mesh,
                       const RVector & response, const RVector & model){

    RVector covModel(coverageDCtrans(S, 1.0 / response, 1.0 / model));
	RVector covMesh(covModel(mesh.cellMarkers()));
    if (mesh.cellCount() == model.size()) {
        covMesh /= mesh.cellSizes();
    } else {
        RVector modelCellSizes(covMesh.size(), 0.0);
        for (Index i = 0; i < mesh.cellCount(); i ++){
            Cell *c = &mesh.cell(i);
            modelCellSizes[c->marker()] += c->shape().domainSize();
        }
        if (min(modelCellSizes) > TOLERANCE){
            covMesh /= modelCellSizes;
        } else {
            log(Error, "Coverage fails:" + str(mesh.cellCount()) + " " + str(model.size()));
        }
    }

    return covMesh;
}
} // namespace GIMLI
