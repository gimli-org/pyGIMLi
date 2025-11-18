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

#include "modellingbase.h"

#include "datacontainer.h"
#include "mesh.h"
#include "regionManager.h"
#include "stopwatch.h"
#include "vector.h"
#include "vectortemplates.h"
#include "sparsematrix.h"

#include "calculateMultiThread.h"

namespace GIMLI{

ModellingBase::ModellingBase(bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
}

ModellingBase::ModellingBase(DataContainer & data, bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
    setData(data);
}

ModellingBase::ModellingBase(const Mesh & mesh, bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
    setMesh(mesh);
}

ModellingBase::ModellingBase(const Mesh & mesh, DataContainer & data, bool verbose)
    : dataContainer_(NULL), verbose_(verbose){
    init_();
    setData(data);
    setMesh(mesh);
}

ModellingBase::~ModellingBase() {
    // __MS("delete: " << this)
    if (ownRegionManager_) delete regionManager_;
    if (mesh_) delete mesh_;
    if (jacobian_ && ownJacobian_) delete jacobian_;
    if (constraints_ && ownConstraints_) delete constraints_;

}

void ModellingBase::init_() {
    regionManager_      = new RegionManager(verbose_);
    regionManagerInUse_ = false;

    // __MS("ModellingBase::init_: " << this <<  " rm:" << regionManager_)

    mesh_               = 0;
    jacobian_           = 0;
    constraints_        = 0;
    dataContainer_      = 0;

    nThreads_           = min(8, numberOfCPU()-2);
    nThreadsJacobian_   = 1;

    ownJacobian_        = false;
    ownConstraints_     = false;
    ownRegionManager_   = true;

    initJacobian();
    initConstraints();
}

void ModellingBase::setVerbose(bool verbose) {
    regionManager_->setVerbose(verbose);
    verbose_=verbose;
}

void ModellingBase::setThreadCount(Index nThreads) {
    nThreads_ = max(1, (int)nThreads);
    GIMLI::setThreadCount(nThreads);
}

Index ModellingBase::threadCount(){
    bool verbose = this->verbose();

    if (verbose){
        std::cout << "J(" << numberOfCPU() << "/" << this->nThreads_;
    #if USE_BOOST_THREAD
            std::cout << "-boost::mt";
    #else
            std::cout << "-std::mt";
    #endif
            std::cout << ") " << std::flush;
    }
    if (verbose){
        std::cout << std::endl;
    }
    return max(1, this->nThreads_);
}

void ModellingBase::setData(DataContainer & data){
    //if (dataContainer_) {
    dataContainer_ = &data;
//     } else {
//         dataContainer_ = new DataContainer(data);
//     }
    updateDataDependency_();
}

DataContainer & ModellingBase::data() const{
    if (dataContainer_ == 0){
        throwError(WHERE_AM_I + " no data defined");
    }
    return * dataContainer_;
}

RVector ModellingBase::startModel() {
    //*! Create startmodel from default builder (abstract may be overloaded)
    //*! Create startmodel from regionManger
    if (startModel_.size() == 0 && regionManager_){
        setStartModel(regionManager_->createStartModel());
    }

    if (startModel_.size() == 0){
        setStartModel(createDefaultStartModel());
    }

    if (startModel_.size() == 0){
        std::cout << "Warning! there is no startmodel defined." << std::endl;
    }
    return startModel_;
}

void ModellingBase::setStartModel(const RVector & startModel){
    startModel_=startModel;

    //** for mesh less inversion
    if (regionManager().parameterCount() != startModel_.size()){
       regionManager().setParameterCount(startModel_.size());
    }
}

void ModellingBase::createRefinedForwardMesh(bool refine, bool pRefine){
    this->initRegionManager();

    if (regionManager_->pMesh()){
        if (refine){
            if (pRefine){
                log(Info, "Create P2 refined mesh for forward tasks.");
                setMesh_(regionManager_->mesh().createP2());
            } else {
                log(Info, "Create H2 refined mesh for forward tasks.");
                setMesh_(regionManager_->mesh().createH2());
            }
        } else {
            log(Info, "copy mesh for forward tasks.");
            setMesh_(regionManager_->mesh());
        }
    } else {
        throwError("Cannot create a refined forward mesh since I have none.");
    }
}

void ModellingBase::setMesh(const Mesh & mesh, bool ignoreRegionManager) {
    Stopwatch swatch(true);
    if (regionManagerInUse_ && !ignoreRegionManager){
        // && holdRegionInfos e.g., just give it a try to ignore the regionmanager if necessary
        // if (ownRegionManager_ == true){
        //     __M
            regionManager_->setMesh(mesh);//#, ignoreRegionManger);
            if (verbose_) std::cout << "ModellingBase::setMesh() switch to regionmanager mesh" << std::endl;
            setMesh_(regionManager_->mesh());
        // } else {
        //     __MS("omiiting")
        // }
    } else {
        if (verbose_) std::cout << "ModellingBase::setMesh() copying new mesh ... ";
        setMesh_(mesh);
        if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
    }

    if (verbose_) std::cout << "FOP updating mesh dependencies ... ";
    startModel_.clear();

    if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
}

void ModellingBase::setMesh_(const Mesh & mesh, bool update){
    this->clearConstraints();

    if (!mesh_) mesh_ = new Mesh();

    if (update) deleteMeshDependency_();
    (*mesh_) = mesh;
    if (update) updateMeshDependency_();
}

void ModellingBase::deleteMesh(){
    if (mesh_) delete mesh_;
    mesh_ = 0;
}

void ModellingBase::initJacobian(){
    if (!jacobian_){
        jacobian_ = new RMatrix();
        ownJacobian_ = true;
    }
}

void ModellingBase::setJacobian(MatrixBase * J){
    if (!jacobian_ && ownJacobian_){
        delete jacobian_;
    }
    jacobian_ = J;
    ownJacobian_ = false;
}

void ModellingBase::setMultiThreadJacobian(Index nThreads){
    nThreadsJacobian_ = max(1, nThreads);
}

class JacobianBaseMT : public GIMLI::BaseCalcMT{
public:
    JacobianBaseMT(MatrixBase * J,
                   const ModellingBase & fop,
                   const RVector & resp,
                   const RVector & model,
                   bool verbose)
    : BaseCalcMT(verbose), J_(J), fop_(&fop), _resp(&resp),
      _model(&model) {
    }

    virtual ~JacobianBaseMT(){}

    virtual void calc(){
        log(Debug, "Thread #" + str(_threadNumber) + ": on CPU " + str(schedGetCPU()) +
                   " slice " + str(start_) + ":" + str(end_));

        RMatrix *J = dynamic_cast< RMatrix * >(J_);

        for (Index i = start_; i < end_; i ++) {
            RVector modelChange(*_model);
            modelChange[i] *= 1.05;
            RVector respChange(fop_->response_mt(modelChange, i));
            J->setCol(i, (respChange - *_resp) / (modelChange[i] - (*_model)[i]));
        }
    }

protected:
    MatrixBase              * J_;
    const ModellingBase     * fop_;
    const RVector           * _resp;
    const RVector           * _model;
};

void ModellingBase::createJacobian_mt(const RVector & model,
                                      const RVector & resp){
    THROW_TO_IMPL
    if (verbose_) std::cout << "Create Jacobian matrix (brute force, mt) ...";

    Stopwatch swatch(true);
    // double fak = 1.05;

    if (!jacobian_){
        this->initJacobian();
    }
    RMatrix *J = dynamic_cast< RMatrix * >(jacobian_);
    if (J->rows() != resp.size()){ J->resize(resp.size(), model.size()); }

    ALLOW_PYTHON_THREADS
    distributeCalc(JacobianBaseMT(jacobian_, *this, resp, model, verbose_),
                   jacobian_->rows(), nThreadsJacobian_, verbose_);
    swatch.stop();
    if (verbose_) std::cout << " ... " << swatch.duration() << " s." << std::endl;
}

void ModellingBase::createJacobian(const RVector & model,
                                   const RVector & resp){
    if (verbose_) std::cout << "Create Jacobian matrix (brute force) ...";

    Stopwatch swatch(true);
    double fak = 1.05;

    if (!jacobian_){
        this->initJacobian();
    }
    RMatrix *J = dynamic_cast< RMatrix * >(jacobian_);
    if (J->rows() != resp.size()){ J->resize(resp.size(), model.size()); }

    for (size_t i = 0; i < model.size(); i++) {
        RVector modelChange(model);
        modelChange[i] *= fak;

        RVector respChange(response(modelChange));

        if (::fabs(modelChange[i] - model[i]) > TOLERANCE){
            J->setCol(i, (respChange - resp) / (modelChange[i] - model[i]));
        } else {
            J->setCol(i, RVector(resp.size(), 0.0));
        }

//         __MS(i << " " << min(J->col(i)) << " " << max(J->col(i)))

//         for (size_t j = 0; j < resp.size(); j++){
//             if (::fabs(modelChange[i] - model[i]) > TOLERANCE){
//                 (*J)[j][i] = (respChange[j] - resp[j]) / (modelChange[i] - model[i]);
//             } else {
//                 (*J)[j][i] = 0.0;
//             }
//         }
    }

    swatch.stop();
    if (verbose_) std::cout << " ... " << swatch.duration() << " s." << std::endl;
}

void ModellingBase::createJacobian(const RVector & model){
    RVector resp;
    if (nThreadsJacobian_ > 1){
        resp = response_mt(model);
    } else {
        resp = response(model);
    }

    if (!jacobian_){
        this->initJacobian();
    } else {
        jacobian_->clear();
    }
    RMatrix *J = dynamic_cast< RMatrix * >(jacobian_);
    if (J->rows() != resp.size()){
        J->resize(resp.size(), model.size());
    }

    if (nThreadsJacobian_ > 1){
        return createJacobian_mt(model, resp);
    } else {
        return createJacobian(model, resp);
    }
}

void ModellingBase::initConstraints(){
    if (constraints_ == 0){
        constraints_ = new RSparseMapMatrix(0, 0, 0);
        ownConstraints_ = true;
    }
}

void ModellingBase::setConstraints(MatrixBase * C){
    if (!constraints_ && ownConstraints_){
        delete constraints_;
    }
    constraints_ = C;
    ownConstraints_ = false;
}

void ModellingBase::createConstraints(){
//     __MS(constraints_->rtti())
    this->regionManager().fillConstraints(constraintsRef());
}

void ModellingBase::clearConstraints(){
    if (constraints_) {
        constraints_->clear();
    }
}

MatrixBase * ModellingBase::constraints() {
    return constraints_;
}

MatrixBase * ModellingBase::constraints() const {
    return constraints_;
}

RSparseMapMatrix & ModellingBase::constraintsRef() const {
    if (!constraints_) throwError(WHERE_AM_I + " constraints matrix is not initialized.");
    return *dynamic_cast < RSparseMapMatrix *>(constraints_);
}

RSparseMapMatrix & ModellingBase::constraintsRef() {
    if (!constraints_) throwError(WHERE_AM_I + " constraints matrix is not initialized.");
    return *dynamic_cast < RSparseMapMatrix *>(constraints_);
}

RVector ModellingBase::createMappedModel(const RVector & model,
                                         double background) const {
    if (mesh_ == 0) throwError("ModellingBase has no mesh for ModellingBase::createMappedModel");

    // __MS("ModellingBase::createMappedModel: " << this << " rm: " << regionManager_
    //             << " " << model.size() << " " <<  mesh_->cellCount())

    if (model.size() == mesh_->cellCount()) {
        IVector cM(mesh_->cellMarkers());

        // test if model are cell values instead of model that needs mapping
        if (unique(sort(cM[cM > -1])).size() != model.size()) return model;

    }

    RVector cellAtts(mesh_->cellCount());

    int marker = -1;
    std::vector< Cell * > emptyList;
    mesh_->createNeighborInfos();

    for (Index i = 0, imax = mesh_->cellCount(); i < imax; i ++){
        marker = mesh_->cell(i).marker();
        if (marker >= 0) {
            if ((size_t)marker >= model.size()){
                mesh_->exportVTK("mapModelfail");
                std::cerr << WHERE_AM_I << std::endl
                          << "Wrong mesh here .. see mapModelfail.vtk" << std::endl
                          << *mesh_ << std::endl
                          << "mesh contains " << unique(sort(mesh_->cellMarkers())).size() << " unique markers. " << std::endl;
                throwLengthError(WHERE_AM_I + " marker >= than model.size() " + str(marker)
                       + " >= " + str(model.size()));
            }
            if (abs(model[marker]) < TOLERANCE){
                // Zero model can be valid
                // emptyList.push_back(&mesh_->cell(i));
            }
            cellAtts[i] = model[marker];

        } else {
            // general background without fixed values, fixed values will be set at the end
            if (marker == -1) {
                cellAtts[i] = 0.0;
                emptyList.push_back(&mesh_->cell(i));
            }
        }
    }

    // if background == 0.0 .. empty cells are allowed
    if (emptyList.size() == mesh_->cellCount() && background != 0.0){
        throwLengthError(WHERE_AM_I + " too many empty cells" + str(emptyList.size())
                       + " == " + str(mesh_->cellCount()));
    }

    if (background != 0.0){
        mesh_->prolongateEmptyCellsValues(cellAtts, background);
    }
    // __MS(regionManager_->regionCount())
    // exit(1);

    bool warned = false;
    // search for fixed regions
    for (Index i = 0, imax = mesh_->cellCount(); i < imax; i ++){
        // if (abs(cellAtts[i]) < TOLERANCE){ // this will never work since the prior prolongation
        if (mesh_->cell(i).marker() <= MARKER_FIXEDVALUE_REGION){
                // setting fixed values
            SIndex regionMarker = -(mesh_->cell(i).marker() -
                                    MARKER_FIXEDVALUE_REGION);

            // __MS(regionManagerInUse_ << " " << this <<  " " <<
            //         regionManager_->region(regionMarker)->fixValue())

            if (regionManagerInUse_){
                double val = regionManager_->region(regionMarker)->fixValue();
                if (warned == false){
                    log(Info, "fixing region: ", regionMarker," to: ", val);
                    warned = true;
                }
                cellAtts[i] = val;
                // __MS(val)
            } else {
                // temporay hack until fixed in modelling.py
                if (warned == false){
                    // __MS(cellAtts[i])
                    log(Warning, "** TMP HACk ** fixing region: ", regionMarker, " to: ", 1/0.058);
                    warned = true;
                    __MS(this << " rm:" << regionManager_)
                    THROW_TO_IMPL
                }
                cellAtts[i] = 1./0.058;
            }
        }
    }

//     mesh_->exportVTK("postmap", cellAtts);
    return cellAtts;
}

void ModellingBase::mapModel(const RVector & model, double background){
    // non readonly version"!!!!!!!!!!!
    mesh_->setCellAttributes(createMappedModel(model, background));
}

void ModellingBase::initRegionManager() {
    // clean this up .. this will fail for second try with mesh

    if (!regionManagerInUse_){
        if (mesh_){
            regionManager_->setMesh(*mesh_);
            this->setMesh_(regionManager_->mesh());
        }
        regionManagerInUse_ = true;
        // __MS("ModellingBase::initRegionManager " << this << " rm: " << regionManager_)
    }
}

void ModellingBase::setRegionManager(RegionManager * reg){
    if (reg){
        regionManagerInUse_ = true;
        delete regionManager_;
        regionManager_ = reg;
        ownRegionManager_ = false;
        // __MS("ModellingBase::setRegionManager reg " << this << " rm: " << regionManager_)
    } else {
        regionManagerInUse_ = false;
        regionManager_      = new RegionManager(verbose_);
        ownRegionManager_   = true;
        // __MS("ModellingBase::setRegionManager new " << this << " rm: " << regionManager_)
    }
    // __MS("MB:setRegionManager: " << this << " rm: " << regionManager_)
    // __MS(regionManagerInUse_ << " " << this)
    // THROW_TO_IMPL
}

const RegionManager & ModellingBase::regionManager() const {
    if (regionManager_ == 0) throwError("No RegionManager initialized");
    return *regionManager_;
}

RegionManager & ModellingBase::regionManager(){
    initRegionManager();
    return *regionManager_;
}


Region * ModellingBase::region(int marker){
    return regionManager().region(marker);
}

RVector ModellingBase::createStartVector() {
    DEPRECATED
    return this->createStartModel();
}

RVector ModellingBase::createStartModel() {
    return regionManager().createStartModel();
}

LinearModelling::LinearModelling(MatrixBase & A, bool verbose)
    : ModellingBase(verbose){//, A_(& A) {
        setJacobian(&A);
        this->regionManager().setParameterCount(A.cols());
}

RVector LinearModelling::response(const RVector & model) {
    if (jacobian_->cols() != model.size()){
        throwLengthError(WHERE_AM_I + " Jacobian col size != model.size()"
                                + str(jacobian_->cols()) + " != " + str(model.size()));
    }
    return *jacobian_ * model;
}

/* dummy function avoiding brute-force Jacobian as J=A is already there */
void LinearModelling::createJacobian(const RVector & model){}

RVector LinearModelling::createDefaultStartModel() {
    return RVector(jacobian_->cols(), 1.0);
}

} // namespace GIMLI
