/******************************************************************************
 *   Copyright (C) 2008-2024 by the GIMLi development team                    *
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

#include "regionManager.h"

#include "mesh.h"
#include "node.h"
#include "shape.h"
#include "sparsematrix.h"
#include "stopwatch.h"
#include "trans.h"
#include "vectortemplates.h"

namespace GIMLI{

Region::Region(SIndex marker, RegionManager * parent, bool single)
    : marker_(marker), parent_(parent),
    isBackground_(false), isSingle_(single), parameterCount_(0)
    , tM_(NULL) {
    init_();
    if (isSingle_) {
        parameterCount_ = 1;
        constraintType_ = 0;
        this->setModelControl(1.0);
    }
}

Region::Region(SIndex marker, const Mesh & mesh, RegionManager * parent)
    : marker_(marker), parent_(parent),
        isBackground_(false), isSingle_(false),
        parameterCount_(0), tM_(NULL) {
    init_();
    this->resize(mesh, marker);
}

Region::Region(SIndex marker, const Mesh & mesh, SIndex cellMarker, RegionManager * parent)
    : marker_(marker), parent_(parent),
        isBackground_(false), isSingle_(false),
        parameterCount_(0), tM_(NULL) {
    init_();
    this->resize(mesh, cellMarker);
}

Region::Region(const Region & region){
    copy_(region);
}

Region::~Region(){
   if (tM_ && ownsTrans_) delete tM_;
}

Region & Region::operator = (const Region & region){
    if (this != & region){
        copy_(region);
    } return *this;
}

void Region::init_() {

    isPermuted_     = false;
    constraintType_ = 1;
    startParameter_ = 0;
    endParameter_   = 0;
    parameterCount_ = 0;
    zWeight_        = 1.0;
    fixValue_       = 0.0;
    lowerBound_     = 0.0;
    upperBound_     = 0.0;
    mcDefault_      = 1.0;
    modelControl_   = 1.0;
    startDefault_   = 0.0;
    ownsTrans_ = true;
    _isInParaDomain = true;

    transString_    = "Log";
    tM_ = new TransLogLU < RVector >;
}

void Region::copy_(const Region & region){
    THROW_TO_IMPL
    marker_ = region.marker();
    isBackground_ = region.isBackground();
    isSingle_ = region.isSingle();
    isPermuted_  = region.isPermuted();
}

void Region::setBackground(bool background){
    // __MS(marker_ << " set "<< background << " is " << isBackground_)
    if (background != isBackground_) {
        isBackground_ = background;
        // __MS(marker_ << " is "<< isBackground_)

        parent_->recountParaMarker_();
        parent_->createParaDomain_();
        bounds_.clear();
        this->constraintWeights_.clear();
    }
}

void Region::setFixValue(double val){
    fixValue_ = val;
    isBackground_ = false;
    this->setBackground(true);
    this->constraintWeights_.clear();
}

void Region::setSingle(bool single){
    if (single != isSingle_) {
        markSingle(single);
        parent_->recountParaMarker_();
        parent_->createParaDomain_();
        bounds_.clear();
        this->constraintWeights_.clear();
    }
}

void Region::resize(const Mesh & mesh, SIndex cellMarker){
    if (cellMarker != marker_) _isInParaDomain = false;
    cells_ = mesh.findCellByMarker(cellMarker);
    bounds_.clear();

    if (!isBackground_ && !isSingle_){

        Index nBounds = mesh.boundaryCount();
        if (nBounds == 0) {
            std::cerr << "WARNING! no boundaries defined! run mesh.createNeighborInfos()" << std::endl;
        }
        bool leftParaId = false;
        bool rightParaId = false;

        for (Index i = 0; i < nBounds; i ++){
            leftParaId = false;
            rightParaId = false;

            if (mesh.boundary(i).leftCell() != NULL) {
                leftParaId = (mesh.boundary(i).leftCell()->marker() == cellMarker);
            }
            if (mesh.boundary(i).rightCell() != NULL){
                rightParaId = (mesh.boundary(i).rightCell()->marker() == cellMarker);
            }
            if (leftParaId && rightParaId) bounds_.push_back(&mesh.boundary(i));
        }
    }
    this->constraintWeights_.clear();
}

void Region::resize(const std::vector < Cell * > & cells){
    cells_ = cells;
    bounds_.clear();
    if (!isSingle_){
        // find new bounds_ without the mesh
        log(Error, WHERE_AM_I, "In use?");

    }
    this->constraintWeights_.clear();
}

void Region::countParameter(Index start){
//     __MS(marker()<< " " << fixValue())
    startParameter_ = start;
    if (isBackground_) {
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
//             __MS(cells_[i]->marker() << " " << fixValue_)
            if (abs(fixValue_) > TOLERANCE) {
                if (cells_[i]->marker() >= -1) {
                    // set only for positive marker that not already fixed regions
                    cells_[i]->setMarker(MARKER_FIXEDVALUE_REGION - marker());
                }
            } else {
                cells_[i]->setMarker(-1);
            }
        }
        bounds_.clear();
        parameterCount_ = 0;
    } else if (isSingle_) {
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
            cells_[i]->setMarker(start);
        }
        bounds_.clear();
        parameterCount_ = 1;
    } else {
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
            cells_[i]->setMarker(start + i);
        }
        parameterCount_ = cells_.size();
    }

    // reset attributes to allow for clean mapping
    for (Index i = 0, imax = cells_.size(); i < imax; i ++) cells_[i]->setAttribute(0.0);

    endParameter_ = start + parameterCount_;
    // __MS(this->marker_ << " "<< this->isSingle_)
    // modelControl_.resize(parameterCount_, mcDefault_);
    startModel_.resize(parameterCount_, startDefault_);
    paraIDs_ = IndexArray(parameterCount_);
    for (Index i = 0; i < paraIDs_.size(); i ++) paraIDs_[i] = start + i;
    // std::cout << WHERE_AM_I << " " << marker_ << " " << parameterCount_ << " " << startParameter_ << " " << endParameter_ <<  std::endl;
}

void Region::permuteParameterMarker(const IndexArray & p){
    for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
        SIndex m = cells_[i]->marker();
        if (m >= 0) {
            ASSERT_RANGE((Index)m, 0, p.size())
//             __MS(cells_[i]->id() << ": " << m << " - "<< p[m] )
            cells_[i]->setMarker(p[m]);
        }
    }
    isPermuted_ = true;
    for (Index i = 0; i < paraIDs_.size(); i ++) paraIDs_[i] = p[paraIDs_[i]];
}

//################ Start values
void Region::setStartModel(const RVector & start){
    if (isBackground_){
        // log(Warning, "Region Nr:", marker_, " is background and should not get a starting model.");
        return;
    }
    if (start.size() == parameterCount_){
       startModel_ = start;
    } else {
        throwLengthError(WHERE_AM_I + " sizes missmatch for region " +
        str(this->marker_) + " " +str(start.size()) + " != " + str(parameterCount_));
    }
}

void Region::setStartModel(double start){
    startDefault_ = start;
    this->setStartModel(RVector(parameterCount_, start));
}

void Region::fillStartModel(RVector & vec){
    if (isBackground_) return;
    if (startModel_.size() != parameterCount_){
        std::cerr << "WARNING! starting value for region " << marker_ << " not set. "<< std::endl;
    } else {
        if (isSingle_){
            vec[startParameter_] = startModel_[0];
        } else {
            for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
                vec[cells_[i]->marker()] = startModel_[i];
            }
        }
    }
}

//################ Model behaviour
void Region::setModelControl(double val){
    if (isBackground_){
        // log(Warning, "Region Nr:", marker_, " is background and should not get a model control.");
        val = 1.0;
    }

    if (val < TOLERANCE) val = 1.0;
    mcDefault_ = val;
    modelControl_ = val;
}

void Region::setModelControl(const RVector & mc){
    log(Error, "don't use it. modelControl is scalar");
    if (isBackground_){
        // log(Warning, "Region Nr:", marker_, " is background and should not get a model control.");
        return;
    }
    modelControl_ = 1.0;
}

void Region::fillModelControl(RVector & vec){
    if (isBackground_) return;
    if (isSingle_){
        vec[startParameter_] = 1.0;
    } else{
        for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
            vec[cells_[i]->marker()] = 1.0;
        }
    }
}

//################ constraints behaviour
Index Region::constraintCount() const {
    if (isBackground_ ) return 0;
    if (isSingle_ && constraintType_ == 0) return 0;
    if (isSingle_ && constraintType_ == 1) return 1;

    if (constraintType_ == 0 || constraintType_ == 2 || constraintType_ == 20) return parameterCount();
    if (constraintType_ == 10) return bounds_.size() + parameterCount();
    return bounds_.size();
}

void Region::fillConstraints(RSparseMapMatrix & C, Index startConstraintsID){
    // __MS(isBackground_ << " " << isSingle_ << " " << constraintType_ <<  " " <<
    //      startConstraintsID << " " << startParameter_)
    if (isBackground_ ) return;

    if (isSingle_ && constraintType_ == 0) return;

    if (isSingle_ && constraintType_ == 1) {
        // __MS(startConstraintsID << " " << startParameter_)
        C[startConstraintsID][startParameter_] = 1.0;
        return;
    }

    double cMixRatio = 1.0; // for mixing 1st or 2nd order with 0th order (constraintTypes 10 and 20)
    if (constraintType_ == 10 || constraintType_ == 20) cMixRatio = 1.0; //**retrieve from properties!!!

    if (constraintType_ == 0 || constraintType_ == 20){ //purely 0th or mixed 2nd+0th
        if (isPermuted_){
            std::set < SIndex > para;
            for (Index i = 0, imax = cells_.size(); i < imax; i ++) {
                para.insert(cells_[i]->marker());
            }
            Index i = 0;
            for (std::set< SIndex >::iterator it = para.begin(); it!= para.end(); it++, i ++){
                C[startConstraintsID + i][(*it)] = cMixRatio;
            }
        } else {
            for (Index i = 0; i < parameterCount_; i++) {
                // this fails after parameter permution
                C[startConstraintsID + i][startParameter_ + i] = cMixRatio;
            }
        }
        if (constraintType_ == 0) return;
    }

    SIndex leftParaId = -1, rightParaId = -1;
    if (constraintType_ == 2 || constraintType_ == 20) { //** 2nd order constraints (opt. mixed with 0th)
        for (auto & it : this->bounds_){
            Boundary * b = (it);
            leftParaId = -1;
            rightParaId = -1;
            if (b->leftCell() ) leftParaId  = b->leftCell()->marker();
            if (b->rightCell()) rightParaId = b->rightCell()->marker();

            if (isPermuted_){
            // better check if cell is part of this region or not ..  mesh.data(regionMarker)
                if (leftParaId != rightParaId){
                    C[leftParaId][rightParaId] = -1;
                    C[rightParaId][leftParaId] = -1;
                    C[leftParaId][leftParaId] += 1;
                    C[rightParaId][rightParaId] += 1;
                }
            } else {
                // unsure if necessary ..bounds_ should be valid
                if (leftParaId  >= (int)startParameter_ && leftParaId  < (int)endParameter_ &&
                    rightParaId >= (int)startParameter_ && rightParaId < (int)endParameter_ &&
                        leftParaId != rightParaId){
                    C[leftParaId][rightParaId] = -1;
                    C[rightParaId][leftParaId] = -1;
                    C[leftParaId][leftParaId] += 1;
                    C[rightParaId][rightParaId] += 1;
                }
            }
        }
        return;
    }
    //** 1st order constraints (opt. combined with 0th order)
    Index cID = startConstraintsID;
    for (auto & it : this->bounds_){
        Boundary * b = (it);

        leftParaId = -1;
        rightParaId = -1;
        if (b->leftCell()) leftParaId  = b->leftCell()->marker();
        if (b->rightCell()) rightParaId = b->rightCell()->marker();

        if (isPermuted_){
            // unsure if necessary ..bounds_ should be valid
            // better check if cell is part of this region or not ..  mesh.data(regionMarker)
            if (leftParaId != rightParaId){
                C[cID][leftParaId] = 1;
                C[cID][rightParaId] = -1;
            }
        } else{
            //if (leftParaId > -1 && rightParaId > -1 && leftParaId != rightParaId) {
            if (leftParaId  >= (int)startParameter_ && leftParaId  < (int)endParameter_ &&
                rightParaId >= (int)startParameter_ && rightParaId < (int)endParameter_ &&
                    leftParaId != rightParaId){
                C[cID][leftParaId] = 1;
                C[cID][rightParaId] = -1;
            }
        }
        cID ++;
    }
    if (constraintType_ == 10) { //** combination with 0th order
        for (Index i = 0; i < parameterCount_; i++) {
            // this fails after parameter permution ..
            //CR: dont know what C==10 is supposed to be to fix it
            if (isPermuted_){
                __MS("Ctype 10 are untested for permuted region.")
            }
            C[cID][startParameter_ + i] = cMixRatio;
            cID++;
        }

    }
}

void Region::setConstraintType(Index type) {
    constraintType_ = type;
}

void Region::setConstraintWeights(double val){
    return setConstraintWeights(RVector(this->constraintCount(), val));
}

void Region::setConstraintWeights(const RVector & cw){
    //std::cout << "Region::setConstraintsWeight(const RVector & sw) " << sw.size() << " " <<  this->constraintCount() << std::endl;
    if (isBackground_){
        // log(Warning, "Region Nr:", marker_, " is background and should not get a constraints weight.");
        return;
    }
    if (cw.size() == this->constraintCount()){
        zWeight_ = 1.0;
        this->constraintWeights_ = cw;
    } else {
        throwLengthError(WHERE_AM_I + " " + str(cw.size()) + " != " + str(constraintCount()));
    }
}

const RVector & Region::constraintWeights(){
    if (this->constraintWeights_.size() != this->constraintCount()){
        this->_createConstraintWeights();
    }
    return this->constraintWeights_;
}

void Region::fillConstraintWeights(RVector & vec, Index cIDStart){
    if (isBackground_) return;
    vec.setVal(this->constraintWeights(), cIDStart);
}

void Region::_createConstraintWeights(){
    if (isBackground_ || (constraintType() == 0) ||
        (constraintType() == 2)) return;

    this->constraintWeights_.resize(constraintCount(), 1.0);

    int dim = parent_->mesh().dim();

    Index cId = 0;
    for (auto & it : this->bounds_){
        Boundary * b = (it);

        double cWeight = this->modelControl_;
        if (b->marker() != 0){
            if (this->parent_->interfaceConstraints().count(b->marker())){
                cWeight = this->parent_->interfaceConstraints().at(b->marker());
            } else {
               cWeight = 0.0;
            }
        }


        double zDir = std::fabs(b->norm()[dim-1]); //! z-component

        this->constraintWeights_[cId] = cWeight*(1.0 + zDir * (zWeight_ - 1.0)); //! rather linear for bigger angles
        cId ++;
    }
}

void Region::fillBoundaryNorm(std::vector< RVector3 > & vnorm, Index boundStart){
    log(Warning, WHERE_AM_I, "Who use this. Is needed?"); //190505
    if (isBackground_ || isSingle_ || (constraintType() == 0)) return;

    for (Index i = 0, imax = bounds_.size(); i < imax; i ++) {
        vnorm[i + boundStart] = bounds_[i]->norm();
    }
}

void Region::fillBoundarySize(RVector & vec, Index boundStart){
    log(Warning, WHERE_AM_I, "Who use this. Is needed?"); //190505
    if (isBackground_ || isSingle_ || (constraintType() == 0)) return;

    for (Index i = 0, imax = bounds_.size(); i < imax; i ++) {
        vec[i + boundStart] = bounds_[i]->shape().domainSize();
    }
}

void Region::setTransModel(Trans< RVector > & tM){

    if (isBackground_){
        // log(Warning, "Region Nr:", marker_, " is background and should not get a model transformation.");
        return;
    }

    if (tM_ && ownsTrans_) delete tM_;
    tM_ = & tM;
    parent_->setLocalTransFlag(true);
    ownsTrans_ = false;
}

void Region::setModelTransStr_(const std::string & val){
    if (isBackground_){
        // log(Warning, "Region Nr:", marker_, " is background and should not get a model transformation.");
        return;
    }
    transString_ = val;
    delete tM_; tM_ = NULL;

    if (val == "lin" || val == "Lin"){
        tM_ = new Trans< RVector >();
    } else if (val == "log" || val == "Log"){
        tM_ = new TransLogLU< RVector >(lowerBound_, upperBound_);
    } else if (val == "cot" || val == "Cot" || val == "tan" || val == "Tan"){
        tM_ = new TransCotLU< RVector >(lowerBound_, upperBound_);
    } else  {
        throwLengthError(WHERE_AM_I + " : " + val + ". Available are: lin, log, cot/tan.");
    }
    parent_->setLocalTransFlag(true);

    ownsTrans_ = true;
}

void Region::setLowerBound(double lb){
    lowerBound_ = lb;
    setModelTransStr_(transString_);
}

void Region::setUpperBound(double ub){
    upperBound_ = ub;
    setModelTransStr_(transString_);
}

void Region::setParameters(double start, double lb, double ub, std::string transString){
    if (lb < ub) {
        if ((start <= lb) | (start >= ub)) {
            std::cout << "WARNING! starting model not within bounds! readjusting" << std::endl;
            setStartModel(std::sqrt(lb * ub));
        } else {
            setStartModel(start);
        }
        lowerBound_ = lb;
        upperBound_ = ub;
        if (transString.size() > 0) { // any given
            setModelTransStr_(transString);
        } else { // use preset otherwise
            setModelTransStr_(transString_);
        }
    } else {
        throwError(WHERE_AM_I + " bounds not matching: " + str(lb) + ">=" + str(ub));
    }
}

//******************************************************************************
RegionManager::RegionManager(bool verbose) : verbose_(verbose), mesh_(NULL){
    paraDomain_ = new Mesh();
    parameterCount_ = 0;
    haveLocalTrans_ = false;
    isPermuted_ = false;
    localTransHaveChanges_ = true;
    interRegionConstraintZWeights_ = 1.0;
}

RegionManager & RegionManager::operator = (const RegionManager & rm){
    if (this != &rm){
        copy_(rm);
    } return *this;
}

void RegionManager::copy_(const RegionManager & rm){
    CERR_TO_IMPL
}

RegionManager::~RegionManager(){
    clear();
    if (paraDomain_) delete paraDomain_;
}

const Mesh & RegionManager::mesh() const {
    if (mesh_== 0){
        throwError("RegionManager knows no mesh.");
    }
    return *mesh_;
}

Region * RegionManager::region(SIndex marker){
    if (regionMap_.count(marker) == 0){
        throwError(WHERE_AM_I + " no region with marker " + str(marker));
    } else {
        return regionMap_[marker];
    }
    return NULL;
}

void RegionManager::clear(){
    for (std::map< SIndex, Region* >::const_iterator it = regionMap_.begin(), end = regionMap_.end();
          it != end; it ++){
        delete it->second;
    }
    regionMap_.clear();
    interRegionInterfaceMap_.clear();
    interRegionConstraints_.clear();
    interfaceConstraints_.clear();
    isPermuted_ = false;
    _cWeights.clear();

    if (paraDomain_) { paraDomain_->clear(); }
    if (mesh_) { delete mesh_; mesh_ = NULL; }
}


void RegionManager::setMesh(const Mesh & mesh, bool holdRegionInfos){

    if (!holdRegionInfos){
        if (verbose_) std::cout << "Reset region parameter" << std::endl;
        this->clear();
    }

    Stopwatch swatch(true);
    if (verbose_) std::cout << "RegionManager copying mesh ...";

    if (mesh_) delete mesh_;
    mesh_ = new Mesh(mesh);

    if (verbose_){
        std::cout << swatch.duration(true) << " s " << std::endl;
        std::cout << "create NeighborInfos ... ";
    }
    mesh_->createNeighborInfos();
    if (verbose_) {
        std::cout << swatch.duration(true) << " s " << std::endl;
        std::cout << "analyzing mesh ... ";
    }

    //** looking for and create regions
    IVector regions(unique(sort(mesh_->cellMarkers())));

    if (verbose_) std::cout << regions.size() << " regions." << std::endl;

    bool singleOnly = false;
    if (regions.size() > 50){
        singleOnly = true;
        log(Info, "More than 50 regions, so we assume single regions only.");
    }

    std::map < int, std::vector< Cell * > > markerCellVectorMap;
    if (singleOnly){
        for(std::vector< Cell * >::const_iterator it = mesh_->cells().begin();
                                               it != mesh_->cells().end(); it++){
            if (! markerCellVectorMap.count((*it)->marker())){
                markerCellVectorMap[(*it)->marker()] = std::vector< Cell * >();
            }
            markerCellVectorMap[(*it)->marker()].push_back((*it));
        }
    }

    for (Index i = 0; i < regions.size(); i ++){
        if (singleOnly){
            createSingleRegion_(regions[i], markerCellVectorMap[regions[i]]);
        } else {
            createRegion_(regions[i], *mesh_, regions[i]);
        }
    }
    //** looking for and create region interfaces

    //** looking for and create inter-region interfaces
    this->findInterRegionInterfaces();

    if (singleOnly){
        log(Info, "Applying *:* interregion constraints.");
        for (Index i = 0; i < regions.size(); i ++){
            for (Index j = 0; j < regions.size(); j ++){
                if (i != j){
                    setInterRegionConstraint(regions[i], regions[j], 1.0);
                                // if (verbose_) std::cout << minRegion[i] << " <-> "
                                //                         << maxRegion[j] << " weight:"
                                //                         << toDouble(row[2]) << std::endl;
                }
            }
        }
    }
    this->recountParaMarker_();
    this->createParaDomain_();
}

Region * RegionManager::createSingleRegion_(SIndex marker, const std::vector < Cell * > & cells){
    Stopwatch swatch(true);
    Region * region = NULL;
    if (regionMap_.count(marker) == 0){
        region = new Region(marker, this, true);
        regionMap_.insert(std::make_pair(marker, region));
    } else {
        THROW_TO_IMPL
        region = regionMap_[marker];
    }
    if (cells.size() > 0){
        region->resize(cells);
    }
    return region;
}

Region * RegionManager::createRegion_(SIndex marker, const Mesh & mesh, SIndex cellMarker){
    Region * region = NULL;

    if (regionMap_.count(marker) == 0){
        __MS(marker)
        region = new Region(marker, mesh, cellMarker, this);
        regionMap_.insert(std::make_pair(marker, region));
    } else {
        region = regionMap_[marker];
        region->resize(mesh, cellMarker);
        //std::cerr << WHERE_AM_I << " Region with marker " << marker << " already exists." << std::endl;
    }
    return region;
}

Region * RegionManager::addRegion(SIndex marker){
    Region * region = createSingleRegion_(marker, std::vector < Cell * > ());
    recountParaMarker_(); //** make sure the counter is right
    return region;
}

Region * RegionManager::addRegion(SIndex marker, const Mesh & mesh, SIndex cellMarker){
    Region * region = createRegion_(marker, mesh, cellMarker);
    recountParaMarker_(); //** make sure the counter is right
    return region;
}

void RegionManager::createParaDomain_(){
    if (verbose_) std::cout << "creating para domain ... ";
    Stopwatch swatch(true);
    IndexArray cellIdx;
    cellIdx.reserve(mesh_->cellCount());

    for (auto & x: this->regionMap_){
        if (!x.second->isBackground() && x.second->isInParaDomain()){
            for (auto & c: x.second->cells()){
                cellIdx.push_back(c->id());
            }
        }
    }

    paraDomain_->createMeshByCellIdx(*mesh_, cellIdx);
    if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
}

void RegionManager::permuteParameterMarker(const IVector & p){
    isPermuted_ = true;
    for (auto & x: this->regionMap_){
        x.second->permuteParameterMarker(p);
    }
    this->createParaDomain_();
}

void RegionManager::recountParaMarker_(){
    Index count = 0;
    for (auto & x: this->regionMap_){
        x.second->countParameter(count);
        count += x.second->parameterCount();
    }
}

void RegionManager::findInterRegionInterfaces(){
    // find all boundaries per inter region boundary
    // and store them < <left, right>, [ptr Boundary]>

    interRegionInterfaceMap_.clear();
    interRegionConstraints_.clear();

    //** having fun with stl
    std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > > ::iterator iRMapIter;

    for (auto & bIter: mesh_->boundaries()){
        Boundary & b = *bIter;

        if (b.leftCell() && b.rightCell()){
            if (b.leftCell()->marker() != b.rightCell()->marker()){
                SIndex minMarker = min(b.leftCell()->marker(),
                                       b.rightCell()->marker());
                SIndex maxMarker = max(b.leftCell()->marker(),
                                       b.rightCell()->marker());

                iRMapIter = interRegionInterfaceMap_.find(std::pair< SIndex, SIndex >(minMarker, maxMarker));

                if (iRMapIter == interRegionInterfaceMap_.end()){
                    interRegionInterfaceMap_.insert(std::pair<
                        std::pair< SIndex, SIndex >,                       std::list < Boundary * > > (std::pair< SIndex, SIndex>(minMarker, maxMarker),
                                        std::list< Boundary* >()));
                }

                interRegionInterfaceMap_[std::pair< SIndex, SIndex > (minMarker, maxMarker)].push_back(&b);
            }
        }
    }
    // if (verbose_ && interRegionInterfaceMap_.size()){
    //     for (std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > >::iterator
    //             it  = interRegionInterfaceMap_.begin();
    //             it != interRegionInterfaceMap_.end(); it ++){
    //         std::cout << "(" << it->first.first << "," << it->first.second << ") "
    //                          << it->second.size() << " boundaries. " << std::endl;
    //     }
    // }
}

void RegionManager::fillStartModel(RVector & vec){
    if (vec.size() != parameterCount()) vec.resize(parameterCount());

    for (auto & x: this->regionMap_){
        x.second->fillStartModel(vec);
    }
}

RVector RegionManager::createStartModel(){
    RVector vec(parameterCount());
    fillStartModel(vec);
    return vec;
}

void RegionManager::fillModelControl(RVector & vec){
    //!** no regions: fill 0th-order constraints
    if (regionMap_.empty()){
        vec.resize(this->parameterCount(), 1.0);
        return;
    }

    if (vec.size() != parameterCount()) vec.resize(parameterCount(), 1.0);

    for (auto & x: this->regionMap_){
        x.second->fillModelControl(vec);
    }
}

RVector RegionManager::createModelControl(){
    RVector vec(parameterCount(), 1.0);
    fillModelControl(vec);
    return vec;
}

RVector RegionManager::constraintWeights(){
    if (this->_cWeights.size() == 0){
        log(Error, "no cWeights defined. You should create constraints matrix first.");
    }
    return this->_cWeights;
}

void RegionManager::fillConstraintWeights(RVector & vec){
    log(Error, WHERE_AM_I, "in use??");
    if (this->_cWeights.size() == 0){
        log(Error, "no cWeights defined. You should create constraints matrix first.");
    }
    vec = this->_cWeights;
    return;
}

void RegionManager::fillBoundarySize(RVector & vec){
    log(Error, WHERE_AM_I, "in use??");
    vec.resize(constraintCount(), 1.0);

    Index boundCount = 0;
    for (auto & x: this->regionMap_){
        x.second->fillBoundarySize(vec, boundCount);
        boundCount += x.second->constraintCount();
    }
}

Index RegionManager::parameterCount() const {
    if (regionMap_.empty()) {
        if (parameterCount_ == 0){
            return parameterCount_;

            // zero should be possible
            throwLengthError(WHERE_AM_I + " neither region defined nor parameterCount set.");
        }
        return parameterCount_;
    }

    Index count = 0;
    for (auto & x: this->regionMap_){
        // __MS(x.second->parameterCount())
        count += x.second->parameterCount();
    }
    return count;
}

Index RegionManager::constraintCount() const {
    if (regionMap_.empty()) {
        return parameterCount_;
    }
    Index cID = 0;
    for (auto & x: this->regionMap_){
        cID += x.second->constraintCount();
    }
    return cID + interRegionConstraintsCount();
}

Index RegionManager::interRegionConstraintsCount() const {
    Index count = 0;

    for (auto & x: this->interRegionConstraints_){
        std::pair< SIndex, SIndex > ab(x.first);

        if (regionMap_.find(ab.first )->second->isSingle() &&
            regionMap_.find(ab.second)->second->isSingle()){
             count += 1;
        } else {
            if (interRegionInterfaceMap_.find(ab) !=
                interRegionInterfaceMap_.end()){
                count += interRegionInterfaceMap_.find(ab)->second.size();
            }
        }
    }
    return count;
}

void RegionManager::fillConstraints(RSparseMapMatrix & C){
    // __M
    Index nModel  = parameterCount();
    //__MS(nModel)
    Index nConstr = constraintCount();

    this->_cWeights.resize(nConstr, 1.0);
    C.clear();

    //!** no regions: fill 0th-order constraints
    if (regionMap_.empty() || nConstr == 0){
        C.setCols(nModel);
        C.setRows(nModel);
        this->_cWeights.resize(nModel, 1.0);
        for (Index i = 0; i < parameterCount(); i++) {
            C[i][i] = 1.0;
        }
        return;
    }

    C.setRows(nConstr);
    C.setCols(nModel);

    Index cID = 0;

    for (auto & x : this->regionMap_){
        x.second->fillConstraints(C, cID);
        x.second->fillConstraintWeights(this->_cWeights, cID);
        cID += x.second->constraintCount();
        // __MS(cID)
    }

    // __MS(interRegionConstraints_.size())
    if (interRegionConstraints_.size() > 0){
        if (verbose_) std::cout << "Creating inter region constraints." << std::endl;

        // Index i = 0;
        for (auto & it : this->interRegionConstraints_){

            std::pair< SIndex, SIndex > ab = it.first;

            double cWeight = it.second;

            // if (verbose_) {
                // std::cout << "\t" << i << ": "
                //                     << ab.first << "< (" << cWeight << ") >"
                //                     << ab.second << std::endl;
                // i ++;
            // }
            Region * regA = regionMap_.find(ab.first)->second;
            Region * regB = regionMap_.find(ab.second)->second;

            if (regA->isBackground() || regB->isBackground()){
                // __MS(ab.first << " " << ab.second << " " << cWeight)
                // log(Warning, "no inter region constraints for background regions.");
                continue;
            }

            double mcA = regA->modelControl();
            double mcB = regB->modelControl();

            Index aStartParam = regA->startParameter();
            Index bStartParam = regB->startParameter();

            std::map< std::pair< SIndex, SIndex >, std::list < Boundary * > >::const_iterator iRMapIter = interRegionInterfaceMap_.find(ab);

            if (iRMapIter != interRegionInterfaceMap_.end()){

                std::list< Boundary * > bounds = iRMapIter->second;

                for (auto & bIter : bounds){

                    Boundary & b = *bIter;

                    Index aMarker = b.leftCell()->marker();
                    Index bMarker = b.rightCell()->marker();

                    if (aMarker < aStartParam || bMarker < bStartParam){ //** interchange left/right
                        std::swap(aMarker, bMarker);
                    }

                    C[cID][aMarker] = +1.0 / mcA;
                    C[cID][bMarker] = -1.0 / mcB;

                    // setting cWeights
                    if (interfaceConstraints_.count(b.marker())) {
                        // setting interface constraints between regions over inter region constraints ..
                        // inner interfaces are allready mentioned in the regions itself
                        this->_cWeights[cID] = interfaceConstraints_[b.marker()];
                        if (debug()){
                            std::cout << b << " Interface weight: " << interfaceConstraints_[b.marker()] << std::endl;
                        }
                    } else {
                        double zWeight = interRegionConstraintZWeights_;
                        if (std::fabs(zWeight - 1.0) < TOLERANCE){
                            if (debug()){
                                std::cout << b << " Inter region weight: " << cWeight  << std::endl;
                            }
                            this->_cWeights[cID] = cWeight;
                        } else {

                            double zDir = std::fabs(b.norm()[this->mesh_->dim() -1]); //! z-component

                            if (debug()){
                                std::cout << b << " Inter region-z-weight: " << zWeight << std::endl;
                            }
                            this->_cWeights[cID] = cWeight * (1.0 + zDir * (zWeight - 1.0));
                        }
                    }

                    cID ++;

                    if (regA->isSingle() && regB->isSingle()){
                        // only 1 weight so we can break the loop
                        break;
                    }
                } // for each boundary off the A B Interface
            } else {
                log(Error, "No boundaries for the inter region interface.");
            }
        } // for each inter region combination with weight > 0
    } // if have inter regions
}

std::vector < RVector3 > RegionManager::boundaryNorm() const {
    log(Warning, WHERE_AM_I, "Who use this. Is needed?"); //190505
    std::vector < RVector3 > vnorm(constraintCount(), RVector3(0.0, 0.0, 0.0));
    Index cID = 0;

    for (auto & x : this->regionMap_){
        x.second->fillBoundaryNorm(vnorm, cID);
        cID += x.second->constraintCount();
    }
    return vnorm;
}

void RegionManager::loadMap(const std::string & fname){

    std::map< std::string, void (Region::*)(const std::string & val) > regionAttributeMap;
    regionAttributeMap[lower("MC")]       = &Region::setModelControlStr_;
    regionAttributeMap[lower("start")]    = &Region::setStartModelStr_;
    regionAttributeMap[lower("zWeight")]  = &Region::setZWeightStr_;
    regionAttributeMap[lower("fix")]      = &Region::setFixValueStr_;
    regionAttributeMap[lower("Ctype")]    = &Region::setConstraintTypeStr_;
    regionAttributeMap[lower("Trans")]    = &Region::setModelTransStr_;
    regionAttributeMap[lower("lBound")]   = &Region::setLowerBoundStr_;
    regionAttributeMap[lower("uBound")]   = &Region::setUpperBoundStr_;
    regionAttributeMap[lower("single")]   = &Region::setSingleStr_;
    regionAttributeMap[lower("background")] = &Region::setBackgroundStr_;

    if (verbose_) std::cout << "Reading region control file" << std::endl;

    std::fstream file; openInFile(fname, &file);
    std::vector < std::string > token;
    std::vector < std::string > row;

    while (!file.eof()){

        row = getNonEmptyRow(file, '@');
        if (row.empty()){
            continue;
//             file.close();
//             return;
        }

        if (row[0][0] == '#'){

            token.clear();

            if (row[0].size() > 1) {
                // tokenline starts with '#XYZ'
                token.push_back(row[0].substr(1, row[0].size()));
            }

            for (Index i = 1; i < row.size(); i ++) {
                token.push_back(row[i]);
            }

            // read next row line to parse them
            continue;
        }

        if (token.size() == 0) {
             std::cerr << WHERE_AM_I << " not a valid region file. looking for leading '#' in " << fname << std::endl;
             file.close();
             return;
        }

        // interpret row as region informations
        if (lower(token[0]) == "no"){

            if (verbose_){
                if (verbose_) std::cout << "Get region property tokens: " << std::endl;
                for (Index i = 0; i < token.size(); i ++) {
                    if (verbose_) std::cout << token[i] << ", ";
                }
                if (verbose_) std::cout << std::endl;
            }

            if (token.size() >= row.size()){
                IVector regionMarker;

                if (row[0] != "*"){
                    //std::cout << row[0] << std::endl;
                    if (this->regionExists(toInt(row[0]))) {
                        regionMarker.push_back(toInt(row[0]));
                    } else {
                        std::cerr << "Region number " << toInt(row[0]) << " " << row[0] << " does not exist!" << std::endl;
                    }
                } else {
                    regionMarker = allRegionMarker_();
                }

                for (Index j = 0; j < regionMarker.size(); j ++){
                    for (Index i = 1; i < row.size(); i ++){
                        if (regionAttributeMap.count(lower(token[i]))){
                            if (verbose_) {
                                if (j < 5){
                                    std::cout << regionMarker[j] << " : " << token[i]
                                                                 << " " << row[i] << std::endl;
                                } else if (j == 5 && j < regionMarker.size()-1){
                                    std::cout << "..." << std::endl;
                                } else if (j == regionMarker.size()-1){
                                    std::cout << regionMarker[j] << " : " << token[i]
                                                                 << " " << row[i] << std::endl;
                                }

                            }
                            (region(regionMarker[j])->*regionAttributeMap[lower(token[i])])(row[i]);
                        } else {
                            std::cerr << WHERE_AM_I << " no region attribute associated with key:"
                                    << token[i] << std::endl;
                            std::cerr << "valid keys are: " << std::endl;
                            for (std::map< std::string, void (Region::*)(const std::string &) >::iterator it =
                                regionAttributeMap.begin(); it != regionAttributeMap.end(); it ++){
                                std::cerr << "\t" << it->first << std::endl;
                            }
                        }
                    }
                }
            } else {
                std::cerr << WHERE_AM_I << " too few tokens defined in region control file: " << fname << std::endl;
            }
        } else if (lower(token[0]) == "inter-region" ||
                   lower(token[0]) == "interregion") {
            if (verbose_){
                std::cout << "Get inter-region properties" << std::endl;
                for (Index i = 0; i < token.size(); i ++) {
                    std::cout << token[i] << " ";
                }
                std::cout << std::endl;
            }

            if (row.size() == 3){
                IVector minRegion;
                IVector maxRegion;
                if (row[0] != "*"){
                    minRegion.push_back(toInt(row[0]));
                } else {
                    minRegion = allRegionMarker_(true);
                }
                if (row[1] != "*"){
                    maxRegion.push_back(toInt(row[1]));
                } else {
                    maxRegion = allRegionMarker_(true);
                }

                for (Index i = 0; i < minRegion.size(); i ++){
                    for (Index j = 0; j < maxRegion.size(); j ++){
                        if (minRegion[i] != maxRegion[j]){
                            setInterRegionConstraint(minRegion[i],
                                                     maxRegion[j],
                                                    toDouble(row[2]));
                        }
                    }
                }
            } else {
                std::cerr << WHERE_AM_I << " too few tokens defined for inter-region constraints 3 != " << row.size() << std::endl;
            }

        } else if (lower(token[0]) == "interface"){
            if (row.size() == 2){ // interface constraint
                double cw = toDouble(row[1]);
                if (row[0] == "*"){
                    for (auto & x: regionMap_){
                        Region * reg = x.second;
                        for (auto & it: reg->boundaries()){
                            Boundary *b = it;
                            if (b->marker() != 0){
                                log(Debug, "Applying interface constaint: ", b->marker(), "=", cw);
                                this->setInterfaceConstraint(b->marker(), cw);
                            }
                        }
                    }
                } else {
                    SIndex iF = toInt(row[0]);
                    log(Debug, "Applying interface constaint: ", iF, "=", cw);
                    this->setInterfaceConstraint(iF, cw);
                }
            } else {
                for (Index i = 0; i < row.size(); i ++){
                    if (verbose_) std::cout << row[i] << " ";
                }
                log(Error, "Format for interface settings unknown.");
            }

        } else {
            std::cerr << WHERE_AM_I << " cannot interpret 1.st token: " << token[0] << std::endl;
            std::cerr << "Available are: " << std::endl
                      << "#no" << std::endl
                      << "#inter-region" << std::endl
                      << "#interface" << std::endl;
            row = getRow(file); if (row.empty()) continue;
        }
    }
    file.close();

    this->recountParaMarker_();
    this->createParaDomain_();
}

void RegionManager::setInterRegionConstraint(SIndex aIn, SIndex bIn, double cw){
    SIndex a = min(aIn, bIn);
    SIndex b = max(aIn, bIn);

    if (regionMap_.count(a) == 0 || regionMap_.count(b) == 0){
        std::cerr << WHERE_AM_I << " ignoring inter-region constraints (no region)"
                << a << " " << regionMap_.count(a) << " "
                << b << " " << regionMap_.count(b) << std::endl;
        return;
    }
    if (this->region(a)->isBackground() || this->region(b)->isBackground()){
            std::cerr << WHERE_AM_I << " ignoring inter-region constraints (is background)"
                << a << " " << this->region(a)->isBackground() << " "
                << b << " " << this->region(b)->isBackground() << std::endl;
        return;
    }

    if (a == b){
        std::cerr << WHERE_AM_I << " ignoring inter-region constraints "
                << a << " == " << b << std::endl;
    } else {
        if (interRegionInterfaceMap_.find(std::pair< SIndex, SIndex >(a, b))
             != interRegionInterfaceMap_.end()){
            interRegionConstraints_[std::pair< SIndex, SIndex >(a, b)] = cw;
            if (debug()){
               std::cout << "Constraining regions: " << a << "<->" << b
               << "(weigth: " << cw << ")" << std::endl;
           }
        }
    }
}

void RegionManager::saveMap(const std::string & fname){
    THROW_TO_IMPL
}

void RegionManager::setLocalTransFlag(bool flag) {
    haveLocalTrans_ = flag;
    localTransHaveChanges_ = true;
}

TransCumulative < RVector > * RegionManager::transModel(){
    if (regionMap_.empty()){
        return NULL;
    }

    if (localTransHaveChanges_) localTrans_.clear();

    if (localTrans_.size() != allRegionMarker_(true).size()){

        for (auto & x: regionMap_){
            //  for (std::map< SIndex, Region* >::const_iterator
            //  it = regionMap_.begin(); it != regionMap_.end(); it ++){

            if (!x.second->isBackground()){
                if (isPermuted_){
                    // __MS(localTrans_.size() << " " << x.second->paraIds())
                    localTrans_.add(*x.second->transModel(),
                                    x.second->paraIds());
                } else {
                    // __MS(localTrans_.size() << " " << x.second<< " "
                    //         << typeid(*x.second->transModel()).name() << " "
                    //         << x.second->transModel() << " "
                    //         << x.second->startParameter() << " "
                    //         << x.second->endParameter())

                    localTrans_.add(*x.second->transModel(),
                                    x.second->startParameter(),
                                    x.second->endParameter());
                }
            }
        }
    }
    return & localTrans_;
}

void RegionManager::setZWeight(double z){
    for (auto & x: regionMap_){
        x.second->setZWeight(z);
    }
    interRegionConstraintZWeights_ = z;
}

void RegionManager::setConstraintType(Index type){
    for (auto & x: regionMap_){
        x.second->setConstraintType(type);
    }
}

IVector RegionManager::allRegionMarker_(bool excludeBoundary) const {
    IVector tmp;
    for (auto & x: regionMap_){
        if (excludeBoundary && x.second->isBackground()){
            continue;
        }
        tmp.push_back(x.first);
    }
    return tmp;
}

} // namespace GIMLI
