/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *   Thomas Günther thomas@resistivity.net                                    *
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

#ifndef _GIMLI_INVERSION__H
#define _GIMLI_INVERSION__H

// #include "vector.h"
#include "inversionBase.h"
#include "mesh.h"
#include "modellingbase.h"
#include "numericbase.h"
#include "regionManager.h"
#include "solver.h"
#include "stopwatch.h"
#include "sparsematrix.h"
#include "trans.h"
#include "vector.h"


namespace GIMLI{

#define DOSAVE if (dosave_)

//#define PLUS_TMP_VECSUFFIX + ".vec"
#define PLUS_TMP_VECSUFFIX

/*! template function for computing L1 norm (robust/blocky) weightings */
template < class Vec > Vec getIRLSWeights(const Vec & a, double locut = 0.0, double hicut = 0.0) {
    double suabs = sum(abs(a));
    double suabsq = dot(a, a);

    Vec tmp(suabsq / suabs / (abs(a) + TOLERANCE));
    for(uint i = 0; i < a.size(); i++) {
        if((locut > 0.0) && (tmp[ i ] < locut)) tmp[ i ] = locut;
        if((hicut > 0.0) && (tmp[ i ] > hicut)) tmp[ i ] = hicut;
    }
    return tmp;
}

/*! template function for computing Lp norm (robust/blocky) weightings */
template < class Vec > Vec getIRLSWeightsP(const Vec & a, int p, double locut = 0.0, double hicut = 0.0) {
    Vec ap = pow(abs(a), p);
    Vec aq = pow(abs(a), 2);
    double suabsp = sum(ap);
    double suabsq = sum(aq);

    Vec tmp((ap + 1e-12) / (aq + 1e-12) * suabsp / suabsq);
    for(uint i = 0; i < a.size(); i++) {
        if((locut > 0.0) && (tmp[ i ] < locut)) tmp[ i ] = locut;
        if((hicut > 0.0) && (tmp[ i ] > hicut)) tmp[ i ] = hicut;
    }
    return tmp;
}

/*! Inversion template using a Generalized Minimization Approach, atm fixed to Gauss-Newton solver
    Inversion(bool verbose, bool dosave
    Inversion(RVector data, FOP f, bool verbose, bool dosave
    Inversion(RVector data, FOP f, transData, transModel, bool verbose, bool dosave */
class DLLEXPORT RInversion : public InversionBase< double > {
public:
    typedef RVector Vec;
    typedef RVector ModelVector;

    /*! Minimal constructor. verbose -- gives some output and status information about the inversion progress; dosave -- advanced debuging, saves alot of stuff and temporary data */
    RInversion(bool verbose=false, bool dosave=false)
        : InversionBase< double >(), verbose_(verbose),
        dosave_(dosave), saveModelHistory_(dosave){
        this->init_();
    }

    /*! Contructor for forward operator only */
    RInversion(ModellingBase & forward,
              bool verbose=false, bool dosave=false)
        : InversionBase< double >(), verbose_(verbose),
        dosave_(dosave), saveModelHistory_(dosave) {
        this->init_();

        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);
    }

    /*! Usual constructor. data -- vector of the given data; forward -- the associated forward operator */
    RInversion(const Vec & data, ModellingBase & forward,
              bool verbose=false, bool dosave=false)
    : InversionBase< double >(), verbose_(verbose),
        dosave_(dosave), saveModelHistory_(dosave) {
        this->init_();

        //** init: data_
        this->setData(data);

        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);
    }

    /*! Intermediate constructor including transforms and data transformation function */
    RInversion(const Vec & data, ModellingBase & forward,
              Trans< Vec > & transData, bool verbose=true, bool dosave=false)
    : InversionBase< double >(), verbose_(verbose),
        dosave_(dosave), saveModelHistory_(dosave) {
        //** set: default values
        this->init_();
        //** set: paraDomain init: modelWeight_, constraintWeights_, ConstraintMatrix
        //** init: data_
        this->setData(data);
        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);

        this->setTransData(transData);
    }

    /*! Full constructor including transforms. transData -- data transformation function; transModel -- model transformation function */
    RInversion(const Vec & data, ModellingBase & forward,
              Trans< Vec > & transData, Trans< Vec > & transModel,
              bool verbose = true, bool dosave = false)
    : InversionBase< double >(), verbose_(verbose), dosave_(dosave), saveModelHistory_(dosave) {
        //** set: default values
        this->init_();
        //** set: paraDomain init: modelWeight_, constraintWeights_, ConstraintMatrix
        //** init: data_
        this->setData(data);
        //** init: model_, modelRef, response_
        this->setForwardOperator(forward);

        this->setTransData(transData);
        this->setTransModel(transModel);
    }

    /*! Destructor. Frees allocated memory */
    virtual ~RInversion(){
        delete transDataDefault_;
        delete transModelDefault_;
    }

private:
    /*! Copyconstructor */
    RInversion(const RInversion & inv){
        THROW_TO_IMPL
    }


protected:
    /*! Internal initialization function, which is called from constructor. Set default paramater and allocate required memory */
    void init_(){
        transDataDefault_   = new Trans< Vec >;
        transModelDefault_  = new Trans< Vec >;
        tD_                 = transDataDefault_;
        tM_                 = transModelDefault_;
        isRobust_           = false;
        isBlocky_           = false;
        useLinesearch_      = true;
        optimizeLambda_     = false;
        recalcJacobian_     = true;
        jacobiNeedRecalc_   = true;
        doBroydenUpdate_    = false;
        localRegularization_= false;
        abort_              = false;
        stopAtChi1_         = true;
        haveReferenceModel_ = false;
        isRunning_          = false;
        fixError_           = true;
        activateFillConstraintWeights_ = true; //jointinv hack!!!

        iter_               = 0;
        maxiter_            = 20;
        maxCGLSIter_        = 200;
        lambda_             = 50.0;
        lambdaFactor_       = 1.0;
        lambdaMin_          = 1.0;
        dPhiAbortPercent_   = 2.0;

        CGLStol_            = -1.0; //** -1 means automatic scaled
    }

public:

    /*! Set data vector */
    inline void setData(const Vec & data) { data_ = data; }

    /*! Get data vector */
    inline const Vec & data() const { return data_; }

    /*! Set the relative data error to the vector err. */
    void setRelativeError(const Vec & e){
        error_ = e;
        checkError();
    }

    /*! Set relative data error to scalar error value relerr.
     * If you force relerr == 0 here. Data will not corrected or weighted. */
    inline void setRelativeError(double relerr) {
        if (relerr == 0.0) fixError_ = false;
        setRelativeError(RVector(data_.size(), relerr));
    }

    /*! Set the absolute data error to the vector abserr */
    inline void setAbsoluteError(const Vec & abserr) {
        setRelativeError(abs(abserr) / abs(fixZero(data_, TOLERANCE)));
    }

    /*! Set absolute data error to the scalar value abserr */
    inline void setAbsoluteError(double abserr) {
        setRelativeError(std::fabs(abserr) / abs(fixZero(data_, TOLERANCE)));
    }

    inline void setError(const Vec & err, bool isRelative=true) {
        if (isRelative){
            this->setRelativeError(err);
        } else {
            this->setAbsoluteError(err);
        }
    }
    inline void setError(double err, bool isRelative=true) {
        if (isRelative){
            this->setRelativeError(err);
        } else {
            this->setAbsoluteError(err);
        }
    }
    /*! Return the used data error */
    inline const Vec & error() const { return error_; }

    /*! Validate the data error. This check will be called at every new run or change of the error;
        Check size and content unequal to zero. If problems found, warn and set the data error to errorDefault_()  */
    void checkError(){
        if (error_.size() != data_.size()) {
            std::cerr << WHERE_AM_I << " Warning error has the wrong size, reset to default. "
                      << error_.size() << " != " << data_.size() << std::endl;
            error_ = errorDefault_();
        } else if (min(abs(error_)) < TOLERANCE && fixError_) {
            std::cerr << WHERE_AM_I << " Warning error contains zero values, reset to default. " << std::endl;
            error_ = errorDefault_();
        }

        dataWeight_ = 1.0 / tD_->error(fixZero(data_, TOLERANCE), error_);

        if (verbose_) std::cout << "min/max(dweight) = " << min(dataWeight_) << "/"
                                  << max(dataWeight_) << std::endl;

        if (haveInfNaN(dataWeight_)){
            DOSAVE save(dataWeight_,          "Nan_dataWeight_dweight");
            DOSAVE save(error_,               "Nan_dataWeight_error");
            DOSAVE save(data_,                 "Nan_dataWeight_data");

            throwError(WHERE_AM_I + " dataWeight_ contains inf or nan");
        }
    }

    /*! Return a copy of the default data error array, [ Relative data error of 1% ] */
    inline Vec errorDefault_(){ return Vec(data_.size(), 0.01); }

    /*! Set the forward operator and the model-transform-function derived from fop-regionManager if not set before. */
    inline void setForwardOperator(ModellingBase & forward) {
        forward_   = & forward;
        forward_->clearConstraints();  //! why is this so strictly necessary???
        //! Always use a region manager
        forward_->initRegionManager();
//         model_ = forward_->startModel();
    }
    /*! Return forward operator.*/
    inline ModellingBase * forwardOperator() { return forward_; }
    /*! Return forward operator. Shortcut for forwardOperator() */
    inline ModellingBase * fop() { return this->forwardOperator(); }

    /*TODO Change interface to Trans < vec > *tD, or copy the trans ,giving non const reference is dangerous.*/
    /*! Set and get data transform.
     * WARNING! we just store a reference to the trans .. U have to ensure that the trans is and stayes valid.
     * TODO change this!!!
     */
    inline void setTransData(Trans< Vec > & tD) { tD_ = & tD; }
    inline Trans< Vec > & transData() { return * tD_; }

    /*TODO Change interface to Trans < vec > *tD, or copy the trans ,giving non const reference is dangerous. */
    /*! Set and get Model transform
     * WARNING! we just store a reference to the trans .. U have to ensure that the trans is and stayes valid.
     * TODO change this!!!
     */
    inline void setTransModel(Trans< Vec > & tM) { tM_ = & tM; }
    inline Trans< Vec > & transModel() { return * tM_; }

    /*! Return number of constraints */
    inline uint constraintsCount() const { return constraintWeights_.size(); }

    /*! Return number of data */
    inline uint dataCount() const { return data_.size(); }

    /*! Set verbose output */
    inline void setVerbose(bool verbose){ verbose_ = verbose; }
    inline bool verbose() const { return verbose_; }

    /*! Set debug output */
    inline void setDoSave(bool d ){ dosave_ = d; }
    inline bool doSave() const { return dosave_; }

    /*! Set and get verbose behaviour */
    inline void saveModelHistory(bool doSaveModelHistory){ 
        saveModelHistory_ = doSaveModelHistory; }

    /*! Set and get line search */
    inline void setLineSearch(bool linesearch) { useLinesearch_ = linesearch; }
    inline bool lineSearch() const { return useLinesearch_; }

    /*! Set and get blocky model behaviour (by L1 reweighting of constraints) */
    inline void setBlockyModel(bool isBlocky) { isBlocky_ = isBlocky; }
    inline bool blockyModel() const { return isBlocky_; }

    /*! Set and get robust data weighting by L1 reweighting scheme */
    inline void setRobustData(bool isRobust) { isRobust_ = isRobust; }
    inline bool robustData() const { return isRobust_; }

    /*! Set and get regularization parameter and its change over iteration */
    inline void setLambda(double lambda) { lambda_ = lambda; }
    inline double lambda() const { return lambda_; }
    inline double getLambda() const { return lambda_; }
    inline void setLambdaFactor(double lambdaFactor) { lambdaFactor_ = lambdaFactor; }
    inline double lambdaFactor() const { return lambdaFactor_; }

    /*! Set the minimum lambda value that can be reached if lambda factor is set.
     *Default is 1. */
    inline void setLambdaMinimum(double l) { lambdaMin_ = l; }
    /*! Return the minimum lambda value that can be reached if lambda factor is set. */
    inline double lambdaMinimum() const { return lambdaMin_; }

    /*! Set whether regularization is global or local (e.g. Marquardt method) */
    void setLocalRegularization(bool localReg){ 
        localRegularization_ = localReg; }
    /*! Return whether regularization is global or local (e.g. Marquardt method) */
    bool localRegularization() const { return localRegularization_; }

    /*! Set and get whether lambda is being optimized by Lcurve */
    inline void setOptimizeLambda(bool opt) { optimizeLambda_ = opt; }
    inline bool optimizeLambda() const { return optimizeLambda_; }

    /*! Set and get maximum iteration number */
    inline void setMaxIter(int maxiter) { maxiter_ = maxiter; }
    inline int maxIter() const { return maxiter_; }

    /*! Set and get maximum iteration number for inverse sub step*/
    inline void setMaxCGLSIter(int iter){ maxCGLSIter_ = iter; }
    inline int maxCGLSIter() const { return maxCGLSIter_; }

    /*! Set and get abort tolerance for cgls solver. -1 means default scaling [default] */
    inline void setCGLSTolerance(double tol){ CGLStol_ = tol; }
    inline double maxCGLSTolerance() const { return CGLStol_; }

    /*! Return curent iteration number */
    inline uint iter() const { return iter_; }

    /*! Return true if inversion is running */
    inline bool isRunning() const { return isRunning_ ; }

    /*! Force abort at next iteration */
    inline void abort() { abort_ = true; }

    /*! Define whether minimization stops at chi^2=1 */
    inline void stopAtChi1(bool stopAtChi1) { stopAtChi1_ = stopAtChi1; }

    /*! Set and get relative objective function decrease */
    inline void setDeltaPhiAbortPercent(double dPhi){ dPhiAbortPercent_ = dPhi; }
    inline double deltaPhiAbortPercent() const { return dPhiAbortPercent_; }

    /*! Marquardt scheme (local damping with decreasing regularization strength */
    void setMarquardtScheme(double lambdaFactor=0.8){
        setLocalRegularization(true);  //! no contribution of regularization to objective function
        setLambdaFactor(lambdaFactor); //! lambda is decreased towards zero
        stopAtChi1(false);             //! let the solution fully converge (important for statistics)
        if (forward_) forward_->regionManager().setConstraintType(0);
    }

    /*! Check if transfunctions are valid, set default if no transfunctions are given or override given functions by regionManager.transfunctions if available*/
    void checkTransFunctions(){
        if (forward_->regionManager().haveLocalTrans()) {
            if (verbose_) std::cout << "use model trans from RegionManager" << std::endl;
            setTransModel(*forward_->regionManager().transModel());
        }
    }


    /*! Define and find out whether Jacobian is recalculated in each iteration */
    void setRecalcJacobian(bool recalc){
        recalcJacobian_ = recalc;
        if (recalc) doBroydenUpdate_ = false;
    }
    bool recalcJacobian() const { return recalcJacobian_; }

    /*! Enable/disable broyden update (which enforces scaling by transform function derivatives */
    void setBroydenUpdate(bool broydenUpdate){
        doBroydenUpdate_ = broydenUpdate;
        if (doBroydenUpdate_) recalcJacobian_ = false;
        //** puts transformation functions into jacobian since it is not changed anymore;
        if (broydenUpdate) {
            THROW_TO_IMPL
            //scaleMatrix(*J_, tD_->deriv(response_), Vec(1.0 / tM_->deriv(model_)));
        }
    }
    inline bool broydenUpdate() const { return doBroydenUpdate_; }

    /*! Set model vector .
     * If you call \ref run() the inversion starts with this model,
     * otherwise it will start with fop.startModel(). */
    void setModel(const Vec & model){
        if (recalcJacobian_ && model != model_) jacobiNeedRecalc_ = true;
        model_ = model;
    }

    /*! Return a const reference to the current model vector */
    inline const ModelVector & model() const { return model_; }

    /*! Set reference model vector */
    inline void setReferenceModel(const Vec & model){
        haveReferenceModel_ = true; modelRef_ = model;
    }

    /*! Set current model response (e.g. to avoid time-consuming forward calculation, but be careful) */
    inline void setResponse(const Vec & response){ response_ = response; }

    /*! Set constraint right-hand side by hand (very special cases, so be careful) */
    inline void setConstraintsH(const Vec & constraintsH){ 
        __MS("who use this. Please note any setting of this will be overwritten in run.")
        constraintsH_ = constraintsH; 
    }  // size check?

    /*! Return a reference to the current response vector */
    inline const Vec & response() const { return response_; }

    /*! Return IRLS function of roughness vector */
    Vec getIRLS() const {
        return getIRLSWeights(Vec(*forward_->constraints() * tM_->trans(model_) * constraintWeights_), 0.0, 1.0);
    }

    /*! Set the constraint weight (boundary control) vector */
    void setCWeight(const Vec & cWeight){
        constraintWeights_ = cWeight;
        activateFillConstraintWeights_ = false; //jointinv hack!!!
        if (verbose_) std::cout << "min/max(cWeight) = " << min(constraintWeights_) << "/" << max(constraintWeights_) << std::endl;
    }

    /*! Return reference to the current constraints weight vector */
    inline const Vec & cWeight() const { return constraintWeights_; }

    /*! Set the model weight vector */
    void setMWeight(const Vec & mweight){
        modelWeight_ = mweight;
        if (verbose_) std::cout << "min/max(mWeight) = " << min(modelWeight_) << "/" << max(modelWeight_) << std::endl;
    }

    /*! Return reference to the current model weight vector */
    inline const Vec & mWeight() const { return modelWeight_; }

    /*! Call robust date reweighting (see also setRobustData) */
    void robustWeighting() {
        if (verbose_) std::cout << "Robust reweighting " << std::endl;
        Vec deltaData((tD_->trans(data_) - tD_->trans(response_)) * dataWeight_);
        //** reweight errors according to IRLS scheme
        // cr: error wird nur noch von getPhiD genutzt und dataweight wird im Laufe der Inversion nicht mehr aktualisiert, soll das so sein? DRY?
        // tg: checkerror ist doch genau das richtige oder?
        error_ /= (getIRLSWeights(deltaData, 0.0, 1.0) + TOLERANCE);
        checkError(); //** macht implizit:
        // dataWeight_ = (1.0 / tD_->error(data_, error_));
        // tg: maybe better just reweight dweight and hold error constant
        // tg: geht nicht weil phiD auf Fehlern beruht (und das auch muss wegen Trafowechsel)
    }

    /*! Apply blocky model constraints (see also setBlockyModel) */
    void constrainBlocky() {
        if (verbose_) std::cout << "Blocky model constraints " << std::endl;
        setCWeight(getIRLSWeights(Vec((*forward_->constraints() * tM_->trans(model_)) * constraintWeights_), 0.0, 1.0));
    }

    /*! Create constraints, check and compare size of constraint matrix with model/boundary control */
    void checkConstraints();

    /*! Check size of Jacobian matrix against data and model number */
    void checkJacobian(bool force=false);

    /*! Shortcut for \ref echoStatus(response_, model_). */
    inline void echoStatus() const { echoStatus(response_, model_); }

    /*! Echo PhiD/PhiM/Phi and Chi^2 for given model and response */
    void echoStatus(const Vec & response, const Vec & model, const std::string & xtra = "") const {
        double chi2 = getPhiD(response) / data_.size();
        std::cout << iter_ << ": " << xtra
                  << "Model: min = " << min(model) << "; max = " << max(model) << std::endl;
        std::cout << iter_ << ": " << xtra
                  << "Response: min = " << min(response) << "; max = " << max(response) << std::endl;
        std::cout << iter_ << ": rms/rrms(data, " << xtra << "Response) = "
                        << rms(data_, response) << "/"
                        << rrms(data_, response) * 100.0 << "%" << std::endl;
        std::cout << iter_ << ": chi^2(data, " << xtra
                  << "Response, error, log) = " << chi2 <<        std::endl;
        std::cout << iter_ << ": Phi = " << getPhiD(response) << "+" << getPhiM(model)
                            << "*" << lambda_ << "=" << getPhi(model, response) << std::endl;
    }

    /*! Compute model cell resolution (specific column of resolution matrix) by an LSCG solver (Guenther, 2004) */
    Vec modelCellResolution(int iModel) {
        Vec resolution(model_.size(), 0.0);
        resolution[ iModel ] = 1.0; //** both for retrieving solumn and as starting model for resolution
        //** retrieve one colomn from the jacobian matrix and scale according to data/model transformation
        Vec sensCol = (*forward_->jacobian()) * resolution * tD_->deriv(response_) / tM_->deriv(model_)[ iModel ];
        //** call inverse substep with sensCol on right hand side (see Gnther, 2004)
        Vec deltaModel0(model_.size());// !!! h   zvariante

        solveCGLSCDWWtrans(*forward_->jacobian(), *forward_->constraints(),
                           dataWeight_, sensCol, resolution,
                           constraintWeights_, modelWeight_,
                           tM_->deriv(model_), tD_->deriv(response_),
                           lambda_, deltaModel0, maxCGLSIter_, false);

        return resolution;
    }

    /*! Compute cell resolution for closest cell to given position */
    Vec modelCellResolution(const RVector3 & pos){
        Index ipos = 0;
        double mindist = 9e99;
        R3Vector cellCenters = forward_->mesh()->cellCenter();

        for (Index i = 0 ; i < model_.size() ; i++){
            double dist = cellCenters[i].distance(pos);
            if (dist < mindist) {
                mindist = dist;
                ipos = i;
            }
        }
        // Achtung Parametermapping und Backgroundzellen!
        return modelCellResolution(ipos);
    }

    /*! Compute whole resolution matrix by single cell resolutions */
    RMatrix modelResolutionMatrix() {
        RMatrix resMat;
        for (size_t i = 0 ; i < model_.size(); i++) resMat.push_back(modelCellResolution(i));
        return resMat;
    }

    /*! Return (C * m * m_w) * c_w */
    RVector roughness(const RVector & model) const {
       RVector r(*forward_->constraints() 
                  * Vec(tM_->trans(model) * modelWeight_) 
                  * constraintWeights_);

        if (haveReferenceModel_) {
            r = r - constraintsH_;
        }
        return r;
    }
    /*! Shortcut for roughness for the current model vector */
    RVector roughness() const {
       return roughness(model_);
    }
    /*! Return (C * m) , i.e. the pure (unweighted) roughness */
    RVector pureRoughness(const RVector & model) const {
       return *forward_->constraints() * Vec(tM_->trans(model));
    }

    /*! Return data objective function (sum of squared data-weighted misfit) */
    double getPhiD(const Vec & response) const;

    /*! Return model objective function (squared model roughness) */
    double getPhiM(const Vec & model) const;

    /*! Return total objective function (data OF plus lambda times model OF)
     DEPRECATED wrong nameing scheme*/
    double getPhi(const Vec & model, const Vec & response) const {
        return getPhiD(response) + getPhiM(model) * lambda_ * (1.0 - double(localRegularization_));
    }

    /*! Shortcut for \ref getPhiD(response_) necessary ?
     DEPRECATED wrong nameing scheme*/
    inline double getPhiD() const { return getPhiD(response_); }

    /*! Shortcut for \ref getPhiM(model_) necessary ?
     DEPRECATED wrong nameing scheme*/
    inline double getPhiM() const { return getPhiM(model_); }

    /*! Shortcut for \ref getPhi(model_, response_) necessary ?
     DEPRECATED wrong nameing scheme*/
    inline double getPhi() const { return getPhiD() + getPhiM() * lambda_ * (1.0 - double(localRegularization_)); }

    /*! Return the chi-squared data misfit
     DEPRECATED wrong nameing scheme*/
    inline double getChi2() const { return getPhiD() / data_.size(); }

    /*! Return the chi-squared data misfit for given forward response.
     DEPRECATED wrong nameing scheme*/
    inline double getChi2(const Vec & response) const { return getPhiD(response) / data_.size(); }

    /*! Return last chi-squared data misfit */
    inline double chi2() const { return getPhiD() / data_.size(); }

    /*! Return last absolute RMS misfit */
    inline double absrms() const { return rms(data_, response_); }

    /*! Return last relative RMS misfit */
    inline double relrms() const { return rrms(data_, response_) * 100.; }

    /*! Start with linear interpolation, followed by quadratic fit if linesearch parameter tau is lower than 0.03. Tries to return values between 0.03 and 1 */
    double linesearch(const Vec & modelNew, const Vec & responseNew) const {
        Vec phiVector(101, getPhi());
        Vec phiDVector(101 , getPhiD());

        Vec dModel(tM_->trans(modelNew)    - tM_->trans(model_));
        Vec dData( tD_->trans(responseNew) - tD_->trans(response_));

        double tau = 0.0, minTau = 0.0;
        double minPhi = phiVector[ 0 ];

        if (localRegularization_) minPhi = phiDVector[ 0 ];

        double thisPhi = minPhi;
        for (int i = 1; i < 101; i++) {
            tau = 0.01 * (double) i;
            Vec appModel (tM_->update(model_, dModel * tau));
            Vec appResponse (tD_->update(response_, dData * tau));
            phiVector[ i ]  = getPhi(appModel, appResponse);
            phiDVector[ i ] = getPhiD(appResponse);
            thisPhi = phiVector[ i ]; //! rigorous minimization of the total objective function
            //** this could also be controlled by another switch (which is enforced by local reg.)
            if (localRegularization_) thisPhi = phiDVector[ i ];
            if (thisPhi < minPhi){
                minPhi = thisPhi;
                minTau = tau;
            }
        }
        DOSAVE save(phiVector,  "linesearchPhi");
        DOSAVE save(phiDVector, "linesearchPhiD");

        tau = minTau;

        /*! parabolic line search using step 0, 1 and tauquad to fit parabola */
        if (tau < 0.03) {
            double tauquad = 0.3;
            if (verbose_) std::cout << "tau = " << tau
                            << ". Trying parabolic line search with step length " << tauquad;
            RVector modelQuad(tM_->update(model_, dModel * tauquad));
            RVector responseQuad (forward_->response(modelQuad));
            tau = linesearchQuad(modelNew, responseNew, modelQuad, responseQuad, tauquad);
            if (verbose_) std::cout << " ==> tau = " << tau;
            if (tau > 1.0) { //! too large
                tau = 1.0;
                if (verbose_) std::cout << " resetting to " << tau;
            }
            if (verbose_) std::cout << std::endl;

            if (tau < 0.03) { //! negative or nearly zero (stucked) -> use small step instead
                tau = 0.03;
                if (verbose_) std::cout << " tau < 0.03 ==> tau = " << tau << std::endl;
            }
        } // else tau > 0.03

        if (tau < 0.95) { //! only if not nearly one (which saves a forward run
            if (verbose_) echoStatus(responseNew, modelNew, "LS new");
        }

        if (verbose_) std::cout << "Performing line search with tau = " << tau << std::endl;

        return tau;
    }

    /*! Compute objective function for old (tau=0), new (tau=1) and another model */
    double linesearchQuad(const Vec & modelNew, const Vec & responseNew,
                           const Vec & modelQuad, const Vec & responseQuad,
                           double tauquad) const {
        double phi0  = getPhi();
        double phi10 = getPhi(modelNew, responseNew)   - phi0;
        double phit0 = getPhi(modelQuad, responseQuad) - phi0;
        double dphit = phit0 - phi10 * tauquad;
        if (abs(dphit) < TOLERANCE) return 0.0;
        double tauopt = (phit0 - phi10 * tauquad * tauquad) / dphit / 2.0;

        DOSAVE std::cout << "LineSearchQuad: Phi = " << phi0 << " - " << phit0 + phi0
                         << " - " << phi10 + phi0 << " -> tau= " << tauopt << std::endl;
        return tauopt;
    }

    /*! Return the single models for each iteration. For debugging.*/
    inline const std::vector < RVector > & modelHistory() const { return modelHist_; }

    /*! Compute model update by solving one inverse sub-step
       \param rhs The right-hand-side vector of the system of equation
    */
    Vec invSubStep(const Vec & rhs) {
        Vec deltaModel0(model_.size());//!!! h-variante
        Vec solution(model_.size());
        solveCGLSCDWWtrans(*forward_->jacobian(), *forward_->constraints(),
                           dataWeight_, rhs, solution, constraintWeights_,
                           modelWeight_,
                           tM_->deriv(model_), tD_->deriv(response_),
                           lambda_, deltaModel0, maxCGLSIter_, dosave_);
        return solution;
    }

    /*! Optimization of regularization parameter by L-curve */
    Vec optLambda(const Vec & deltaData, const Vec & deltaModel0); //!!! h-variante

    /*! One iteration step. Return true if the step can be calculated successfully else false is returned. */
    bool oneStep();

    /*! Start the inversion with specific data.*/
    const Vec & invert(const Vec & data);

    /*! Start the inversion procedure from starting model.*/
    const Vec & start();

    /*! Start the inversion procedure from the last model and return the final model vector.*/
    const Vec & run();

    /*! Specialized run function that tries to reach a datafit chi^2=1 by varying the regularization paramater lambda */
    Vec runChi1(double acc = 0.01, int maxiter = 50){
        stopAtChi1(false);
        Vec model = run();

        double lambda = lambda_;
        double chi2 = getChi2();
        double fak = 2.0, oldchi2 = chi2;
        double dir = 0, olddir = 0;
        bool verbose = verbose_;
//        setVerbose(false);
        forward_->setVerbose(false);
        if (verbose) std::cout << "Optimizing lambda subject to chi^2=1." << std::endl;
        if (verbose) std::cout << "chi^2 = " << chi2 << " lambda = " << lambda << " dir = " << dir << std::endl;
        int iter = 0;
        while (std::fabs(chi2 - 1.0) > acc and iter < maxiter){
            if (dir < 0 && chi2 > oldchi2*1.001) break;
            dir = - sign(chi2 - 1.0);                           //** direction: up (1) or down (-1)
            if (dir * olddir == -1) fak = std::pow(fak, 0.6); //** change direction: decrease step
            lambda *= std::pow(fak, dir);                       //** increase or decrease lambda
            setLambda(lambda);
            model = run();
            chi2 = getChi2();
            if(verbose) std::cout << "chi^2 = " << chi2 << " lambda = " << lambda << " dir = " << dir << std::endl;
            olddir = dir;                                         //** save old direction for step length
            oldchi2 = chi2;
            iter++;
        }
        return model;
    }

    const RVector & modelWeight() const { return modelWeight_; }
    const RVector & modelRef() const { return modelRef_; }
    const RVector & dataWeight() const { return dataWeight_; }

    const RVector & deltaDataIter() const { return deltaDataIter_; }
    const RVector & deltaModelIter() const { return deltaModelIter_; }

    /*! Resets this inversion to the given startmodel. */
    void reset(){
        this->setModel(forward_->startModel());
    }

protected:

    Vec                   data_;
    ModellingBase       * forward_;

    Trans< Vec >        * tD_;
    Trans< Vec >        * tM_;
    Trans< Vec >        * transModelDefault_;
    Trans< Vec >        * transDataDefault_;

    bool verbose_;
    bool dosave_;
    bool saveModelHistory_;

    Vec  error_;
    Vec  response_;

    Vec  model_;
    Vec  modelRef_;

    Vec  constraintsH_;
    Vec  constraintWeights_;
    Vec  modelWeight_;
    Vec  dataWeight_;

    Vec  deltaDataIter_;
    Vec  deltaModelIter_;

    int maxiter_;
    int iter_;
    int maxCGLSIter_;

    double lambda_;
    double lambdaFactor_;
    double lambdaMin_;
    double dPhiAbortPercent_;
    double CGLStol_;

    bool isBlocky_;
    bool isRobust_;
    bool isRunning_;
    bool useLinesearch_;
    bool optimizeLambda_;
    bool abort_;
    bool stopAtChi1_;
    bool doBroydenUpdate_;
    bool localRegularization_;
    bool haveReferenceModel_;
    bool recalcJacobian_;
    bool jacobiNeedRecalc_;
    bool activateFillConstraintWeights_; //jointinv hack!!!

    /*! Set this to zero if u want to use absolute errors == zero*/
    bool fixError_;

    /*! Hold old models, for debuging */
    std::vector < RVector > modelHist_;

};

} // namespace GIMLI

#endif
