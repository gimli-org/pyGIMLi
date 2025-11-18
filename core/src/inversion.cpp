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

#include "inversion.h"

namespace GIMLI{

void RInversion::checkConstraints() {
    if (forward_->constraints()->cols() == 0 ||
        forward_->constraints()->rows() == 0){
        if (verbose_) std::cout << "Building constraints matrix" << std::endl;
        //forward_->regionManager().fillConstraints(forward_->constraints());
        forward_->createConstraints();
    } else {
        if (verbose_) std::cout << " found valid constraints matrix. omit rebuild" << std::endl;
    }
    Index nModelC = forward_->constraints()->cols();
    Index nCWeightC = forward_->constraints()->rows();

    if (verbose_){
        std::cout << "constraint matrix of size(nBounds x nModel) "
                  << nCWeightC << " x " << nModelC << std::endl;
    }

    if (model_.size() != nModelC){
        //** vielleicht irgendwo ne schöne resize funktion
        std::cout << WHERE_AM_I << " resize model " << model_.size()
                    << " to fit constrain size: "
                    << nModelC << std::endl;
        model_.resize(nModelC);
    }

    forward_->regionManager().fillModelControl(modelWeight_);

    if (modelWeight_.size() != model_.size()){
        modelWeight_.resize(model_.size(), 1.0);
    }
    if (activateFillConstraintWeights_) {
        constraintWeights_ = forward_->regionManager().constraintWeights();
    }
    if (constraintWeights_.size() != nCWeightC) {
        constraintWeights_.resize(nCWeightC, 1.0);
    }
}

void RInversion::checkJacobian(bool force) {

    if ((forward_->jacobian()->rows() == data_.size() &&
        forward_->jacobian()->cols() == model_.size()) && !force) return;

    if (verbose_ &&
            (forward_->jacobian()->rows() != data_.size() ||
                forward_->jacobian()->cols() != model_.size())
        ){
        std::cout << "check Jacobian: wrong dimensions: "
                    << "(" << forward_->jacobian()->rows()  << "x" << forward_->jacobian()->cols() << ") should be "
                    << "(" << data_.size() << "x" << model_.size()  << ") " << " force: " << force << std::endl;
        std::cout << "jacobian size invalid, forced recalc: " << force << std::endl;
    }

    Stopwatch swatch(true);
    if (verbose_) std::cout << "Calculating Jacobian matrix (forced=" << force << ")...";
    forward_->createJacobian(model_);
    jacobiNeedRecalc_ = false;
    if (verbose_) std::cout << "... " << swatch.duration(true) << " s" << std::endl;
}

double RInversion::getPhiD(const Vec & response) const {
    Vec deltaData((tD_->trans(data_) - tD_->trans(response)) /
                    tD_->error(fixZero(data_, TOLERANCE), error_)) ;

    double ret = dot(deltaData, deltaData);
    if (isnan(ret) || isinf(ret)){
        save(tD_->trans(data_),          "Nan_PhiD_tD_data");
        save(response,                     "Nan_PhiD_response");
        save(tD_->trans(response),       "Nan_PhiD_tD_response");
        save(tD_->error(data_, error_), "Nan_PhiD_tD_error");

        throwError(WHERE_AM_I + " getPhiD == " + str(ret));
    }
    return ret;
}

double RInversion::getPhiM(const Vec & model) const {
//        Vec dModel(tM_->trans(model));
//        if (haveReferenceModel_) dModel = dModel - tM_->trans(modelRef_);
//        Vec roughness(Vec(forward_->constraints() * dModel) * constraintWeights_);
    Vec rough(this->roughness(model));

    double ret = dot(rough, rough);
    if (isnan(ret) || isinf(ret)){
        DOSAVE std::cerr << "haveReferenceModel_: " << haveReferenceModel_<< std::endl;
        DOSAVE save(model,       "Nan_PhiM_model");
        DOSAVE save(modelRef_,  "Nan_PhiM_modelref");
//            DOSAVE save(dModel,      "Nan_PhiM_tM_dmodel");
        DOSAVE save(rough,   "Nan_PhiM_roughness");
        DOSAVE save(constraintWeights_,   "Nan_PhiM_cweight");

        throwError(WHERE_AM_I + " getPhiM == " + str(ret));
    }
    return ret;
}

const RVector & RInversion::invert(const RVector & data){
    this->reset();
    this->setData(data);
    return run();
}

/*! Start inversion from starting model. */
const RVector & RInversion::start(){
    this->reset();
    return run();
}

/*! Run inversion with current model. */
const RVector & RInversion::run(){ ALLOW_PYTHON_THREADS

    if (model_.size() == 0) setModel(forward_->startModel());

    if (data_.size() == 0) {
        throwError(WHERE_AM_I + " no data given");
    }

    abort_ = false;

    //** check if transfunctions are valid
    this->checkTransFunctions();

    //! calculation of initial modelresponse
    response_ = forward_->response(model_);
    //response_ = forward_->response(forward_->startModel());

    //! () clear the model history
    modelHist_.clear();

    //! validate and rebuild the data error if necessary
    this->checkError();

    //! validate and rebuild the constraint matrix if necessary
    this->checkConstraints();

    if (haveReferenceModel_) {
        //! compute roughness constraint and correct it for inter-region constraints
        Index cc = forward_->regionManager().constraintCount();
        constraintsH_ = (*forward_->constraints() *
                          Vec(tM_->trans(modelRef_) * modelWeight_)
                        ) * constraintWeights_;
        Index ircc = forward_->regionManager().interRegionConstraintsCount();
        if (ircc > 0) constraintsH_.setVal(0.0, cc - ircc, (SIndex)cc);
    }

    //! validate and rebuild the jacobian if necessary
    this->checkJacobian(jacobiNeedRecalc_);

    //** End preparation

    if (saveModelHistory_) { save(model_    , "model_0"  ); }
    DOSAVE save(response_ , "response_0");
    DOSAVE save(modelRef_ , "modelRef_0");
    DOSAVE save(RVector(response_ / data_ -1.0), "deltaData_0");
    DOSAVE forward_->constraints()->save("constraint.matrix");
    DOSAVE save(constraintWeights_, "cweight_0");
    DOSAVE save(modelWeight_, "mweight_0");
    DOSAVE save(*forward_->jacobian(), "sens.bmat");

    DOSAVE std::cout << "C size: " << forward_->constraints()->cols()
                     << " x " << forward_->constraints()->rows() << std::endl;

    double phiD = getPhiD();

    if (verbose_) {
        echoMinMax(data_,  "data");
        echoMinMax(error_,  "error");
        echoMinMax(response_,  "response");
        if (haveReferenceModel_) {
            echoMinMax(modelRef_,  "reference model");
        } else {
            std::cout << "calc without reference model" << std::endl;
        }

        std::cout << 0 << ": rms/rrms(data, response) = " << rms(data_, response_)
                  << "/" << rrms(data_, response_) * 100.0 << "%" << std::endl;
        std::cout << 0 << ": chi^2(data, response, error, log) = "
        << phiD / data_.size() << std::endl;
        std::cout << 0 << ": Phi = " << phiD << " + " << getPhiM() << " * "
                                    << lambda_ << " = " << getPhi() << std::endl;
    }

    //** Start iteration
    iter_ = 0;
    double oldphi = phiD;

    //** store initial model
    modelHist_.push_back(model_);

    isRunning_ = true;

    while (iter_ < maxiter_ && !abort_){
        if (verbose_) std::cout << "Iter: " << iter_ << std::endl;

        if (!oneStep()) break;
        //** no idea why this should be saved
        //DOSAVE save(*forward_->jacobian() * model_, "dataJac_"  + str(iter_) PLUS_TMP_VECSUFFIX);
        DOSAVE save(response_,    "response_" + str(iter_) PLUS_TMP_VECSUFFIX);

        modelHist_.push_back(model_);

        double phiD = getPhiD();

        if (stopAtChi1_ && (phiD < data_.size())) {
            if (verbose_) std::cout << "Reached data fit criterion (chi^2 <= 1). Stop." << std::endl;
            break;
        }
        double phi = getPhi();
        if (phi / oldphi > (1.0 - dPhiAbortPercent_ / 100.0) && iter_ > 2) {
            if (verbose_) std::cout << "Reached data fit criterion (delta phi < "
                                    << dPhiAbortPercent_ << "%). Stop." << std::endl;
            break;
        }

        oldphi = phi;

        if (isRobust_) robustWeighting();
        if (isBlocky_) constrainBlocky();
        if (lambdaFactor_ > 0.0) max(lambdaMin_, lambda_ *= lambdaFactor_);

    } //** while iteration;
    isRunning_ = false;
    return model_;
} //** run

bool RInversion::oneStep() {
    iter_++;
    deltaModelIter_.resize(model_.size());
    deltaModelIter_ *= 0.0;
    deltaDataIter_ = (tD_->trans(data_) - tD_->trans(response_));

    if (sum(abs(deltaDataIter_)) < TOLERANCE) {
        if (verbose_) std::cout << "sum(abs(deltaDataIter_)) == Zero" << std::endl;
        return false;
    }
    //    Vec deltaModel0(model_.size());
    Vec modelNew(   model_.size());
    Vec responseNew( data_.size());
    Vec roughness(this->constraintsCount(), 0.0);

    this->checkJacobian((recalcJacobian_ && iter_ > 1) || jacobiNeedRecalc_);

    // if ((recalcJacobian_ && iter_ > 1) || jacobiNeedRecalc_ ) {
        //     Stopwatch swatch(true);
        //     if (verbose_) std::cout << "Recalculating Jacobian matrix ...";
        //     forward_->createJacobian(model_);
        //     if (verbose_) std::cout << swatch.duration(true) << " s" << std::endl;
        // }

    if (!localRegularization_) {
        DOSAVE echoMinMax(model_, "model: ");
        roughness = this->roughness();
    } else {
        if (verbose_) std::cout << "use local regularization" << std::endl;
    }

    if (iter_ == 1 && optimizeLambda_) { //optimize regularization strength using L-curve
//        deltaModelIter_ = optLambda(deltaDataIter_, deltaModel0); //!!! h-variante

        /////////////*********************
        // fix this!!!!!!!!!!!!!1 constraintsH != deltaModel0
        /////////////*********************
        deltaModelIter_ = optLambda(deltaDataIter_, constraintsH_);
    } else {
        DOSAVE save(deltaDataIter_, "dd_" + str(iter_) PLUS_TMP_VECSUFFIX);
        DOSAVE echoMinMax(data_,      "data");
        DOSAVE echoMinMax(dataWeight_,  "dW");
        DOSAVE echoMinMax(deltaDataIter_,  "dd");
        DOSAVE echoMinMax(deltaModelIter_, "dm");
//         save(constraintWeights_, "cw.tmp");
        DOSAVE echoMinMax(constraintWeights_,  "cW");
        DOSAVE echoMinMax(modelWeight_,  "mW");
        DOSAVE echoMinMax(model_,    "mod");
        DOSAVE echoMinMax(response_, "resp");
        DOSAVE echoMinMax(tD_->deriv(response_), "dtD");
        DOSAVE echoMinMax(tM_->deriv(model_), "dtM");
//        DOSAVE echoMinMax(deltaModel0, "dM0");
        DOSAVE echoMinMax(constraintsH_, "constraintsH");
        DOSAVE save(constraintsH_, "constraintsH");
        DOSAVE save(tM_->deriv(model_), "modelTrans");
        DOSAVE save(tD_->deriv(response_), "responseTrans");

        if (verbose_) std::cout << "solve CGLSCDWWtrans with lambda = " << lambda_ << std::endl;
        try{
            solveCGLSCDWWhtrans(*forward_->jacobian(), *forward_->constraints(),
                                dataWeight_,
                                deltaDataIter_,
                                deltaModelIter_,
                                constraintWeights_, modelWeight_,
                                tM_->deriv(model_),
                                tD_->deriv(response_),
                                lambda_, roughness, maxCGLSIter_, CGLStol_,
                                dosave_);
        } catch(...){
            // __MS("Debug halt! 1700")
            // forward_->mesh()->save("S1700.bms");
            // forward_->regionManager().mesh().save("R1700.bms");

            // __MS("Debug halt! 17708")
            // forward_->mesh()->save("S17708.bms");
            // forward_->regionManager().mesh().save("R17708.bms");

            // std::cout<<forward_->mesh()->boundary(2121).marker() << std::endl;
            // std::cout<<forward_->mesh()->boundary(2122).marker() << std::endl;

            // exit(1);
            log(Error, "solveCGLSCDWWhtrans throws unknown exception.");
        }
    } // else no optimization
    //restrictMax(deltaModelIter_, 50.0); // only work for log trans
    DOSAVE echoMinMax(deltaModelIter_, "dm");
    modelNew = tM_->update(model_, deltaModelIter_);

    if (haveInfNaN(modelNew)){
        save(deltaModelIter_, "dmodel_Nan" PLUS_TMP_VECSUFFIX);
        save(model_, "mmodel_Nan" PLUS_TMP_VECSUFFIX);
        save(modelNew, "newmodel_Nan" PLUS_TMP_VECSUFFIX);
        log(Error, "Model contains nan values.");
    }

    DOSAVE save(model_, "oldmodel");
    DOSAVE save(deltaModelIter_, "deltaModel");

    if (dosave_) {
        save(modelNew, "model_" + str(iter_) PLUS_TMP_VECSUFFIX);
    } else {
        if (saveModelHistory_) save(modelNew, "modelLS");
    }

    Vec responseLast(response_);
    responseNew = forward_->response(modelNew);

    double tau = 1.0;
    if (useLinesearch_){
        tau = linesearch(modelNew, responseNew);
    }

    if (tau >= 0.95){ //! full step possible;
        response_ = responseNew;
    } else { //! normal line search parameter between 0.03 and 0.94
        modelNew = tM_->update(model_, deltaModelIter_ * tau);
        response_ = forward_->response(modelNew);
    }

    model_ = modelNew;
    if (saveModelHistory_) save(model_, "model_" + str(iter_) PLUS_TMP_VECSUFFIX);

    if (verbose_) echoStatus();

    if (doBroydenUpdate_) { //** perform Broyden update;
    // did not yet reflect moving jacobian into modellingbase
    THROW_TO_IMPL
//         if (verbose_) std::cout << "perform Broyden update" << std::endl;
//         Vec u(tD_->trans(response_) - tD_->trans(responseLast) - (*J_ * deltaModelIter_));
//         Vec v(deltaModelIter_ / dot(deltaModelIter_, deltaModelIter_));
//         rank1Update(*J_, u, v);
    }

    //!** temporary stuff
    if (forward_->mesh()){
// this forces pygimli/generatecode.py to create a ugly log10 declaration, which overwrites the valid log10 declarion
        // DOSAVE forward_->mesh()->addData("F-op-model(log10)",
        //                                         log10(forward_->mesh()->cellAttributes()));
        DOSAVE forward_->mesh()->addData("F-op-model",
                                                forward_->mesh()->cellAttributes());
        DOSAVE forward_->mesh()->exportVTK("fop-model" + str(iter_));
        DOSAVE forward_->mesh()->clearData();
    }
    return true;
}

RVector RInversion::optLambda(const RVector & deltaData,
                              const RVector & deltaModel0) {

ALLOW_PYTHON_THREADS

    std::vector< double > phiM, phiD;
    double ys = 0.0, yss = 0.0, curv = 0.0, oldcurv = -1e10;
    Vec oldDModel(model_.size());
    Vec uroldDModel(oldDModel);
    Vec deltaModel(model_.size());
    Vec tModel(tM_->trans(model_));
    Vec tResponse(tM_->trans(response_));
    Vec roughness(constraintWeights_.size(), 0.0);

    if (!localRegularization_) {
        DOSAVE echoMinMax(model_, "model: ");
        roughness = this->roughness();
    } else {
        if (verbose_) std::cout << "use local regularization" << std::endl;
    }

    DOSAVE echoMinMax(modelWeight_,  "mW");
    DOSAVE echoMinMax(constraintsH_, "constraintsH");
    DOSAVE save(constraintsH_, "constraintsH");

//    solveCGLSCDWWtrans(*J_, forward_->constraints(), dataWeight_, deltaData, deltaModel, constraintWeights_,
//                        modelWeight_, tM_->deriv(model_), tD_->deriv(response_),
//                        lambda_, deltaModel0, maxCGLSIter_, verbose_);
    solveCGLSCDWWhtrans(*forward_->jacobian(), *forward_->constraints(),
                        dataWeight_, deltaDataIter_, deltaModel,
                        constraintWeights_, modelWeight_,
                        tM_->deriv(model_), tD_->deriv(response_),
                        lambda_, roughness, maxCGLSIter_, dosave_);

    Vec appModelStart(tM_->invTrans(tModel + deltaModel));
    DOSAVE save(appModelStart, "appModel");
    double phiMNorm = (getPhiM(appModelStart));
    double phiDNorm = (getPhiD());

    //* normalization
//    phiM.push_back(getPhiM() / phiMNorm);
//    phiD.push_back(1.0);
    //* override normalization
    phiM.push_back(std::log(getPhiM(appModelStart))); phiMNorm = 1.0;
    phiD.push_back(std::log(getPhiD())); phiDNorm = 1.0;
    if(verbose_) std::cout << "lambda(0) = inf" << " PhiD = " << phiD.back() << " PhiM = " << phiM.back()  << std::endl;

    int lambdaIter = 0;
    while (lambdaIter < 30) {
        lambdaIter++;
        if(verbose_) std::cout << lambdaIter << "lambda = " << lambda_ << std::endl;
//        solveCGLSCDWWtrans(*J_, forward_->constraints(), dataWeight_, deltaData, deltaModel, constraintWeights_,
//                          modelWeight_, tM_->deriv(model_), tD_->deriv(response_),
//                          lambda_, deltaModel0, maxCGLSIter_, verbose_);
        solveCGLSCDWWhtrans(*forward_->jacobian(), *forward_->constraints(),
                            dataWeight_, deltaDataIter_, deltaModel,
                            constraintWeights_, modelWeight_,
                            tM_->deriv(model_), tD_->deriv(response_),
                            lambda_, roughness, maxCGLSIter_, dosave_);

        Vec appModel(tM_->invTrans(tModel + deltaModel));
        Vec appResponse(tD_->invTrans(tResponse + *forward_->jacobian() * deltaModel));

        if (lambdaIter == 1) { //* normalize on 1st iteration
//      	 phiMNorm = getPhiM(appModel);
//      	 phiDNorm = getPhiD(appResponse);
        }

        phiM.push_back(std::log(getPhiM(appModel)) / phiMNorm);
        phiD.push_back(std::log(getPhiD(appResponse)) / phiDNorm);

        if (lambdaIter > 1) {
            ys = (phiD[ lambdaIter ] - phiD[ lambdaIter - 2 ]) / (phiM[ lambdaIter ] - phiM[ lambdaIter - 2 ]);
            yss = ((phiD[ lambdaIter ] - phiD[ lambdaIter - 1 ]) / (phiM[ lambdaIter ] - phiM[ lambdaIter - 1]) -
              (phiD[ lambdaIter - 1 ] - phiD[ lambdaIter - 2 ]) / (phiM[ lambdaIter - 1 ] - phiM[ lambdaIter - 2 ])) /
              (phiM[ lambdaIter ] - phiM[ lambdaIter - 2 ]) * 2.0;
            curv = yss / std::pow(1 + ys * ys, 1.5);

            if(verbose_) std::cout << " lambda(" << lambdaIter << ") = " << lambda_ << " PhiD = "
                                     << phiD.back() << " PhiM = " << phiM.back() << " curv = " << curv << std::endl;
            if ((curv < oldcurv) && (lambdaIter > 4)) {
                deltaModel = uroldDModel;
                lambda_ /= (0.8 * 0.8);
                if (verbose_) std::cout << lambdaIter << ": lambdaIter -- " << "Curvature decreasing, choosing lambda = " << lambda_ << std::endl;
                break;
            }
        oldcurv = curv;
        } else { //** lambdaIter > 1
            if(verbose_) std::cout << "lambda(" << lambdaIter << ") = " << lambda_
                                     << " PhiD = " << phiD.back() << " PhiM = " << phiM.back()  << std::endl;
        }
        uroldDModel = oldDModel;
        oldDModel = deltaModel;
        lambda_ *= 0.8;
    }  //** while loop;
    DOSAVE save(phiM, "phiM");
    DOSAVE save(phiD, "phiD");

    return deltaModel;
}

} // namespace GIMLI

