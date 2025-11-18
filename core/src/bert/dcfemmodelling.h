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

#ifndef _BERT_DCFEMMODDELING__H
#define _BERT_DCFEMMODDELING__H

#include "bert.h"
#include "datamap.h"

#include <modellingbase.h>
#include <sparsematrix.h>
#include <pos.h>
#include <matrix.h>

#include <vector>

namespace GIMLI{

/*! if fix is set. Matrix will check and fix singularities. Do not fix the matrix if you need it for the rhs while singulariety removal calculation. */
DLLEXPORT void dcfemDomainAssembleStiffnessMatrix(RSparseMatrix & S, const Mesh & mesh,
                                                  double k=0.0, bool fix=true);

DLLEXPORT void dcfemDomainAssembleStiffnessMatrix(CSparseMatrix & S, const Mesh & mesh,
                                                  double k=0.0, bool fix=true);

// DLLEXPORT void assembleStiffnessMatrixHomogenDirichletBC(RSparseMatrix & S,
//                                                          const IndexArray & nodeID);

DLLEXPORT double mixedBoundaryCondition(const Boundary & boundary,
                                        const RVector3 & sourcePos,
                                        double k=0.0);

DLLEXPORT void assembleCompleteElectrodeModel(RSparseMatrix & S,
                                              const std::vector < ElectrodeShape * > & elecs,                                               
                                              uint oldMatSize, bool lastIsReferenz,
                                              const RVector & contactImpedances);

/*!ERT utility function for the handling of complex resistivity
 * values vs. amplitude/phase data.
 * Data are usually given in amplitude(Ohm m) and phase(mRad).
 * Internal we calculate with real and complex resistivity values.
 * There is no complex data type for cell attributes so we put the complex
 * values into the mesh data map as AttributeReal and AttributeImag.
 * Don't do this directly until you know what u do.
 * Use this utility functions to apply the complex resistivity values.
*/
DLLEXPORT void setComplexResistivities(Mesh & mesh, const CVector & res);

/*! Apply a rho map of Amplitude(Ohm m) and Phase(rad) to
 * set desired complex resistivity values.
 * res=am * cos(ph) - i am * sin(abs(ph))
 * See \ref setComplexResistivities(Mesh & mesh, const CVector & res).
 */
DLLEXPORT void setComplexResistivities(Mesh & mesh,
                                      const std::map < float, Complex > & aMap);

/*! Apply the vectors of Amplitude(Ohm m) and Phase(rad) to
 * set desired complex resistivity values res.
 * res=am * cos(ph) - i am * sin(abs(ph))
 * See \ref setComplexResistivities(Mesh & mesh, const CVector & res).
 */
DLLEXPORT void setComplexResistivities(Mesh & mesh,
                                       const RVector & amp,
                                       const RVector & phase);

/*! Return CVector of the complex resistivity values.
 */
DLLEXPORT CVector getComplexResistivities(const Mesh & mesh);

/*! Set the complex resistivity data values into a DataContainer.
 * as data('u') in Volt and data('ip') in mRad.
 */
DLLEXPORT void setComplexData(DataContainer & data, const CVector & u);

/*! Set the complex resistivity data values into a DataContainer.
 * Complex values are transformed into data('u') in Volt and data('ip') in mRad.
 * u=abs(re -i im), ip=phase(re -i im) * 1000.
 */
DLLEXPORT void setComplexData(DataContainer & data,
                              const RVector & re, const RVector & im);

/*! Return CVector of the complex resistivity values transformed from
 * data('rhoa') and data('ip')
 * z=am * cos(ph) - i am * sin(abs(ph)),
 * with am=data('rhoa') and ph=data('ip') /1000
 */
DLLEXPORT CVector getComplexData(const DataContainer & data);


class DLLEXPORT DCMultiElectrodeModelling : public GIMLI::ModellingBase {
public:
    DCMultiElectrodeModelling(bool verbose=false);

    DCMultiElectrodeModelling(Mesh & mesh,
                              bool verbose=false);

    DCMultiElectrodeModelling(DataContainerERT & dataContainer,
                              bool verbose=false);

    DCMultiElectrodeModelling(Mesh & mesh, DataContainerERT & dataContainer,
                              bool verbose=false);

    virtual ~DCMultiElectrodeModelling();

    /*!Assemble bypasses into the system matrix, i.e., defined shorting circuits that my be useful in some cases, e.g., long electrodes.
     * You can provide a file 'bypass.map' in your your working path. \n
     * This can be overwritten by setting a alternative name for the bypass maps. \ref setBypassMapFile. \n
     * Nodes with a marker leq than -10000 will be recognized and combined into a single bypass electrode with bypass resistance 1e-6 Ohm between the first and all other nodes.\n
     * You can provide a file electrodeBypassResistances.map in the working path to specify the bypass resistance in Ohm for all bypass nodes.
     * e.g,
     * -10000 resistance
     * -10001 resistance
     */

    void assembleStiffnessMatrixDCFEMByPass(RSparseMatrix & S);

    void assembleStiffnessMatrixDCFEMByPass(CSparseMatrix & S);

    virtual void createJacobian(const RVector & model);

    virtual void createConstraints();

    DataContainerERT & dataContainer() const ;

    /*! Essential for abstract class ModellingBase,
     * This method will run by \ref startModel() if no
     * startmodel vector is predefined. */
    virtual RVector createDefaultStartModel();

//     RVector operator () (const RVector & model, double background=0) {
//         return response(model, background); }

    /*! Map resistivity model to the cell attributes for further calculation.
     * Default is resistivity=model[cell.marker()].
     * Empty cells will set to background or prolongated for background==-9e99.
     * If model size equals the mesh().cellCount() then the values are just
     * copied. i.e., resistivity=model[cell.id()]. */
    void mapERTModel(const CVector & model, Complex background);

    /*! Map resistivity model to the cell attributes for further calculation.
     * Default is resistivity=model[cell.marker()].
     * Empty cells will set to background or prolongated for background==-9e99.
     * If model size equals the mesh().cellCount() then the values are just
     * copied. i.e., resistivity=model[cell.id()]. */
    void mapERTModel(const RVector & model, double background);

    /*! Calculate response for a given resistivity model.
     * Either cell based or marker based. See \ref mapERTModel */
    virtual RVector response(const RVector & model){return response(model, -9e99);}

    /*! Calculate response for a given resistivity model.
     * Either cell based or marker based. See \ref mapERTModel */
    RVector response(const RVector & model, double background);

    void createCurrentPattern(std::vector < ElectrodeShape * > & eA,
                              std::vector < ElectrodeShape * > & eB,
                              bool reciprocity);

    virtual void calculate(DataMap & dMap);

    virtual void calculate(DataContainerERT & data, bool reciprocity=false);

    virtual void calculate(const std::vector < ElectrodeShape * > & eA,
                           const std::vector < ElectrodeShape * > & eB);

    virtual void calculateK(const std::vector < ElectrodeShape * > & eA,
                            const std::vector < ElectrodeShape * > & eB,
                            RMatrix & solutionK, int kIdx);

    virtual void calculateK(const std::vector < ElectrodeShape * > & eA,
                            const std::vector < ElectrodeShape * > & eB,
                            CMatrix & solutionK, int kIdx);

    template < class ValueType >
    void calculateKAnalyt(const std::vector < ElectrodeShape * > & eA,
                          const std::vector < ElectrodeShape * > & eB,
                          Matrix < ValueType > & solutionK,
                          double k, int kIdx) const;

    void setTopography(bool topography) { topography_=topography; }
    bool topography() const { return topography_; }

    bool neumann() const { return neumannDomain_; }

    /*! Force complex resistivity calculation if it not set in prior to the mesh.*/
    void setComplex(bool c);

    /*! Return true if the valuetype is complex.*/
    bool complex() const { return complex_; }

    /*! Use dipole-current pattern instead of default pol-current pattern.
    The needed current pattern are collected from the data container.
    Note, that is different to calculate with reference electrode since
    this can be also handled as pol-current pattern while data collection.
    If use dipole-current the reciprocity data collection is disabled. */
    void setDipoleCurrentPattern(bool dipole) {dipoleCurrentPattern_=dipole;}

    /*! Return if dipole current pattern is used. */
    bool dipoleCurrentPattern() const { return dipoleCurrentPattern_; }

    /*! Force analytical calculation. */
    inline void setAnalytical(bool ana){ analytical_=ana; }
    inline bool analytical() const { return analytical_; }

    void collectSubPotentials(RMatrix & subSolutions){
        subSolutions_=& subSolutions;
    }

    /*! Return dipole current pattern map.
     *  corresponds to < CurrentPattern,Idx of Potentialmatrix > */
    inline const std::map < Index, Index > & currentPatternIdxMap() const {
        return currentPatternIdxMap_;
    }

    /*! Provide a file with a list of pairs:\n
     * electrodeID electrodeID resistance \n
     * or\n
     * nodeMarker resistance \n
     * that are short circuited. */
    inline void setBypassMapFile(const std::string & fileName){
        byPassFile_=fileName;
    }

    inline void setkValues(const RVector & v) { kValues_=v; }
    inline const RVector & kValues() const { return kValues_; }

    inline void setWeights(const RVector & v) { weights_=v; }
    inline const RVector & weights() const { return weights_; }

    inline const std::vector< ElectrodeShape * > & electrodes() const {
        return electrodes_;
    }

    void setContactImpedances(const RVector & zi);

    virtual RVector calcGeometricFactor(const DataContainerERT & data,
                                        Index nModel=0);

    virtual void preCalculate(const std::vector < ElectrodeShape * > & eA,
                              const std::vector < ElectrodeShape * > & eB){
    }

    /*! Estimate singular value for the calulated potential if the
     * source electrode position lies on an node.*/
    void setSingValue(bool s=true) { setSingValue_=s;}

    /*! Return true if singular value estimation is switched on.*/
    bool isSetSingValue() const { return setSingValue_;}

    /*! Set a custom solver if you don't want the default Choldmod or UMFPACK. */
    void setSolver(SolverWrapper *solver){ solver_ = solver; }
    
private:
    void init_();

protected:

    template < class ValueType >
    void calculateK_(const std::vector < ElectrodeShape * > & eA,
                     const std::vector < ElectrodeShape * > & eB,
                     Matrix < ValueType > & solutionK, int kIdx);

    template < class ValueType >
    void assembleStiffnessMatrixDCFEMByPass_(SparseMatrix < ValueType > & S);

    template < class ValueType >
    DataMap response_(const Vector < ValueType > & model,
                                   ValueType background);

    template < class ValueType >
    Matrix < ValueType > * prepareJacobianT_(const Vector< ValueType > & model);

    RMatrix * prepareJacobian_(const RVector & model);
    CMatrix * prepareJacobian_(const CVector & model);

    void createJacobian_(const RVector & model, const RMatrix & u, RMatrix * J);
    void createJacobian_(const CVector & model, const CMatrix & u, CMatrix * J);

    virtual void deleteMeshDependency_();
    virtual void updateMeshDependency_();
    virtual void updateDataDependency_();

    virtual void searchElectrodes_();

    MatrixBase * subSolutions_;

    bool complex_;

    bool JIsRMatrix_;
    bool JIsCMatrix_;

    bool analytical_;
    bool topography_;
    bool neumannDomain_;
    bool subpotOwner_;
    bool lastIsReferenz_;
    bool setSingValue_;

    std::string byPassFile_;

    RVector kValues_;
    RVector weights_;

    double surfaceZ_;

    IndexArray calibrationSourceIdx_;
    IndexArray bypassNodeIdx_;

    std::vector< ElectrodeShape * > electrodes_;
    ElectrodeShape * electrodeRef_;

    std::vector< ElectrodeShape * > passiveCEM_;

    RVector3 sourceCenterPos_;

    bool buildCompleteElectrodeModel_;

    bool dipoleCurrentPattern_;
    std::map < Index, Index > currentPatternIdxMap_;

    RMatrix potentialsCEM_;
    RVector vContactImpedance_;

    DataMap * primDataMap_;

    SolverWrapper *solver_;
};

class DLLEXPORT DCSRMultiElectrodeModelling : public DCMultiElectrodeModelling {
public:
    DCSRMultiElectrodeModelling(bool verbose=false)
        : DCMultiElectrodeModelling(verbose) {
        init_();
    }

    DCSRMultiElectrodeModelling(Mesh & mesh, bool verbose=false)
        : DCMultiElectrodeModelling(mesh, verbose) {
        init_();
    }

    DCSRMultiElectrodeModelling(Mesh & mesh, DataContainerERT & dataContainer, bool verbose=false)
        : DCMultiElectrodeModelling(mesh, dataContainer, verbose){
        init_();
    }

    virtual ~DCSRMultiElectrodeModelling(){
        if (primPot_  && primPotOwner_)  delete primPot_;
        if (primMesh_ && primMeshOwner_) delete primMesh_;
    };

    virtual void calculateK(const std::vector < ElectrodeShape * > & eA,
                            const std::vector < ElectrodeShape * > & eB,
                            RMatrix & solutionK, int kIdx);

    inline void setPrimaryPotFileBody(const std::string & primPotFileBody){
        primPotFileBody_=primPotFileBody;
    }

    inline void setPrimaryPotential(RMatrix & primPot) { primPot_=& primPot; }

    inline RMatrix & primaryPotential() { return *primPot_; }

    inline void setPrimaryMesh(Mesh & mesh){ primMesh_=&mesh; }

    void setPrimaryMesh(const std::string & meshname);

    inline Mesh & primaryMesh() {return *primMesh_; }

    //const DataMap & primDataMap() const { return ; }

    virtual void preCalculate(const std::vector < ElectrodeShape * > & eA,
                               const std::vector < ElectrodeShape * > & eB);

private:
    void init_(){
        primPotFileBody_= NOT_DEFINED;

        primPotOwner_  =false;
        primPot_       =NULL;

        primMeshOwner_ =false;
        primMesh_      =NULL;
    }

protected:
    virtual void updateMeshDependency_();
    virtual void updateDataDependency_();

    void checkPrimpotentials_(const std::vector < ElectrodeShape * > & eA,
                               const std::vector < ElectrodeShape * > & eB);

    std::string primPotFileBody_;

    bool primPotOwner_;
    RMatrix * primPot_;

    bool primMeshOwner_;
    Mesh * primMesh_;
    Mesh mesh1_;
};

} //namespace BERT

#endif // _BERT_DCFEMMODDELING__H
