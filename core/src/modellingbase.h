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

#ifndef _GIMLI_MODELLINGBASE__H
#define _GIMLI_MODELLINGBASE__H

#include "gimli.h"
#include "matrix.h"
//#include "blockmatrix.h"

namespace GIMLI{

//class H2SparseMapMatrix;

/*! General modelling interface class.*/
class DLLEXPORT ModellingBase{
public:

    ModellingBase(bool verbose=false) ;

    ModellingBase(DataContainer & dataContainer, bool verbose=false);

    ModellingBase(const Mesh & mesh, bool verbose=false);

    ModellingBase(const Mesh & mesh, DataContainer & dataContainer, bool verbose=false);

    virtual ~ModellingBase();

    /*! Set verbose state. */
    void setVerbose(bool verbose);

    /*! Get verbose state. */
    inline bool verbose() const {return verbose_;}

    virtual RVector response(const RVector & model) {
        throwError(WHERE_AM_I + " you need to implement a response "
        "function.");
        return RVector(0);
    }

    /*!Read only response function for multi threading purposes.
     * Index i is use thread counter. */
    virtual RVector response_mt(const RVector & model, Index i=0) const {
        throwError(WHERE_AM_I + " if you want to use read only response "
        "function to use in multi threading environment .. you need to "
        " implement me");
        return RVector(0);
    }

    inline RVector operator() (const RVector & model){ return response(model); }

    /*! Change the associated data container */
    void setData(DataContainer & data);

    /*! Return the associated data container. */
    DataContainer & data() const;

    /*! Interface to define custom start model generators. */
    virtual RVector createDefaultStartModel() { return RVector(0); }

    /*! Return the start model.
     * If there is no model defined createDefaultStartModel is called which can
     * be overloaded by child classes.
     !*/
    virtual RVector startModel();

    virtual void setStartModel(const RVector & startModel);

    /*! Set new mesh to the forward operator, optionally hold region parameter
     * for the new mesh (i.e. for roll a long) */
    void setMesh(const Mesh & mesh, bool ignoreRegionManager=false);

    inline Mesh * mesh() { return mesh_; }

    void createRefinedForwardMesh(bool refine=true, bool pRefine=false);

    /*! Delete the actual mesh. */
    void deleteMesh();

    /*! Set external Jacobian matrix*/
    virtual void setJacobian(MatrixBase * J);

    /*! Create and fill the Jacobian matrix with a given model vector. */
    virtual void createJacobian(const RVector & model);

    /*! Create and fill the Jacobian matrix with a given model vector and
     * a prior calculated response for this model.
     * The Jacobian matrix need to by initialized and resized.
     * If unsure take createJacobian(const RVector & model).*/
    virtual void createJacobian(const RVector & model, const RVector & resp);

    /*! Create and fill the Jacobian matrix with a given model vector.
     * Multi-threaded version.
     * Will be called if setMultiThreadJacobian has been set > 1.
     * For thread safe reasons response_mt, a read only variant of the response
     * method need to be implemented. */
    virtual void createJacobian_mt(const RVector & model, const RVector & resp);

    /*! Here you should initialize your Jacobian matrix. Default is RMatrix()*/
    virtual void initJacobian();

    /*! Return the pointer to the Jacobian matrix associated with this forward operator. */
    MatrixBase * jacobian() { return jacobian_; }

    /*! Return the pointer to the Jacobian matrix associated with this forward operator. */
    MatrixBase * jacobian() const { return jacobian_; }

    /*! Return the Jacobian Matrix (read only) associated with this forward operator.
     *  Throws an exception if the jacobian is not initialized.
     * Cannot yet be overloaded by pyplusplus (return virtual reference)(Warning 1049). */
    virtual RMatrix & jacobianRef() const {
        if (! jacobian_) {
            throwError(WHERE_AM_I + " Jacobian matrix is not initialized.");
        }
        return *dynamic_cast< RMatrix * >(jacobian_);
    }

    virtual RMatrix & jacobianRef() {
        if (! jacobian_) {
            throwError(WHERE_AM_I + " Jacobian matrix is not initialized.");
        }
        return *dynamic_cast< RMatrix * >(jacobian_);
    }

    /*! Clear Jacobian matrix. */
    virtual void clearJacobian(){ jacobian_->clear(); }

    /*! Set and get constraint matrix */
//     inline void setConstraintMatrix(const RSparseMapMatrix & C){ C_ = C; }
//     inline const RSparseMapMatrix & constraintMatrix() const { return C_; }
//         /*! Clear Constraints matrix */
    virtual void setConstraints(MatrixBase * C);

    virtual void clearConstraints();

    virtual void initConstraints();

    virtual void createConstraints();

    virtual MatrixBase * constraints();

    virtual MatrixBase * constraints() const;

    virtual RSparseMapMatrix & constraintsRef() const;

    virtual RSparseMapMatrix & constraintsRef();

    const RMatrix & solution() const { return solutions_; }

    void mapModel(const RVector & model, double background=0);

    /*! Read only extrapolation of model values given per cell marker to
     values given per cell. Exterior values will be set to background or
     prolongated for background != -9e99.
     */
    RVector createMappedModel(const RVector & model, double background=-9e99) const;

    void setRegionManager(RegionManager * reg);

    const RegionManager & regionManager() const;

    RegionManager & regionManager();

    RegionManager & regionManagerRef(){ return regionManager(); }

    bool verbose() { return verbose_; }

    /*! Syntactic sugar for this->regionManager().region(marker). */
    Region * region(int marker);


    RVector createStartModel();

    /*!DEPRECATED use createStartModel*/
    RVector createStartVector();

    void initRegionManager();

    /*! Set the maximum number of allowed threads for MT calculation.
     * Have to be greater than 0.
     * Will also set ENV(OPENBLAS_NUM_THREADS) .. if used.  */
    void setThreadCount(Index nThreads);

    /*! Return the maximum number of allowed threads for MT calculation based on local setting. GIMLI_NUM_THREADS will be used first. Give some verbose output. */
    Index threadCount();

    /*! Set number of threads used for brute force Jacobian generation.
     *1 is default. If nThreads is greater than 1 you need to implement
     * \ref response_mt with a read only response function.
     * Maybe its worth set the single setThreadCount to 1 than,
     * that you don't find yourself in a threading overkill.*/
    void setMultiThreadJacobian(Index nThreads);

    /*! Return number of threads used for Jacobian generation. */
    inline Index multiThreadJacobian() const { return nThreadsJacobian_; }


protected:

    virtual void init_();

    virtual void deleteMeshDependency_(){}

    virtual void updateMeshDependency_(){}

    virtual void updateDataDependency_(){}

    void setMesh_(const Mesh & mesh, bool update=true);

    Mesh                    * mesh_;

    DataContainer           * dataContainer_;

    MatrixBase              * jacobian_;
    bool                    ownJacobian_;

    MatrixBase              * constraints_;
    bool                    ownConstraints_;

    RMatrix                 solutions_;

    RVector                 startModel_;

    bool                    verbose_;

    bool                    regionManagerInUse_;
    bool                    ownRegionManager_;

    Index                   nThreads_;
    Index                   nThreadsJacobian_;

private:
    RegionManager            * regionManager_;

};

/*! Linear modelling class (multiplication by a matrix).
    Calculate response from derivative matrix (Jacobian,Sensitivity) * model
*/
class DLLEXPORT LinearModelling : public ModellingBase {
public:
    LinearModelling(Mesh & mesh, MatrixBase & A, bool verbose=false)
        : ModellingBase(mesh, verbose){
            setJacobian(&A);
        }

    LinearModelling(Mesh & mesh, MatrixBase * A, bool verbose=false)
        : ModellingBase(mesh, verbose){
            setJacobian(A);
        }

    LinearModelling(MatrixBase & A, bool verbose=false);

    virtual ~LinearModelling() { }

    /*! Calculate */
    RVector response(const RVector & model);

    virtual void createJacobian(const RVector & model);

    RVector createDefaultStartModel();

protected:

};



} // namespace GIMLI

#endif // MODELLINGBASE__H
