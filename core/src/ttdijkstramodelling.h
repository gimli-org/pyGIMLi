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

#ifndef _GIMLI_TTDIJKSTRAMODDELING__H
#define _GIMLI_TTDIJKSTRAMODDELING__H

#include "gimli.h"
#include "modellingbase.h"
#include "mesh.h"

namespace GIMLI {

class GraphDistInfo{
public:    
    GraphDistInfo()
        :  time_(0.0), dist_(0.0){
    }
    GraphDistInfo(double t, double d)
        :  time_(t), dist_(d){
    }
    GraphDistInfo(double t, double d, Index cellID)
        :  time_(t), dist_(d){
        cells_.insert(cellID);
    }
    
    /*! Travel time for the current way element. Depends on the graph initialization.*/
    void setTime(double t) { time_ = t; }

    /*! Travel time for the current way element. Depends on the graph initialization.*/
    double time() const { return time_; }

    /*! Distance of the way element.*/
    double dist() const { return dist_; }

    /*! Index for all cells containing this wayelement*/
    std::set < Index > & cellIDs() { return cells_; }

    /*! Index for all cells containing this wayelement*/
    const std::set < Index > & cellIDs() const { return cells_; }

protected:

    double time_;
    double dist_;
    std::set < Index > cells_;
};

//** sorted vector
typedef std::map< Index, GraphDistInfo > NodeDistMap;
//** sorted matrix
typedef std::map< Index, NodeDistMap > Graph;


/*! Dijkstra's shortest path finding*/
class DLLEXPORT Dijkstra {
public:
    Dijkstra();

    Dijkstra(const Graph & graph);

    ~Dijkstra(){}

    void setGraph(const Graph & graph);

    void setStartNode(Index startNode);

    /*!Get the shortest way from root to node. Inline version.*/
    void shortestPathTo(Index node, IndexArray & way) const;

    /*!Get the shortest way from root to node.*/
    IndexArray shortestPathTo(Index node) const;

    /*!Get the shortest way from node index start to end.*/
    IndexArray shortestPath(Index start, Index end);

    /*!Distance from root to node.*/
    double distance(Index root, Index node);
    
    /*!Distance to node to the last known root.*/
    double distance(Index node); 

    /*!All distances to root.*/
    RVector distances(Index root);
    
    /*!All distances from to last known root.*/
    RVector distances() const;
    
    Graph & graph() {
        return graph_;
    }

    const Graph & graph() const {
        return graph_;
    }

    GraphDistInfo graphInfo(Index na, Index nb) {
        return graph_[na][nb];
    }

    class Edge_ : std::pair< Index, Index > {
    public:
        Edge_() : start(0), end(0) {}
        Edge_(Index a, Index b) : start(a), end(b) {}
        Index start;
        Index end;
    };

    /*! Definition for the priority queue */
    class DistancePair_ : std::pair< double, Edge_ > {
    public:
        DistancePair_() : first(0.0) {}
        DistancePair_(double f, Edge_ & s) : first(f), second(s) {}

        double first;
        Edge_ second;
    };

    template < class T > class ComparePairsClass_ {
    public:
        bool operator() (const T & lhs, const T & rhs) {
            return lhs.first > rhs.first;
        }
    };

protected:
    std::vector < Edge_ > pathMatrix_;
    
    NodeDistMap distances_;

    Graph graph_;
    Index _root;
};

//! Modelling class for travel time problems using the Dijkstra algorithm
/*! TravelTimeDijkstraModelling(mesh, datacontainer) */
class DLLEXPORT TravelTimeDijkstraModelling : public ModellingBase {
public:
    TravelTimeDijkstraModelling(bool verbose=false);

    TravelTimeDijkstraModelling(Mesh & mesh,
                                DataContainer & dataContainer,
                                bool verbose=false);

    virtual ~TravelTimeDijkstraModelling() { }

    virtual RVector createDefaultStartModel();

    RVector createGradientModel(double lBound, double uBound);

    /*! Interface. Calculate response */
    virtual RVector response(const RVector & slowness);

    /*! Interface. */
    virtual void createJacobian(const RVector & slowness);

    /*! Interface. */
    virtual void initJacobian();

    Graph createGraph(const RVector & slownessPerCell) const;

//     RVector calculate();

    double findMedianSlowness() const;

    RVector getApparentSlowness() const;

    void createJacobian(RSparseMapMatrix & jacobian, const RVector & slowness);

    /*! Returns the mesh node indieces for the way from shot to receiver, 
    respective the data for the last Jacobian calculation. 
    If you want further infos about the way element. 
    You can ask the Dijkstra about the graph infos. */
    const IndexArray & way(Index sht, Index rec) const;

    /*! Read only access to the recent dijktra. */
    const Dijkstra & dijkstra() const { return dijkstra_; };

protected:

    /*! Automatically looking for shot and receiver points if the mesh is changed. */
    virtual void updateMeshDependency_();

    Dijkstra dijkstra_;
    double background_;

    /*! Nearest nodes for the current mesh for all shot points.*/
    IndexArray shotNodeId_;

    /*! Map shot id to sequential shot node number of shotNodeId_ */
    std::map< Index, Index > shotsInv_;

    /*! Nearest nodes for the current mesh for all receiver points.*/
    IndexArray receNodeId_;

    /*! Map receiver id to sequential receiver node number of receNodeId_ */
    std::map< Index, Index > receiInv_;

    /*! Way matrix of the last full jacobian generation. */
    std::vector < std::vector < IndexArray > > wayMatrix_;

};

/*! New Class derived from standard travel time modelling */
class DLLEXPORT TTModellingWithOffset: public TravelTimeDijkstraModelling{
public:
    TTModellingWithOffset(Mesh & mesh, DataContainer & dataContainer, bool verbose);

    virtual ~TTModellingWithOffset();

    virtual RVector createDefaultStartModel();

    virtual RVector response(const RVector & model);

    void initJacobian();

    virtual void createJacobian(const RVector & slowness);

    size_t nShots(){ return shots_.size(); }

protected:
    RVector                  shots_;
    std::map< Index, Index > shotMap_;
    Mesh                     offsetMesh_;
};


} //namespace GIMLI

#endif
