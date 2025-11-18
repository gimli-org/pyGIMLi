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

#include "triangleWrapper.h"

#include "mesh.h"
#include "meshentities.h"
#include "node.h"

#if TRIANGLE_FOUND
    extern "C"{
        #define REAL double
        #define VOID void
        #define TRILIBRARY
        #define ANSI_DECLARATORS

        #include <triangle.h>
        #define USE_LIBTRIANGLE TRUE
    }
#elif HAVE_LIBTRIANGLE
    #define USE_LIBTRIANGLE TRUE

extern "C"{

#define REAL double
#define VOID void
#define TRILIBRARY
#define ANSI_DECLARATORS

struct triangulateio {
  REAL *pointlist;                                               /* In / out */
  REAL *pointattributelist;                                      /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  int numberofpoints;                                            /* In / out */
  int numberofpointattributes;                                   /* In / out */

  int *trianglelist;                                             /* In / out */
  REAL *triangleattributelist;                                   /* In / out */
  REAL *trianglearealist;                                         /* In only */
  int *neighborlist;                                             /* Out only */
  int numberoftriangles;                                         /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftriangleattributes;                                /* In / out */

  int *segmentlist;                                              /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  int numberofsegments;                                          /* In / out */

  REAL *holelist;                        /* In / pointer to array copied out */
  int numberofholes;                                      /* In / copied out */

  REAL *regionlist;                      /* In / pointer to array copied out */
  int numberofregions;                                    /* In / copied out */

  int *edgelist;                                                 /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  REAL *normlist;                /* Used only with Voronoi diagram; out only */
  int numberofedges;                                             /* Out only */
};

void triangulate(char *, struct triangulateio *, struct triangulateio *, struct triangulateio *);

}
#else // HAVE_LIBTRIANGLE
    #define USE_LIBTRIANGLE 0
#endif

namespace GIMLI{

TriangleWrapper::TriangleWrapper(const Mesh & inMesh) : inMesh_(& inMesh){
    init_();
}

TriangleWrapper::TriangleWrapper(const Mesh & inMesh, Mesh & outMesh,
                                 const std::string & triSwitches)
    : inMesh_(& inMesh) {
    init_();
    switches_ = triSwitches;
    generate(outMesh);
}

TriangleWrapper::~TriangleWrapper(){
    freeMemory_();
#if USE_LIBTRIANGLE
    if (mesh_input_)          delete mesh_input_;
    if (mesh_output_)         delete mesh_output_;
    if (mesh_voronoi_output_) delete mesh_voronoi_output_;
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

void TriangleWrapper::init_(){
    switches_ = "-pze";
#if USE_LIBTRIANGLE
    mesh_input_           = new struct triangulateio;
    mesh_output_          = new struct triangulateio;
    mesh_voronoi_output_  = new struct triangulateio;
    allocateOutMemory_();
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

void TriangleWrapper::setSwitches(const std::string & s){
    switches_ = s;
}

void TriangleWrapper::generate(Mesh & mesh){
#if USE_LIBTRIANGLE

    if (inMesh_->nodeCount() < 3){
        throwError(WHERE_AM_I + " input mesh must have at least 3 nodes ");
    }
    if (mesh_output_->pointlist) {
        freeMemory_();
        allocateOutMemory_();
    }
    transformMeshToTriangle_(*inMesh_, * mesh_input_);
    triangulate((char *)switches_.c_str(), mesh_input_, mesh_output_, mesh_voronoi_output_);
    transformTriangleToMesh_(* mesh_output_, mesh);
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

Mesh TriangleWrapper::generate(){
    Mesh mesh(2);
    generate(mesh);
    return mesh;
}

void TriangleWrapper::allocateOutMemory_(){
#if USE_LIBTRIANGLE

    mesh_output_->pointlist             = (REAL *) NULL;
    mesh_output_->pointattributelist    = (REAL *) NULL;
    mesh_output_->pointmarkerlist       = (int *) NULL;

    mesh_output_->trianglelist          = (int *) NULL;
    mesh_output_->triangleattributelist = (REAL *) NULL;

    mesh_output_->segmentlist           = (int *) NULL;
    mesh_output_->segmentmarkerlist     = (int *) NULL;

    mesh_output_->edgelist              = (int *) NULL;
    mesh_output_->edgemarkerlist        = (int *) NULL;
    mesh_output_->normlist              = (REAL *) NULL;
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

void TriangleWrapper::freeMemory_(){
#if USE_LIBTRIANGLE
    if (mesh_output_->pointlist != (REAL*)NULL)           free(mesh_output_->pointlist);
    if (mesh_output_->pointattributelist != (REAL*)NULL)  free(mesh_output_->pointattributelist);
    if (mesh_output_->pointmarkerlist != (int*)NULL )     free(mesh_output_->pointmarkerlist) ;

    if (mesh_output_->trianglelist != (int*)NULL)         free(mesh_output_->trianglelist);
    if (mesh_output_->triangleattributelist != (REAL*)NULL) free(mesh_output_->triangleattributelist);

    if (mesh_output_->segmentlist != (int*)NULL)          free(mesh_output_->segmentlist);
    if (mesh_output_->segmentmarkerlist != (int*)NULL)    free(mesh_output_->segmentmarkerlist);

    if (mesh_output_->edgelist != (int*)NULL)             free(mesh_output_->edgelist);
    if (mesh_output_->edgemarkerlist != (int*)NULL)       free(mesh_output_->edgemarkerlist);
    if (mesh_output_->normlist != (REAL *) NULL)          free(mesh_output_->normlist);

    delete [] mesh_input_->pointlist;
    delete [] mesh_input_->pointmarkerlist;
    delete [] mesh_input_->segmentlist;
    delete [] mesh_input_->segmentmarkerlist;
    delete [] mesh_input_->holelist;
    delete [] mesh_input_->regionlist;
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

void TriangleWrapper::transformTriangleToMesh_(const triangulateio & trimesh, Mesh & mesh){
#if USE_LIBTRIANGLE
    mesh.clear();

    //!  Nodes;
    int marker = 0;
    for (int i = 0; i < trimesh.numberofpoints; i ++){
        if (trimesh.pointmarkerlist != (int*)NULL) marker = trimesh.pointmarkerlist[i];
        mesh.createNode(trimesh.pointlist[i * 2], trimesh.pointlist[i * 2 + 1], 0.0, marker);
    }

    //!  Edges / Segments
    if (trimesh.numberofedges == 0){
        for (int i = 0; i < trimesh.numberofsegments; i ++){
            if (trimesh.segmentmarkerlist != (int*)NULL) marker = trimesh.segmentmarkerlist[i];
            mesh.createEdge(mesh.node(trimesh.segmentlist[i * 2]), mesh.node(trimesh.segmentlist[i * 2 + 1]), marker);
        }
    } else {
        if (trimesh.edgelist != NULL){
            for (int i = 0; i < trimesh.numberofedges; i ++){
                if (trimesh.edgemarkerlist != (int*)NULL) marker = trimesh.edgemarkerlist[i];
                mesh.createEdge(mesh.node(trimesh.edgelist[i * 2]), mesh.node(trimesh.edgelist[i * 2 + 1]), marker);
            }
        } else {
            std::cout << "Warning! " << WHERE_AM_I << " edges are not exported. Append -e flag to the triangle command." << std::endl;
        }
    }

    //! Triangles;
    int a = 0, b = 0, c = 0;
    double attribute = 0.0;
    for (int i = 0; i < trimesh.numberoftriangles; i ++){
        if (trimesh.triangleattributelist != (REAL*)NULL) attribute = trimesh.triangleattributelist[i];

        a = trimesh.trianglelist[i * 3];
        b = trimesh.trianglelist[i * 3 + 1];
        c = trimesh.trianglelist[i * 3 + 2];

        mesh.createTriangle(mesh.node(a), mesh.node(b), mesh.node(c), (size_t)attribute);
    }

// //** RegionMarker
//   for (int i = 0; i < trimesh.numberofregions * 4; i += 4){
//     mesh.createRegionMarker(trimesh.regionlist[i], trimesh.regionlist[i + 1],
// 			      trimesh.regionlist[i + 2], trimesh.regionlist[i + 3]);
//   }
//   for (int i = 0; i < trimesh.numberofholes * 2; i += 2){
//     mesh.createRegionMarker(trimesh.holelist[i], trimesh.holelist[i + 1], -1, -1);
//   }
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

void TriangleWrapper::transformMeshToTriangle_(const Mesh & mesh,
                                               triangulateio & trimesh){
#if USE_LIBTRIANGLE
    //! node section
    Index nVerts = mesh.nodeCount();

    trimesh.numberofpoints = nVerts;
    trimesh.numberofpointattributes = 0;  // only one Parameter

    trimesh.pointlist = new double[2 * nVerts];
    trimesh.pointmarkerlist = new int[nVerts];

    for (Index i = 0; i < nVerts; i ++){
//         __MS(mesh.node(i).x() << " " << mesh.node(i).y())
        trimesh.pointlist[i * 2]      = mesh.node(i).x();
        trimesh.pointlist[i * 2 + 1]  = mesh.node(i).y();
        trimesh.pointmarkerlist[i]    = mesh.node(i).marker();
    }

    //! edge section;
    Index nEdges      = mesh.boundaryCount();
    Index edgeCounter = 0;

    trimesh.numberofsegments    = nEdges;
    trimesh.segmentlist         = new int[2 * nEdges];
    trimesh.segmentmarkerlist   = new int[nEdges];

    for (Index i = 0; i < nEdges; i ++){
        Boundary * edge = &mesh.boundary(i);
        trimesh.segmentlist[edgeCounter * 2]      = edge->node(0).id();
        trimesh.segmentlist[edgeCounter * 2 + 1]  = edge->node(1).id();
        trimesh.segmentmarkerlist[edgeCounter]    = edge->marker();
        edgeCounter++;
    }

    //! Holes;
    Index nHoles = mesh.holeMarker().size();//domain.holeCount();
    trimesh.numberofholes = nHoles;
    trimesh.holelist = new double[2 * nHoles + 1];

    for (Index i = 0; i < nHoles; i ++){
        trimesh.holelist[i * 2] = mesh.holeMarker()[i].x();
        trimesh.holelist[i * 2 + 1] = mesh.holeMarker()[i].y();
     }

    if (mesh.cellCount() > 0){
        //! Cell Regions;
        Index nRegions = mesh.cellCount();//domain.regionCount();
        trimesh.numberofregions = nRegions;
        trimesh.regionlist = new double[4 * nRegions + 1];

        Index count = 0;
        for (Index i = 0; i < mesh.cellCount(); i ++ ){
            trimesh.regionlist[count * 4] = mesh.cell(i).center().x();
            trimesh.regionlist[count * 4 + 1] = mesh.cell(i).center().y();
            trimesh.regionlist[count * 4 + 2] = mesh.cell(i).marker();
            trimesh.regionlist[count * 4 + 3] = 0.0;
            count ++;
        }
    } else {
        //! Regions;
        Index nRegions = mesh.regionMarkers().size();//domain.regionCount();
        trimesh.numberofregions = nRegions;
        trimesh.regionlist = new double[4 * nRegions + 1];

        Index count = 0;
        for (auto & r : mesh.regionMarkers()){
    //         std::cout << it->pos().x() << " " << it->pos().y() << " "
    //                   << it->marker() << " "
    //                   << it->area() << std::endl;
            trimesh.regionlist[count * 4] = r.x();
            trimesh.regionlist[count * 4 + 1] = r.y();
            trimesh.regionlist[count * 4 + 2] = r.marker();
            trimesh.regionlist[count * 4 + 3] = r.area();
            count ++;
        }
    }
#else
    std::cerr << WHERE_AM_I << " Triangle not installed" << std::endl;
#endif
}

} // namespace GIMLI
