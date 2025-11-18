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

#include "bertMisc.h"

#include "bertDataContainer.h"
#include "electrode.h"

#include <mesh.h>
#include <node.h>
#include <numericbase.h>

namespace GIMLI{

void initKWaveList(const Mesh & mesh, RVector & kValues, RVector & weights,
                   bool verbose){
    R3Vector sources;
    initKWaveList(mesh, kValues, weights, sources, verbose);
}

void initKWaveList(const Mesh & mesh, RVector & kValues, RVector & weights,
                   const R3Vector & sources, bool verbose){
    kValues.clear();
    weights.clear();

    if (mesh.dim() == 3){
        kValues.resize(1, 0.0);
        weights.resize(1, 1.0);
        return ;
    }

    R3Vector sourcesPos;

    if (sources.empty()){
        sourcesPos = mesh.positions(mesh.findNodesIdxByMarker(MARKER_NODE_ELECTRODE));
        if (verbose) std::cout << "No sensors given for initializing the k wave list .. "
                                  " mesh contains " << sourcesPos.size() << " nodes with marker=-99" << std::endl;
    } else sourcesPos = sources;

    int nElecs = sourcesPos.size();

    double rMin = MAX_DOUBLE, rMax = -MAX_DOUBLE, dist = 0;

    if (nElecs < 2) {
        std::cerr << WHERE_AM_I << " Warning! No sources found, taking some defaults" << std::endl;
        verbose = true;
        rMin = 1.0;
        rMax = (mesh.xMax() - mesh.xMin()) / 20.0;
    } else {
        for (int i = 0; i < nElecs; i ++){
            for (int j = i + 1; j < nElecs; j ++){
                dist = sourcesPos[i].dist(sourcesPos[j]);
                rMin = std::min(dist, rMin);
                rMax = std::max(dist, rMax);
            }
        }
    }

    rMin /= 2.0;
    rMax *= 2.0;
    if (verbose) std::cout << "rMin = " << rMin << ", rMax = " << rMax << std::endl;

    uint nGauLegendre = std::max(static_cast< int >(floor(6.0 * std::log10(rMax / rMin))), 4) ;
    uint nGauLaguerre = 4;
    //nGauLegendre = 10; nGauLaguerre = 10;

    initKWaveList(rMin, rMax, nGauLegendre, nGauLaguerre, kValues, weights);

    if (verbose) std::cout << "NGauLeg + NGauLag for inverse Fouriertransformation: " << nGauLegendre << " + " << nGauLaguerre << std::endl;

}

void initKWaveList(double rMin, double rMax, int nGauLegendre, int nGauLaguerre,
                   RVector & kValues, RVector & weights){

    RVector k, w;
    double k0 = 1.0 / (2.0 * rMin);

    GaussLegendre(0.0, 1.0, nGauLegendre, k, w);
    RVector kLeg(k0 * k * k);
    RVector wLeg(2.0 * k0 * k * w / PI);

    GaussLaguerre(nGauLaguerre, k, w);
    RVector kLag(k0 * (k + 1.0));
    RVector wLag(k0 * exp(k) * w / PI);

    kValues = cat(kLeg, kLag);
    weights = cat(wLeg, wLag);

//     kValues.resize(nGauLegendre + nGauLaguerre);
//     weights.resize(nGauLegendre + nGauLaguerre);

//    GaussLegendre(0.0, 1.0, nGauLegendre, k, w);
//     for (int i = 0; i < nGauLegendre; i++){
//         //** BackSubstitution der Gauss Integration
//         kValues[i] = (k0 * k[i] * k[i]);
//         //** BackSubst. WeightTERM Gauss Integration
//         weights[i] = (2.0 * k0 * k[i] * w[i] / PI);
//     }

//     GaussLaguerre(nGauLaguerre, k, w);
//     for (int i = 0; i < nGauLaguerre; i++){
//         //** BackSubstitution der Gauss Integration
//         kValues[nGauLegendre + i] = (k0 * (k[i] + 1));
//         //** BackSubst. WeightTERM Gauss Integration
//         weights[nGauLegendre + i] = (k0 * std::exp(k[i]) * w[i] / PI);
//     }
//
//     std::cout << kValues << std::endl;
//     std::cout << cat(kLeg, kLag) << std::endl;
//     std::cout << weights << std::endl;
//     std::cout << cat(wLeg, wLag) << std::endl;

}

RVector geometricFactors(const DataContainerERT & data, int dim, bool forceFlatEarth){
    RVector k(data.size());

    if (dim == 2){
        DataContainerERT tmp(data);
        for (Index i = 0; i < data.sensorCount(); i ++){
            RVector3 p(tmp.sensorPosition(i));
            // if y  differ from 0 .. we assume its x|y -> 2D .. so switch y->z
            if (p[1] != 0.0){
                p[2] = p[1];
                p[1] = 0.;
            }
            tmp.setSensorPosition(i, p);
        }
        return geometricFactors(tmp, 3, forceFlatEarth);
    }

    if (forceFlatEarth){
        DataContainerERT tmp(data);
        for (Index i = 0; i < data.sensorCount(); i ++){
            RVector3 p(tmp.sensorPosition(i)); p[2] = 0.0;
            tmp.setSensorPosition(i, p);
        }
        return geometricFactors(tmp, false);
    } else {

        for (Index i = 0; i < k.size(); i ++){
            double uam = 0.0, ubm = 0.0, uan = 0.0, ubn = 0.0;

            int a = data("a")[i];
            int b = data("b")[i];
            int m = data("m")[i];
            int n = data("n")[i];

            if (a > -1 && m > -1) uam = exactDCSolution(data.sensorPosition(a), data.sensorPosition(m));
            if (b > -1 && m > -1) ubm = exactDCSolution(data.sensorPosition(b), data.sensorPosition(m));
            if (a > -1 && n > -1) uan = exactDCSolution(data.sensorPosition(a), data.sensorPosition(n));
            if (b > -1 && n > -1) ubn = exactDCSolution(data.sensorPosition(b), data.sensorPosition(n));
//         std::cout << data(i).eA() << "   " << data(i).eB() << " " << data(i).eM()
//         << " " << data(i).eN() << std::endl;
//         std::cout << uam << " " << ubm << " " << uan << " " << ubn << std::endl;
            k[i] = 1.0 / (uam - ubm - uan + ubn);
        }
    }
    return k;
}

double exactDCSolution(const RVector3 & pot, const RVector3 & src){
    //pot, src, k, surfaceZ, fallback
    return exactDCSolution(pot, src, 0.0, 0.0, 1.0);
}

//#include <boost/math/special_functions/detail/bessel_k0.hpp>
//#include <boost/math/tr1.hpp>

double exactDCSolution(const RVector3 & v, const RVector3 & source,
                       double k, double surfaceZ, double fallback){

    double r = v.dist(source);

    if (r < TOLERANCE) return fallback;

    uint dim = 3;
    if (k > 0) dim = 2;

    RVector3 sourceMirror(source);
    sourceMirror[dim - 1] = 2.0 * surfaceZ - source[dim - 1];

    if (surfaceZ == -MAX_DOUBLE || std::isnan(surfaceZ)){// ** asuming fullspace
        if (k == 0.0) {
            //** 3D fullspace
            return 1.0 / (4.0 * PI * r);
        } else {
//             std::cout<< "Fullspace" << " " << surfaceZ << std::endl;
            return (besselK0(r * k)) / (2.0 * PI);
        }
    } else { // halfspace with mirrorspace
        if (k == 0.0) {
            //** 3D mirrorspace
//             std::cout<< "Mirrorspace" << " " << sourceMirror
//                      << " s:" << source
//                      << " v:" << v
//                      << " r:" << r << "rm:" << v.distance(sourceMirror) << std::endl;
            return (1.0 / r + 1.0 / v.distance(sourceMirror)) / (4.0 * PI) ;
        } else {
            //** just for performance issues read this again

            if (source == sourceMirror){
//                 std::cout<< "FlathEarth" << std::endl;
                //return std::tr1::cyl_bessel_k(0, r*k)/ (PI); // fac.3 slower
                return (besselK0(r * k)) / (PI);
            } else {
               // std::cout<< "Mirror" << std::endl;
                return (besselK0(r * k) +
                        besselK0(v.distance(sourceMirror) * k)) / (2.0 * PI);
            }
        }
    }
    return fallback;
}

RVector exactDCSolution(const Mesh & mesh, const RVector3 & src, double k, double surfaceZ){
    RVector solution(mesh.nodeCount());

    double fallback = 0.0;
    uint i = 0;
    for (std::vector< Node * >::const_iterator it = mesh.nodes().begin();
                it != mesh.nodes().end(); it ++){
        solution[i] = exactDCSolution((*it)->pos(), src, k, surfaceZ, fallback);
        i++;
    }
/*    for (uint i = 0; i < mesh.nodeCount(); i ++){
        solution[i] = exactDCSolution(mesh.node(i).pos(), src, k, surfaceZ, fallback);
    }*/
    return solution;
}

RVector exactDCSolution(const Mesh & mesh, const Node * nA, const Node * nB, double k, double surfaceZ){
    RVector solution;
    solution  = exactDCSolution(mesh, nA->pos(), k, surfaceZ);
    solution -= exactDCSolution(mesh, nB->pos(), k, surfaceZ);
    return solution;
}

RVector exactDCSolution(const Mesh & mesh, const ElectrodeShape * elec,
                        double k, double surfaceZ, bool setSingValue){
    RVector solution(exactDCSolution(mesh, elec->pos(), k , surfaceZ));
    if (setSingValue) elec->setSingValue(solution, 0.0, k);
    return solution;
}

RVector exactDCSolution(const Mesh & mesh, int aID, int bID,
                        double k, double surfaceZ){
    RVector solution(exactDCSolution(mesh, aID, k, surfaceZ));
    if (bID > -1) solution -= exactDCSolution(mesh, bID, k, surfaceZ);
    return solution;
}

RVector exactDCSolution(const Mesh & mesh, int aID, double k, double surfaceZ){
    return exactDCSolution(mesh, mesh.node(aID).pos(), k, surfaceZ);
}

int countKWave(const Mesh & mesh){
    RVector kValues, weights;
    initKWaveList(mesh, kValues, weights, false);
    return kValues.size();
}

void DCErrorEstimation(DataContainerERT & data, double errPerc, double errVolt, double defaultCurrent, bool verbose){
    if (verbose) std::cout << "Estimate error: " << errPerc << "% + " << errVolt << "V" << std::endl;
    RVector voltage(abs(data("u")));
    if (min(voltage) == 0.0) {
        voltage = abs(RVector(data("rhoa") / data("k")));
        if (min(data("i")) > 0.0) voltage = voltage * data("i");
        else voltage = voltage * defaultCurrent;
    }
    if (verbose) std::cout << "u min = " << min(voltage) << " V max = " << max(voltage) << " V" << std::endl;
    data.set("err", errVolt / voltage + errPerc / 100.0);
}

double DCParaDepth(const DataContainerERT & data){

    double quasiInf = 1e5;
    RVector amq(data.size(), quasiInf), anq(data.size(), quasiInf);
    RVector bmq(data.size(), quasiInf), bnq(data.size(), quasiInf);

    for (size_t i = 0; i < data.size(); i ++){
        int a = data("a")[i];
        int b = data("b")[i];
        int m = data("m")[i];
        int n = data("n")[i];

        if (a > -1 && m > -1) amq[i] = data.sensorPosition(a).distSquared(data.sensorPosition(m));
        if (b > -1 && m > -1) bmq[i] = data.sensorPosition(b).distSquared(data.sensorPosition(m));
        if (a > -1 && n > -1) anq[i] = data.sensorPosition(a).distSquared(data.sensorPosition(n));
        if (b > -1 && n > -1) bnq[i] = data.sensorPosition(b).distSquared(data.sensorPosition(n));
    }

    RVector scale((2.0 * PI) / geometricFactors(data));
    RVector tmp(data.size());
    double depth = 0.1, sens = 0.0;
    double depthQ = 4.0 * depth* depth;
    while (sens < 0.95){
        depth = depth * 1.1;
        depthQ = 4.0 * depth * depth;
        tmp = 1.0 / sqrt(bnq + depthQ) - 1.0 / sqrt(bmq + depthQ) -
              1.0 / sqrt(anq + depthQ) + 1.0 / sqrt(amq + depthQ);
        sens = sum (1.0 - tmp / scale) / data.size();
    }
    return depth;
}

void setDefaultBERTBoundaryConditions(Mesh & mesh){
    BoundingBox bbox(mesh.boundingBox());
    mesh.createNeighborInfos();
    for (uint i = 0; i < mesh.boundaryCount(); i ++){
        RVector3 cen(mesh.boundary(i).center());

        for (uint dim = 0; dim < mesh.dimension(); dim++){
            if (cen[dim] == bbox.min()[dim] || cen[dim] == bbox.max()[dim]){
                 mesh.boundary(i).setMarker(MARKER_BOUND_MIXED);
            }
        }

        if (cen[mesh.dimension() -1] == bbox.max()[mesh.dimension()-1]){
            mesh.boundary(i).setMarker(MARKER_BOUND_HOMOGEN_NEUMANN);
        }
    }
}

void setAllNeumannBoundaryConditions(Mesh & mesh){
    BoundingBox bbox(mesh.boundingBox());
    mesh.createNeighborInfos();
    for (uint i = 0; i < mesh.boundaryCount(); i ++){
        RVector3 cen(mesh.boundary(i).center());

        for (uint dim = 0; dim < mesh.dimension(); dim++){
            if (cen[dim] == bbox.min()[dim] || cen[dim] == bbox.max()[dim]){
                 mesh.boundary(i).setMarker(MARKER_BOUND_HOMOGEN_NEUMANN);
            }
        }

    }
}

RVector prepExportPotentialData(const RVector & data, double logdrop){
    RVector tmp(data);

    for (uint i = 0; i < tmp.size(); i ++) {
        tmp[i] =std::fabs(tmp[i] / logdrop);
        if (tmp[i] < 1.0) tmp[i] = 1.0;
    }

    tmp = log10(tmp);
    tmp /= max(abs(tmp)) * sign(data);
    return tmp;
}

} //namespace GIMLI

