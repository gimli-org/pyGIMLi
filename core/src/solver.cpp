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
#include "gimli.h"
#include "matrix.h"
#include "vectortemplates.h"
#include "solver.h"
#include "sparsematrix.h"

namespace GIMLI{

int solveCGLSCDWWhtrans(const MatrixBase & S, const MatrixBase & C,
                        const Vec & dWeight, 
                        const Vec & b, // deltaData
                        Vec & x,       // deltaModel
                        const Vec & wc, const Vec & wm,
                        const Vec & tm, // m/d mod
                        const Vec & td, // d/d resp
                        double lambda, const Vec & roughness,
                        int maxIter, double tol,
                        bool verbose){ //ALLOW_PYTHON_THREADS
    
    uint nData = b.size();
    uint nModel = x.size();
    uint nConst = C.rows();
    Vec bR(b);
    Vec dW(dWeight);
    // bR.round(1e-10);
    // dW.round(1e-10);
    

//     __MS(S.rtti())
//     __MS(C.rtti())

    if (S.rows() != nData)  std::cerr << "J.rows != nData " << S.rows() << " / " << nData << std::endl;
    //if (S.cols() != nModel) std::cerr << "J.cols != nModel " << S.cols() << " / " << nModel << std::endl;
    if (C.cols() != nModel) std::cerr << "C.cols != nModel " << C.cols() << " / " << nModel << std::endl;
    if (dW.size() != nData) std::cerr << "dW.size() != nData" << dW.size() << " != " << nData << std::endl;
    if (wc.size() != nConst) std::cerr << "wc.size() != nConst " << wc.size() << " / " << nConst << std::endl;
    if (wm.size() != nModel) std::cerr << "wm.size() != nModel" << wm.size() << " / " << nModel << std::endl;
    if (tm.size() != nModel) std::cerr << "tm.size() != nModel " << tm.size() << " / " << nModel << std::endl;
    if (td.size() != nData) std::cerr << "td.size() != nData " << td.size() << " / " << nData << std::endl;
    if (roughness.size() != nConst) std::cerr << "roughness.size != nConst " << roughness.size() << " / " << nConst << std::endl;

//Ch  Vec cdx(transMult(C, Vec(wc * wc * (C * Vec(wm * deltaX)))) * wm * lambda); // nModel
    Vec cdx(transMult(C, Vec(wc * roughness)) * wm * lambda); // nModel
    Vec z((bR - S * Vec(x / tm) * td) * dW); // nData
    Vec p(transMult(S, Vec(z * dW * td)) / tm   
          - transMult(C, Vec(wc * wc * (C * Vec(wm * x)))) * wm * lambda
          - cdx );// nModel
    Vec r(transMult(S, Vec(bR * dW * dW * td)) / tm - cdx); // nModel
    
    // p = round(p, 1e-10);
    // r = round(r, 1e-10);
// std::cout.precision(14);
// __MS("tol " << tol)
// __MS("dw " << min(dW) << " " << max(dW) << " " << mean(dW))
// __MS("b " << min(bR) << " " << max(bR) << " " << mean(bR))
// __MS("x " << min(x) << " " << max(x) << " " << mean(x))
// __MS("wc " << min(wc) << " " << max(wc) << " " << mean(wc))
// __MS("wm " << min(wm) << " " << max(wm) << " " << mean(wm))
// __MS("tm " << min(tm) << " " << max(tm) << " " << mean(tm))
// __MS("td " << min(td) << " " << max(td) << " " << mean(td))
// __MS("cdx " << min(cdx) << " " << max(cdx) << " " << mean(cdx))
// __MS("(x/tm) " << min((x/tm)) << " " << max((x/tm)) << " " << mean((x/tm)))
// __MS("S*Vec(x/tm) " << min(S*Vec(x/tm)) << " " << max(S*Vec(x/tm)) << " " << mean(S*Vec(x/tm)))
// __MS("z " << min(z) << " " << max(z) << " " << mean(z))
// __MS("p " << min(p) << " " << max(p) << " " << mean(p))
// __MS("r" << min(r) << " " << max(r) << " " << mean(r))
// __MS(min(Vec(bR * dW * dW * td)) << " " 
//      << max(Vec(bR * dW * dW * td)) << " " 
//      << mean(Vec(bR * dW * dW * td)))
// __MS(min(transMult(S,Vec(bR.size(), 1))) << " " 
//      << max(transMult(S,Vec(bR.size(), 1))) << " " 
//      << mean(transMult(S,Vec(bR.size(), 1))))
// // if (z.size() > 100 )exit(1);

//     RVector sv(dynamic_cast< const RSparseMapMatrix *> (&S)->values());
//     __MS("S:" << " " << min(sv) << " " << max(sv) << " " << mean(sv))
//     RVector cv(dynamic_cast< const RSparseMapMatrix *> (&C)->values());
//     __MS("C:" << " " << min(cv) << " " << max(cv) << " " << mean(cv))
 
    double accuracy = tol;
    if (accuracy < 0.0) accuracy = max(TOLERANCE, 1e-08 * dot(r, r));
    r = p;

    double normR2 = dot(r, r), normR2old = 0.0;
    double alpha = 0.0, beta = 0.0;
    double aQ = 0.0;
    int count = 0;

    Vec q(nData, 0.0);
    Vec wcp(nConst, 0.0); // nBounds

//     std::cout.precision(14);
    // std::cout << "###############################################" << std::endl;
    // std::cout << 0 << "  " << accuracy << std::endl;

    // hack until matMult with expression templates
    Vec pwm(nConst, 0.0);
    Vec p_tm(nModel, 0.0);
    Vec x_tm(nModel, 0.0);
    Vec zdWtd(nData, 0.0);
    Vec wmx(nConst, 0.0);
    Vec wcwcCwmx(nConst, 0.0);
    Vec Cpwm(nConst, 0.0);

    while (count < maxIter && normR2 > accuracy){
        count ++;
        p_tm.assign(p / tm);
        q = (S * p_tm) * (dW * td);
        // q = S * Vec(p / tm) * dW * td;
        // std::cout << "q " << min(q) << " " << max(q) << " " << mean(q) << std::endl;
        // q = round(q, 1e-10);
        // std::cout << "q " << min(q) << " " << max(q) << " " << mean(q) << std::endl;
        pwm.assign(p * wm);
        Cpwm.assign(C * pwm);
        
        //** try to avoid accuracy problems with unsorted C
        Cpwm.round(1e-10);

        // Vec Cpwm((C * pwm));
        wcp = wc * Cpwm;
        
        // std::cout << "pwm " << min(pwm) << " " << max(pwm) << " " << mean(pwm) << std::endl;
        // std::cout << "Cpwm " << min(Cpwm) << " " << max(Cpwm) << " " << mean(Cpwm) << std::endl;
        // std::cout << "wcp " << min(wcp) << " " << max(wcp) << " " << mean(wcp) << std::endl;

        if (count == 1 && std::fabs(mean(Cpwm)+17.00029) < 1e-2){
            // C.save("C17000.bmat");
            // pwm.save("pwm1700.bmat");
            //throwError("debug halt");
__MS("##################################################################")
        }
        // if (count == 1 && std::fabs(mean(Cpwm)+17.708) < 1e-2){
        //     // C.save("C17708.bmat");
        //     // pwm.save("pwm17708.bmat");
        //     throwError("debug halt");
        // }
        // wcp = wc * (C * pwm);
        // wcp = round(wcp, 1e-10);


        aQ = (dot(q, q) + lambda * dot(wcp, wcp));
        // if (aQ < 1e-10){
        //     aQ = 1e-10;
        // }
        
        // std::cout << "aQ " << aQ<< std::endl;
        alpha = normR2 / aQ;
        //alpha = std::round(alpha * 1e10) / 1e10;

        x += p * alpha;

        if ((count % 10) == -1) { // TOM
            x_tm.assign(x / tm);
            z = bR - (S * x_tm) * td; // z exakt durch extra Mult.
            z *= dW;
        } else {
            z -= q * alpha;
        }
        //z = round(z, 1e-10);

        zdWtd.assign(z * dW * td);
        wmx.assign(wm * x);
        wcwcCwmx.assign(wc * wc * (C * wmx));

        r =   transMult(S, zdWtd) / tm 
            - transMult(C, wcwcCwmx) * wm * lambda - cdx;
        // r = transMult(S, Vec(z * dW * td)) / tm 
        //     - transMult(C, Vec(wc * wc * (C * Vec(wm * x)))) * wm * lambda 
        //     - cdx;
        //r = round(r, 1e-8);

        normR2old = normR2;
        normR2 = dot(r, r);
        beta = normR2 / normR2old;
        //beta = std::round(beta * 1e10) / 1e10;

        p = r + p * beta;

//         std::cout << "\n" << count << " +++++"  << " " << normR2<< std::endl;
//         std::cout << alpha << " " << beta<< " " << alpha << std::endl;
//         std::cout << "q " << min(q) << " " << max(q) << " " << mean(q) << std::endl;
//         std::cout << "wcp "  << min(wcp) << " " << max(wcp) << " " << mean(wcp) << std::endl;
//         std::cout << "x "  << min(x) << " " << max(x) << " " << mean(x) << std::endl;
//         std::cout << "r "  << min(r) << " " << max(r) << " " << mean(r) << std::endl;
//         std::cout << "p "  << min(p) << " " << max(p) << " " << mean(p) << std::endl;
        
// __MS("---" << count << "-------" << normR2 << "-------------------------")
// __MS(min(q) << " " << max(q) << " " << mean(q))
// __MS(min(wcp) << " " << max(wcp) << " " << mean(wcp))
// __MS(alpha)

//         if((count < 200)) {
//             std::cout << count << "  " << normR2 << std::endl;
//         }

#ifndef _WIN32
        if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
#endif
    }
#ifdef _WIN32
    if (verbose) std::cout << "[ " << count << "/" << normR2 << "]\t" << std::endl;
#endif

    // if (z.size() > 100 )exit(1);
    return 1;
}

int solveCGLSCDWWtrans(const MatrixBase & S, const MatrixBase & C,
                       const Vec & dWeight,  const Vec & b, Vec & x,
                       const Vec & wc, const Vec & mc, const Vec & tm,
                       const Vec & td, double lambda, const Vec & deltaX,
                       int maxIter, bool verbose){ //ALLOW_PYTHON_THREADS
// __M
  uint nData = b.size();
  uint nModel = x.size();
  uint nConst = C.rows();

  if (S.rows() != nData)  std::cerr << "J.rows != nData " << S.rows() << " / " << nData << std::endl;
  // not needed
  // if (S.cols() != nModel) std::cerr << "J.cols != nModel " << S.cols() << " / " << nModel << std::endl;
  if (C.rows() != nConst) std::cerr << "C.rows != nConst " << C.rows() << " / " << nConst << std::endl;
  if (C.cols() != nModel) std::cerr << "C.cols != nModel " << C.cols() << " / " << nModel << std::endl;
  if (dWeight.size() != nData) std::cerr << "dWeight.size() != nData" << dWeight.size() << " != " << nData << std::endl;
  if (wc.size() != nConst) std::cerr << "wc.size() != nConst " << wc.size() << " / " << nConst << std::endl;
  if (mc.size() != nModel) std::cerr << "mc.size() != nModel" << mc.size() << " / " << nModel << std::endl;
  if (tm.size() != nModel) std::cerr << "tm.size() != nModel " << tm.size() << " / " << nModel << std::endl;
  if (td.size() != nData) std::cerr << "td.size() != nData " << td.size() << " / " << nData << std::endl;

  Vec cdx(transMult(C, Vec(wc * wc * (C * Vec(mc * deltaX)))) * mc * lambda); // nModel
  Vec z((b - S * Vec(x / tm) * td) * dWeight); // nData
  Vec p(transMult(S, Vec(z * dWeight * td)) / tm - cdx - transMult(C, Vec(wc * wc * (C * Vec(mc * x)))) * mc * lambda);// nModel
  Vec r(transMult(S, Vec(b * dWeight * dWeight * td)) / tm - cdx); // nModel

  double accuracy = max(TOLERANCE, 1e-08 * dot(r, r));
  r = p;

  double normR2 = dot(r, r), normR2old = 0.0;
  double alpha = 0.0, beta = 0.0;

  int count = 0;

  Vec q(nData);
  Vec wcp(nConst); // nBounds

  while (count < maxIter && normR2 > accuracy){

    count ++;
    q = S * Vec(p / tm) * dWeight * td;
    wcp = wc * (C * Vec(p * mc));

    alpha = normR2 / (dot(q, q) + lambda * dot(wcp, wcp));
    x += p * alpha;
    if((count % 10) == -1) { // TOM
      z = b - S * Vec(x / tm) * td; // z exakt durch extra Mult.
      z *= dWeight;
    } else {
      z -= q * alpha;
    }
    r = transMult(S, Vec(z * dWeight * td)) / tm - transMult(C, Vec(wc * wc * (C * Vec(mc * x)))) * mc * lambda - cdx;

    normR2old = normR2;
    normR2 = dot(r, r);
    beta = normR2 / normR2old;
    p = r + p * beta;
    #ifndef _WIN32
    if (verbose) std::cout << "\r[ " << count << "/" << normR2 << "]\t";
    #endif
  }
  #ifdef _WIN32
  if (verbose) std::cout << "[ " << count << "/" << normR2 << "]\t" << std::endl;
  #endif
  return 1;
}

} //namespace GIMLI

