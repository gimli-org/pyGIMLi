/******************************************************************************
 *   Copyright (C) 2009-2024 by the GIMLi development team                    *
 *   Thomas Günther thomas@resistivity.net                                    *
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

#include "gimli.h"
#include "em1dmodelling.h"
#include "meshgenerators.h"

#include <math.h>

namespace GIMLI {

RVector MT1dModelling::rhoaphi(const RVector & rho, const RVector & thk) { // after mtmod.c by R.-U. B�rner
    size_t nperiods = periods_.size();
    RVector rhoa(nperiods), phi(nperiods);

    RVector::ValType my0 = PI * 4e-7;
    Complex i_unit(0.0 , 1.0), adm, alpha, tanalpha;
    CVector z(nlay_);
    for (size_t i = 0 ; i < nperiods ; i++) {
        RVector::ValType omega = 2.0 * PI / periods_[i];
        z[nlay_ - 1] = sqrt(i_unit * omega * rho[nlay_ - 1] / my0);
        for (int k = nlay_ - 2 ; k >= 0 ; k--) {
            adm = sqrt(my0 / (rho[k] * i_unit * omega));
            alpha = thk[k] * sqrt(i_unit * my0 * omega / rho[k]);
            tanalpha = sinh(alpha) / cosh(alpha);
            z[k] = (adm * z[k + 1] + tanalpha) / (adm * z[k + 1] * tanalpha + (RVector::ValType)1.0);
            z[k] /= adm;
        }
        rhoa[i] = abs(z[0]) * abs(z[0]) * my0 / omega;
        phi[i] = std::atan(imag(z[0]) / real(z[0]));
    }
    return cat(rhoa, phi);
}

RVector MT1dModelling::rhoa(const RVector & model){ //! app. res. for thk/res vector
    if (model.size() != nlay_ * 2 - 1) {
	__M
return 0;
}
    RVector thk(model, 0, nlay_ - 1), rho(model, nlay_ - 1, 2 * nlay_ - 1);
    return rhoa(rho, thk);
}

    /*! the actual (full) forward operator returning app.res.+phase for thickness+resistivity */
RVector MT1dModelling::response(const RVector & model) {
    if (model.size() != nlay_ * 2 - 1) {
__M
return 0;
}
    RVector thk(model, 0, nlay_ - 1), rho(model, nlay_ - 1, 2 * nlay_ - 1);
    return rhoaphi(rho, thk);
}

RVector MT1dModelling::rhoa(const RVector & rho, const RVector & thk) { // after mtmod.c by R.-U. B�rner
    size_t nperiods = periods_.size();
    return rhoaphi(rho, thk)(0,nperiods);
    //CR not needed?
//     RVector rhoa(nperiods);
//     static double my0 = PI * 4e-7;
//     Complex i_unit(0.0 , 1.0), adm, alpha, tanalpha;
//     CVector z(nlay_);
//     for (size_t i = 0 ; i < nperiods ; i++) {
//         double omega = 2.0 * PI / periods_[i];
//         z[nlay_ - 1] = sqrt(i_unit * omega * rho[nlay_ - 1] / my0);
//         for (int k = nlay_ - 2 ; k >= 0 ; k--) {
//             adm = sqrt(my0 / (rho[k] * i_unit * omega));
//             alpha = thk[k] * sqrt(i_unit * my0 * omega / rho[k]);
//             tanalpha = sinh(alpha) / cosh(alpha);
//             z[k] = (adm * z[k + 1] + tanalpha) / (adm * z[k + 1] * tanalpha + 1.0);
//             z[k] /= adm;
//         }
//         rhoa[i] = abs(z[0]) * abs(z[0]) * my0 / omega;
//     }
//     return rhoa;
}

FDEM1dModelling::FDEM1dModelling(size_t nlay,
                                 const RVector & freqs,
                                 const RVector & coilspacing,
                                 double z, bool verbose)
    : ModellingBase(verbose), nlay_(nlay), freqs_(freqs),
      coilspacing_(coilspacing),
      zs_(-std::fabs(z)), ze_(-std::fabs(z)) {
    init();
}

FDEM1dModelling::FDEM1dModelling(size_t nlay,
                                 const RVector & freqs,
                                 double coilspacing,
                                 double z, bool verbose)
    : ModellingBase(verbose), nlay_(nlay), freqs_(freqs),
      zs_(-std::fabs(z)), ze_(-std::fabs(z)) {

    coilspacing_ = RVector(freqs.size(), coilspacing);
    init();
}

void FDEM1dModelling::init(){
    setMesh(createMesh1DBlock(nlay_));
    nfr_ = freqs_.size();
    // double zp = ze_ + zs_;  // indeed a bad mistake for airborne
    double zpq = (ze_ - zs_) * (ze_ - zs_);
    RVector rpq(coilspacing_ * coilspacing_ + zpq);
    freeAirSolution_ = (rpq - zpq * 3.0) / rpq / rpq / sqrt(rpq) / 4.0 / PI;
}

Complex btp(double u, double f, RVector rho, RVector d){
    size_t nl = rho.size();
    double mu0 = 4e-7 * PI;
    Complex c(0.0, mu0 * 2. * PI * f);

    Complex b(std::sqrt(c / rho[nl-1] + u*u));
    if(nl > 1) {
        for(int nn=nl-2; nn>=0 ; nn--){
            Complex alpha(std::sqrt(c/rho[nn] + u*u));
            Complex cth(std::exp(alpha * d[nn] * -2.0));
            cth=(Complex(1.0) - cth) / (cth + 1.0);
            b=(alpha * cth + b) / (cth * b / alpha + 1.0);
        }
    }
    return b;
}

RVector FDEM1dModelling::calc(const RVector & rho, const RVector & thk){
    double hankelJ0[100]={
        2.89878288E-07,3.64935144E-07,4.59426126E-07,5.78383226E-07,
        7.28141338E-07,9.16675639E-07,1.15402625E-06,1.45283298E-06,
        1.82900834E-06,2.30258511E-06,2.89878286E-06,3.64935148E-06,
        4.59426119E-06,5.78383236E-06,7.28141322E-06,9.16675664E-06,
        1.15402621E-05,1.45283305E-05,1.82900824E-05,2.30258527E-05,
        2.89878259E-05,3.64935186E-05,4.59426051E-05,5.78383329E-05,
        7.28141144E-05,9.16675882E-05,1.15402573E-04,1.45283354E-04,
        1.82900694E-04,2.30258630E-04,2.89877891E-04,3.64935362E-04,
        4.59424960E-04,5.78383437E-04,7.28137738E-04,9.16674828E-04,
        1.15401453E-03,1.45282561E-03,1.82896826E-03,2.30254535E-03,
        2.89863979E-03,3.64916703E-03,4.59373308E-03,5.78303238E-03,
        7.27941497E-03,9.16340705E-03,1.15325691E-02,1.45145832E-02,
        1.82601199E-02,2.29701042E-02,2.88702619E-02,3.62691810E-02,
        4.54794031E-02,5.69408192E-02,7.09873072E-02,8.80995426E-02,
        1.08223889E-01,1.31250483E-01,1.55055715E-01,1.76371506E-01,
        1.85627738E-01,1.69778044E-01,1.03405245E-01,-3.02583233E-02,
        -2.27574393E-01,-3.62173217E-01,-2.05500446E-01,3.37394873E-01,
        3.17689897E-01,-5.13762160E-01,3.09130264E-01,-1.26757592E-01,
        4.61967890E-02,-1.80968674E-02,8.35426050E-03,-4.47368304E-03,
        2.61974783E-03,-1.60171357E-03,9.97717882E-04,-6.26275815E-04,
        3.94338818E-04,-2.48606354E-04,1.56808604E-04,-9.89266288E-05,
        6.24152398E-05,-3.93805393E-05,2.48472358E-05,-1.56774945E-05,
        9.89181741E-06,-6.24131160E-06,3.93800058E-06,-2.48471018E-06,
        1.56774609E-06,-9.89180896E-07,6.24130948E-07,-3.93800005E-07,
        2.48471005E-07,-1.56774605E-07,9.89180888E-08,-6.24130946E-08};

    //** extract resistivity and thickness
    RVector inph(nfr_), outph(nfr_);//** inphase and quadrature components
    int nc = 100, nc0 = 60; // number of coefficients
    RVector::ValType q=0.1*std::log(10.0);
    for (size_t i = 0 ; i < nfr_ ; i++) {
        Complex aux(0.0, 0.0);
        for (int ii = 0 ; ii < nc ; ii++) {
            RVector::ValType ui=std::exp(q * (nc - ii - nc0)) / coilspacing_[i];
            Complex bti(btp(ui, freqs_[i], rho, thk));
            Complex delta((bti - ui) / (bti + ui) * std::exp(ui * ze_) * std::exp(ui * zs_));
            aux += delta * ui * ui * hankelJ0[nc - ii - 1];
        }

        aux /= PI * 4.0 * coilspacing_[i];
        // normalize by free air solution in per cent
        inph[i]  = real(aux) / freeAirSolution_[i] * 100.0;
        outph[i] = imag(aux) / freeAirSolution_[i] * 100.0;
    }
    //** paste together both components and
    return cat(inph, outph);
}

RVector FDEM1dModelling::response(const RVector & model){
    RVector thk(model, 0, nlay_ - 1), rho(model, nlay_ - 1, 2 * nlay_ - 1);
    return calc(rho, thk);
}

RVector MRSModelling::response(const RVector & model) {
    RVector outreal(*KR_ * model);
    RVector outimag(*KI_ * model);
    return RVector(sqrt(outreal * outreal + outimag * outimag));
}

void MRSModelling::createJacobian(const RVector & model) {
    RVector ddr(*KR_ * model);
    RVector ddi(*KI_ * model);
    RVector dda(sqrt(ddr * ddr + ddi * ddi));

    RMatrix * jacobian = dynamic_cast < RMatrix * >(jacobian_);
    jacobian->resize(dda.size(), model.size());

    for (size_t i = 0 ; i < KR_->rows() ; i++) {
        (*jacobian)[i] = ((*KR_)[i] * ddr[i] + (*KI_)[i] * ddi[i]) / dda[i];
    }
}

RVector MRS1dBlockModelling::response(const RVector & model){
        //! extract water content and thickness from model vector
        RVector wc(model, nlay_ - 1 , nlay_ * 2 - 1);
        RVector thk(model, 0 , nlay_ - 1);
        //! fill vector of original size wit last layer water content
        RVector wcvec(nvec_, wc[nlay_ - 1]);
        size_t iz1 = 0, iz2 = 0;
        double zthk = 0;
        //! run through layers and fill water content
        for (size_t i = 0 ; i < nlay_ - 1 ; i++){
            zthk += thk[i];
            iz2 = 0;
            while (iz2 < zvec_.size() && zvec_[iz2] < zthk) iz2++;
            if (iz2 > nvec_) iz2 = nvec_;
            for (size_t j = iz1 ; j < iz2 ; j++) wcvec[j] = wc[i];
            if (iz2 + 1 >= zvec_.size()) break; // end reached
            wcvec[iz2] = ((zthk - zvec_[iz2]) * wc[i]
                           + (zvec_[iz2 + 1] - zthk) * wc[i + 1])
                         / (zvec_[iz2 + 1] - zvec_[iz2]);
            iz1 = iz2 + 1;
        }

        if (verbose_) save(wcvec, "wctmp.vec");
        //! call original forward response and return;
        return MRSModelling::response(wcvec);
    }


} // namespace GIMLI{
