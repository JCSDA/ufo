/*
 * (C) Copyright 2024-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Author: Rong Kong (NCAR/MMM) - C++ rewrite of the operator from Fortran;
 *         integration into the PPRO library (2025).
 *
 * Original Fortran author: Tzu-Chin Tsai (CWA, original implementation).
 */

#include "ufo/operators/radarpolarimetric/ppro/TCWA2Forward.h"

#include <algorithm>
#include <cmath>

namespace ufo {
namespace ppro {

namespace {
// Shared gamma-PSD afa fitting coefficients.
constexpr double afamin = 0.0, afamax = 30.0;
constexpr double afa1 = -0.35897435,   afa2 = 0.010041278;
constexpr double afa3 = 1.9349749e-3,  afa4 = -0.067453607;
constexpr double afa5 = -5.5693923e-3, afa6 = 0.06376154;
constexpr double afb1 = -6.5100404,    afb2 = 0.039597245;
constexpr double afb3 = 2.3169628e-3,  afb4 = -0.45856646;
constexpr double afb5 = -0.04537631,   afb6 = 0.25713479;
constexpr double thrd = 1.0 / 3.0;
constexpr double zxmin = 1.0e-21;
}  // namespace

// ----------------------------------------------------------------------------
// ln(Gamma(xx)) for xx > 0.  Uses the C++ standard library log-gamma; all
// arguments in the TCWA2 routines are strictly positive (afa + {1,4,7}).
double tcwa2Gamln(double xx) {
  return std::lgamma(xx);
}

// ----------------------------------------------------------------------------
void dualpolOpRainTcwa2(const int iband, const double rho, const double qc0,
                        const double qr0, const double nr0,
                        double &zh, double &zv, double &kdp) {
  const double drmin = 1.0e-4, drmax = 6.0e-3;
  zh = zxmin;
  zv = zxmin;
  kdp = 0.0;
  double qc = qc0 / 1.0e3;   // g/kg -> kg/kg
  double qr = qr0 / 1.0e3;   // g/kg -> kg/kg
  double nr = nr0 / rho;     // 1/m^3 -> 1/kg
  const double dbzwl = (iband == 1) ? 0.11 : 0.053;

  if (qr >= 1.0e-14 && nr >= 1.0e-2) {
    double mvrr = std::min(drmax / 2.0,
                  std::max(drmin / 2.0, std::pow(qr / nr / 4.18879e3, thrd)));
    double afar;
    if (qc >= 1.0e-14 || nr < 0.5) {
      afar = std::min(afamax, std::max(afamin,
             (afa1 + afa2 * std::log(nr) + afa3 * std::pow(std::log(nr), 2.0)
              + afa4 * std::log(mvrr))
             / (1.0 + afa5 * std::log(nr) + afa6 * std::log(mvrr))));
    } else {
      afar = std::min(afamax, std::max(afamin,
             (afb1 + afb2 * std::log(nr) + afb3 * std::pow(std::log(nr), 2.0)
              + afb4 * std::log(mvrr))
             / (1.0 + afb5 * std::log(nr) + afb6 * std::log(mvrr))));
    }
    if ((nr * std::pow(mvrr, 6.0)) >= 1.0e-16) {
      if (nr >= 3.8e3) {
        afar = afamax;
      } else {
        afar = 0.5 * afar;
      }
    }
    const double gr1 = tcwa2Gamln(afar + 1.0);
    const double gr4 = tcwa2Gamln(afar + 4.0);
    double lamr = std::pow(std::exp(gr4 - gr1) * nr * 523.5987756 / qr, thrd);
    const double lamrmin = std::pow(std::exp(gr4 - gr1), thrd) / drmax;
    const double lamrmax = std::pow(std::exp(gr4 - gr1), thrd) / drmin;
    if (lamr < lamrmin) {
      lamr = lamrmin;
      nr = qr * 1.909859e-3 * std::exp(gr1 - gr4 + 3.0 * std::log(lamr));
    } else if (lamr > lamrmax) {
      lamr = lamrmax;
      nr = qr * 1.909859e-3 * std::exp(gr1 - gr4 + 3.0 * std::log(lamr));
    }
    const double mvdr = std::pow(std::exp(gr4 - gr1), thrd) / lamr;
    if (qr >= 1.0e-6) {
      const double llmr = std::log(lamr);
      const double lar2 = std::log(afar + 2.0);
      const double llmr2 = llmr * llmr;
      const double llmr3 = llmr2 * llmr;
      const double lar22 = lar2 * lar2;
      const double lar23 = lar22 * lar2;
      const double fdbzr = 1950.2082 + 29.685692 * lar2 - 17.761215 * std::sqrt(lar2)
                         + 328.61414 * std::pow(std::log(llmr), 2.0)
                         - 1122.8963 * std::sqrt(llmr);
      zh = nr * std::pow(10.0, (fdbzr - 200.0) / 10.0 - 18.0);
      const double dbzr = 10.0 * std::log10(1.0e18 * zh);
      double fzdrr;
      if (afar >= 3.0) {
        fzdrr = std::exp(27.639873 + 3.7680268 * lar2 - 4.3382556 * llmr
                + 0.36498072 * lar22 + 0.29286255 * llmr2 - 0.58159197 * lar2 * llmr
                + 9.968596e-3 * lar23 - 7.1806791e-3 * llmr3
                + 0.023642573 * lar2 * llmr2 - 0.028767891 * lar22 * llmr) / 1.0e4;
      } else {
        fzdrr = std::exp(10.243103 + 0.46561642 * lar2 + 1.7270962 * llmr
                + 5.8974032e-3 * lar22 - 0.34268382 * llmr2 + 0.10193409 * lar2 * llmr
                - 0.11005788 * lar23 + 0.015622039 * llmr3
                - 0.038623263 * lar2 * llmr2 + 0.10413007 * lar22 * llmr) / 1.0e4;
      }
      zv = std::pow(10.0, dbzr / 10.0 - 18.0 - fzdrr / 10.0);
      double fkdpr;
      if (mvdr < 7.0e-4) {
        fkdpr = std::min(1.5 / std::max(1.0, nr),
                std::exp(14.837602 - 18.786496 / lar2 + 72.211021 / lar22
                - 124.16879 / lar23 + 98.39131 / std::pow(lar2, 4.0)
                - 28.355229 / std::pow(lar2, 5.0) - 0.092813834 * llmr) / 1.0e10);
      } else {
        if (afar <= 3.0) {
          fkdpr = std::exp(35.450715 + 2.4505311 * lar2 - 1.4651261 * llmr
                  - 0.22974058 * lar22 - 0.22139098 * llmr2
                  + 0.29382343 * lar2 * llmr) / 1.0e10;
        } else {
          fkdpr = std::exp((38.35738 + 8.8170679 * lar2 - 8.971269 * llmr
                  + 0.49634922 * lar22 + 0.51897478 * llmr2 - 1.0162854 * lar2 * llmr)
                  / (1.0 + 0.16321213 * lar2 - 0.1676006 * llmr
                  + 4.0072333e-3 * lar22 + 4.9701831e-3 * llmr2
                  - 9.1512549e-3 * lar2 * llmr)) / 1.0e10;
        }
        fkdpr = std::min(10.0 / std::max(1.0, nr), fkdpr);
      }
      kdp = nr * fkdpr * 0.11 / dbzwl;
      if (zv > zh) zv = 0.9999 * zh;
    }
  }
}

// ----------------------------------------------------------------------------
void dualpolOpIceTcwa2(const int iband, const double tk, const double rho,
                       const double qc0, const double qi0, const double ni0,
                       double &zh, double &zv, double &kdp) {
  const double dimin = 6.0e-6, dimax = 5.0e-3;
  static const double itble[121] = {
    1.000000, 0.979490, 0.959401, 0.939723, 0.920450, 0.899498,
    0.879023, 0.857038, 0.833681, 0.810961, 0.783430, 0.755092,
    0.703072, 0.537032, 0.467735, 0.524807, 0.630957, 0.812831,
    1.096478, 1.479108, 1.905461, 2.089296, 2.290868, 2.398833,
    2.454709, 2.426610, 2.371374, 2.290868, 2.137962, 1.995262,
    1.862087, 1.737801, 1.621810, 1.513561, 1.396368, 1.288250,
    1.188502, 1.096478, 1.000000, 0.922571, 0.851138, 0.785236,
    0.724436, 0.668344, 0.616595, 0.575440, 0.537032, 0.501187,
    0.467735, 0.436516, 0.407380, 0.380189, 0.354813, 0.331131,
    0.316228, 0.301995, 0.291743, 0.285102, 0.281838, 0.278612,
    0.275423, 0.278612, 0.281838, 0.285102, 0.291743, 0.298538,
    0.309030, 0.319890, 0.331131, 0.346737, 0.367282, 0.393550,
    0.426580, 0.457088, 0.489779, 0.524807, 0.562341, 0.609537,
    0.660693, 0.716143, 0.785236, 0.860994, 0.954993, 1.047129,
    1.148154, 1.258925, 1.380384, 1.496236, 1.603245, 1.698244,
    1.778279, 1.840772, 1.883649, 1.905461, 1.905461, 1.883649,
    1.862087, 1.840772, 1.798871, 1.737801, 1.698244, 1.640590,
    1.584893, 1.548817, 1.513561, 1.475707, 1.452112, 1.428894,
    1.412538, 1.393157, 1.377209, 1.361445, 1.348963, 1.336596,
    1.327394, 1.318257, 1.309182, 1.303167, 1.294196, 1.288250,
    1.279381};

  zh = zxmin;
  zv = zxmin;
  kdp = 0.0;
  double qc = qc0 / 1.0e3;
  double qi = qi0 / 1.0e3;
  double ni = ni0 / rho;
  const double dbzwl = (iband == 1) ? 0.11 : 0.053;

  if (qi >= 1.0e-14 && ni >= 1.0e-2) {
    const double rhoi = 900.0;
    double mvri = std::min(dimax / 2.0,
                  std::max(dimin / 2.0, std::pow(1.90986 * qi / ni / rhoi, thrd)));
    double afai;
    if (qc >= 1.0e-14 || ni < 0.5) {
      afai = std::min(afamax, std::max(afamin,
             (afa1 + afa2 * std::log(ni) + afa3 * std::pow(std::log(ni), 2.0)
              + afa4 * std::log(mvri))
             / (1.0 + afa5 * std::log(ni) + afa6 * std::log(mvri))));
    } else {
      afai = std::min(afamax, std::max(afamin,
             (afb1 + afb2 * std::log(ni) + afb3 * std::pow(std::log(ni), 2.0)
              + afb4 * std::log(mvri))
             / (1.0 + afb5 * std::log(ni) + afb6 * std::log(mvri))));
    }
    const double gi1 = tcwa2Gamln(afai + 1.0);
    const double gi4 = tcwa2Gamln(afai + 4.0);
    const double i3m = 1.90986 * qi / rhoi;
    double lami = std::pow(std::exp(gi4 - gi1) * ni / i3m, thrd);
    const double lamimin = std::pow(std::exp(gi4 - gi1), thrd) / dimax;
    const double lamimax = std::pow(std::exp(gi4 - gi1), thrd) / dimin;
    if (lami < lamimin) {
      lami = lamimin;
      ni = i3m * std::exp(gi1 - gi4 + 3.0 * std::log(lami));
    } else if (lami > lamimax) {
      lami = lamimax;
      ni = i3m * std::exp(gi1 - gi4 + 3.0 * std::log(lami));
    }
    const double mvdi = std::pow(std::exp(gi4 - gi1), thrd) / lami;
    double adagr;
    if (mvdi >= 6.0e-6) {
      const double lbdi = std::log(mvdi * 1.0e6);
      const double ttr = std::min(1.0, std::max(0.01,
                         9.207e-3 * std::pow(lbdi, 3.0) - 0.150212 * std::pow(lbdi, 2.0)
                         + 0.565632 * lbdi + 0.401116));
      const double tc = tk - 273.15;
      const int hid = std::max(std::min(
          static_cast<int>(std::lround(std::fabs(tc) / 0.25)), 120), 0);
      const double inhgr = itble[hid];
      adagr = std::pow(inhgr, ttr);
    } else {
      adagr = 1.0;
    }
    const double zeta = (adagr - 1.0) / (adagr + 2.0);
    if (qi >= 1.0e-6) {
      const double llmi = std::log(lami);
      const double lai2 = std::log(afai + 2.0);
      const double lroi = std::log(rhoi);
      const double lmdi = std::log(1.0e6 * mvdi);
      const double llmi2 = llmi * llmi;
      const double lai22 = lai2 * lai2;
      const double lai23 = lai22 * lai2;
      const double lmdi2 = lmdi * lmdi;
      const double lmdi3 = lmdi2 * lmdi;
      const double lroi2 = lroi * lroi;
      const double lroi3 = lroi2 * lroi;
      const double zeta2 = zeta * zeta;
      const double zeta3 = zeta2 * zeta;
      if ((adagr - 1.0) >= 1.0e-2) {
        const double fdbzi = 389.04082 + 18.161704 * lai2 + 1.3401764 * lai22
                           - 0.07360894 * lai23 - 26.137103 * llmi;
        const double zhrox = 1.9382e-6 * std::exp(1.9260732373 * lroi);
        const double zhzet = 1.0;
        zh = ni * std::pow(10.0, (fdbzi - 200.0) / 10.0 - 18.0) * zhrox * zhzet;
        const double dbzi = 10.0 * std::log10(1.0e18 * zh);
        const double zdrox = 0.0585053 * lroi3 - 0.80223023 * lroi2
                           + 3.75756297 * lroi - 5.86526003;
        const double zdzet = std::max(0.0, -0.01644484 + 0.010454643 * lmdi
                           + 9.4367973 * zeta - 0.0026137203 * lmdi2 + 5.7453232 * zeta2
                           + 0.3234453 * lmdi * zeta + 2.0622912e-4 * lmdi3
                           - 1.3945604 * zeta3 - 2.2879493 * lmdi * zeta2
                           - 0.013856582 * lmdi2 * zeta);
        const double fzdri = std::exp((15.181145 + 0.34794935 * lai2
                           + 0.11677701 * lai22 - 6.3942628e-3 * lai23 - 1.0984758 * llmi)
                           / (1.0 + 0.022930408 * lai2 + 7.5639466e-3 * lai22
                           - 4.146688e-4 * lai23 - 0.071648235 * llmi)) / 1.0e6
                           * zdzet * zdrox;
        zv = std::pow(10.0, dbzi / 10.0 - 18.0 + fzdri / 10.0);
        const double fkdpi = std::exp((54.728473 + 7.4097498 * lai2 - 7.5413134 * llmi
                           + 0.25580463 * lai22 + 0.25257362 * llmi2 - 0.50222586 * lai2 * llmi)
                           / (1.0 + 0.081702748 * lai2 - 0.085895343 * llmi
                           + 7.6712618e-4 * lai22 + 4.5412963e-4 * llmi2
                           - 9.3521894e-4 * lai2 * llmi)) / 2.0e18;
        const double kdrox = 1.75e-6 * std::exp(1.93990711 * lroi);
        const double kdzet = std::max(0.0, -0.0067191022 + 0.030489666 * lmdi
                           + 7.6470878 * zeta - 0.01097138 * lmdi2 + 9.8497711 * zeta2
                           + 0.83287327 * lmdi * zeta + 9.1590015e-4 * lmdi3
                           - 7.2917335 * zeta3 - 2.2235781 * lmdi * zeta2
                           - 0.053116593 * lmdi2 * zeta);
        kdp = ni * fkdpi * kdrox * kdzet * 0.11 / dbzwl;
        if (zv < zh) zv = 1.0001 * zh;
      } else if ((1.0 - adagr) >= 1.0e-2) {
        const double fdbzi = 390.47265 + 18.157805 * lai2 + 1.3558212 * lai22
                           - 0.074452587 * lai23 - 26.228751 * llmi;
        const double zhrox = 1.4172e-6 * std::exp(1.9687416987 * lroi);
        const double zhzet = 1.1531946 - 0.10474946 * lmdi + 1.3268439 * zeta
                           + 0.013483871 * lmdi2 + 3.0239757 * zeta2 - 0.53774565 * lmdi * zeta
                           - 7.6839476e-4 * lmdi3 + 4.3808336 * zeta3
                           + 0.080046674 * lmdi * zeta2 + 0.015810946 * lmdi2 * zeta;
        zh = ni * std::pow(10.0, (fdbzi - 200.0) / 10.0 - 18.0) * zhrox * zhzet;
        const double dbzi = 10.0 * std::log10(1.0e18 * zh);
        const double zdrox = 0.05683556 * lroi3 - 0.77585835 * lroi2
                           + 3.62207807 * lroi - 5.63826259;
        const double zdzet = std::max(0.0, 0.12567123 - 0.077409338 * lmdi
                           - 8.3245922 * zeta + 0.013440196 * lmdi2 + 8.7068519 * zeta2
                           - 0.53878869 * lmdi * zeta - 7.0609359e-4 * lmdi3
                           + 6.371616 * zeta3 - 1.6563927 * lmdi * zeta2
                           + 0.032682807 * lmdi2 * zeta);
        const double fzdri = std::exp((15.587817 + 0.23578544 * lai2
                           + 0.045497555 * lai22 - 1.0862071 * llmi + 3.4688369e-3 * llmi2)
                           / (1.0 + 0.015556043 * lai2 + 2.1779895e-3 * lai22
                           + 7.2733931e-5 * lai23 - 0.064401818 * llmi)) / 1.0e6
                           * zdzet * zdrox;
        zv = std::pow(10.0, dbzi / 10.0 - 18.0 - fzdri / 10.0);
        const double fkdpi = std::exp((55.096057 + 7.4047265 * lai2 - 7.6277374 * llmi
                           + 0.25324833 * lai22 + 0.25735457 * llmi2 - 0.50528828 * lai2 * llmi)
                           / (1.0 + 0.078240442 * lai2 - 0.083742581 * llmi
                           + 5.4080817e-4 * lai22 + 2.9610344e-4 * llmi2
                           - 5.9621769e-4 * lai2 * llmi)) / 2.0e18;
        const double kdrox = 1.89e-6 * std::exp(1.92963337 * lroi);
        const double kdzet = std::max(0.0, 0.044182106 - 0.020578987 * lmdi
                           - 8.8659688 * zeta + 0.0022174675 * lmdi2 + 5.2415098 * zeta2
                           - 0.40918383 * lmdi * zeta - 3.696556e-5 * lmdi3
                           + 3.5504491 * zeta3 - 1.4353113 * lmdi * zeta2
                           + 0.022796706 * lmdi2 * zeta);
        kdp = ni * fkdpi * kdrox * kdzet * 0.11 / dbzwl;
        if (zv > zh) zv = 0.9999 * zh;
      } else if (std::fabs(adagr - 1.0) < 1.0e-2) {
        const double fdbzi = 390.47265 + 18.157805 * lai2 + 1.3558212 * lai22
                           - 0.074452587 * lai23 - 26.228751 * llmi;
        const double zhrox = 1.9382e-6 * std::exp(1.9260732373 * lroi);
        const double zhzet = 1.0;
        zh = ni * std::pow(10.0, (fdbzi - 200.0) / 10.0 - 18.0) * zhrox * zhzet;
        const double dbzi = 10.0 * std::log10(1.0e18 * zh);
        zv = std::pow(10.0, dbzi / 10.0 - 18.0);
        kdp = 0.0;
      }
    }
  }
}

// ----------------------------------------------------------------------------
void dualpolOpSnowTcwa2(const int iband, const double tk, const double rho,
                        const double qc0, const double qr0, const double qs0,
                        const double ns0, const double smlf,
                        double &zh, double &zv, double &kdp) {
  const double dsmin = 2.0e-5, dsmax = 1.0e-2;
  zh = zxmin;
  zv = zxmin;
  kdp = 0.0;
  double qs = qs0 / 1.0e3;
  double qc = qc0 / 1.0e3;
  double qr = qr0 / 1.0e3;
  double ns = ns0 / rho;
  const double dbzwl = (iband == 1) ? 0.11 : 0.053;
  (void)qr;

  if (qs >= 1.0e-14 && ns >= 1.0e-2) {
    const double ltk = std::log(tk);
    const double ltk2 = ltk * ltk;
    const double lqs = -1.0 * std::log(qs);
    const double lqs2 = lqs * lqs;
    double rhos;
    if (tk < 273.15) {
      rhos = std::min(900.0, std::exp(75.401034 - 25.730216 * ltk + 0.81467519 * lqs
             + 2.339118 * ltk2 + 1.9976667e-3 * lqs2 - 0.14086434 * ltk * lqs));
    } else {
      rhos = std::max(100.0, std::exp(-64808.666 + 23113.508 * ltk - 36.46632 * lqs
             - 2060.6024 * ltk2 - 0.005729458 * lqs2 + 6.5057411 * ltk * lqs));
    }
    rhos = std::min(std::max(rhos, 100.0), 910.0);
    double mvrs = std::min(dsmax / 2.0,
                  std::max(dsmin / 2.0, std::pow(1.90986 * qs / ns / rhos, thrd)));
    double afas;
    if (qc >= 1.0e-14 || ns < 0.5) {
      afas = std::min(afamax, std::max(afamin,
             (afa1 + afa2 * std::log(ns) + afa3 * std::pow(std::log(ns), 2.0)
              + afa4 * std::log(mvrs))
             / (1.0 + afa5 * std::log(ns) + afa6 * std::log(mvrs))));
    } else {
      afas = std::min(afamax, std::max(afamin,
             (afb1 + afb2 * std::log(ns) + afb3 * std::pow(std::log(ns), 2.0)
              + afb4 * std::log(mvrs))
             / (1.0 + afb5 * std::log(ns) + afb6 * std::log(mvrs))));
    }
    if (tk > 273.15) afas = std::min(afamax, 0.01 * std::exp(-1.0 * std::log(mvrs)));
    double gs1 = tcwa2Gamln(afas + 1.0);
    double gs4 = tcwa2Gamln(afas + 4.0);
    double s3m = 1.90986 * qs / rhos;
    double lams = std::pow(std::exp(gs4 - gs1) * ns / s3m, thrd);
    double lamsmin = std::pow(std::exp(gs4 - gs1), thrd) / dsmax;
    double lamsmax = std::pow(std::exp(gs4 - gs1), thrd) / dsmin;
    if (lams < lamsmin) {
      lams = lamsmin;
      ns = s3m * std::exp(gs1 - gs4 + 3.0 * std::log(lams));
    } else if (lams > lamsmax) {
      lams = lamsmax;
      ns = s3m * std::exp(gs1 - gs4 + 3.0 * std::log(lams));
    }
    double mvds = std::pow(std::exp(gs4 - gs1), thrd) / lams;
    if (qs >= 1.0e-6) {
      const double llms = std::log(lams);
      const double lbds = std::log(mvds * 1.0e6);
      double lros = std::log(rhos);
      double saspr = std::max(0.2, std::min(1.0, 1.0207483 - 0.40684481 * lros
                   + 0.043489251 * std::pow(lros, 2.0) + 0.020462134 * lbds
                   + 0.014926001 * std::pow(lbds, 2.0) - 1.1252149e-3 * std::pow(lbds, 3.0)));
      if (tk < 238.16) saspr = 0.7;
      if (smlf >= 0.01 && tk >= 273.16) {
        const double tc = tk - 273.15;
        const double epsin = 5.27137 + 2.16474e-2 * tc - 1.31198e-3 * std::pow(tc, 2.0);
        const double epss = 78.54 * (1.0 - 4.579e-3 * (tc - 25.0)
                          + 1.19e-5 * std::pow(tc - 25.0, 2.0) - 2.8e-8 * std::pow(tc - 25.0, 3.0));
        const double alphe = -16.8129 / tk + 6.09265e-2;
        const double lamdz = 3.3836e-4 * std::exp(2513.98 / tk) * 1.0e-2;
        const double nerel = 1.0 + 2.0 * std::pow(lamdz / dbzwl, 1.0 - alphe)
                           * std::sin(1.570796327 * alphe)
                           + std::pow(lamdz / dbzwl, 2.0 - 2.0 * alphe);
        const double scawa = epsin + ((epss - epsin) * (std::pow(lamdz / dbzwl, 1.0 - alphe)
                           * std::sin(1.570796327 * alphe) + 1.0)) / nerel;
        const double scawb = ((epss - epsin) * (std::pow(lamdz / dbzwl, 1.0 - alphe)
                           * std::cos(1.570796327 * alphe))) / nerel + 1.25664 * dbzwl / 1.88496;
        const double scaw0 = std::pow(scawa, 2.0) + 4.0 * scawa + 4.0 + std::pow(scawb, 2.0);
        const double scakw =
            std::pow((std::pow(scawa, 2.0) + scawa - 2.0 + std::pow(scawb, 2.0)) / scaw0, 2.0)
            + std::pow(3.0 * scawb / scaw0, 2.0);
        s3m = 1.90986 * qs / rhos;
        mvrs = std::min(dsmax / 2.0,
               std::max(dsmin / 2.0, std::pow(1.90986 * qs / ns / rhos, thrd)));
        afas = std::min(afamax, std::max(afamin,
               std::exp(-18.420681 - 0.5 * std::log(ns) - 2.94 * std::log(mvrs))));
        gs1 = tcwa2Gamln(afas + 1.0);
        gs4 = tcwa2Gamln(afas + 4.0);
        lams = std::pow(std::exp(gs4 - gs1) * ns / s3m, thrd);
        lamsmin = std::pow(std::exp(gs4 - gs1), thrd) / dsmax;
        lamsmax = std::pow(std::exp(gs4 - gs1), thrd) / dsmin;
        if (lams < lamsmin) {
          lams = lamsmin;
          ns = s3m * std::exp(gs1 - gs4 + 3.0 * std::log(lams));
        } else if (lams > lamsmax) {
          lams = lamsmax;
          ns = s3m * std::exp(gs1 - gs4 + 3.0 * std::log(lams));
        }
        mvds = std::pow(std::exp(gs4 - gs1), thrd) / lams;
        const double smlv = std::pow(rhos / (1.0e3 / smlf - 1.0e3 + rhos), 0.3);
        const double dsmm = mvds * 1.0e3;
        saspr = smlf * std::max(0.4, 0.9951 + 2.51e-2 * dsmm - 3.644e-2 * std::pow(dsmm, 2.0)
              + 5.303e-3 * std::pow(dsmm, 3.0) - 2.492e-4 * std::pow(dsmm, 4.0))
              + (1.0 - smlf) * saspr;
        rhos = std::min(910.0, std::max(50.0, smlf * 1.0e3 + (1.0 - smlf) * rhos));
        const double ora0 = 10.0 * smlf + (1.0 - smlf) * 40.0;
        const double orar = std::exp(-2.0 * std::pow(ora0 * 0.017453293, 2.0));
        const double ora1 = (3.0 + 4.0 * orar + std::pow(orar, 4.0)) / 8.0;
        const double ora2 = (3.0 - 4.0 * orar + std::pow(orar, 4.0)) / 8.0;
        const double ora3 = (1.0 - std::pow(orar, 4.0)) / 8.0;
        const double sca0 = std::max(1.0e-6, std::sqrt(std::pow(saspr, -2.0) - 1.0));
        const double scala = (1.0 + std::pow(sca0, 2.0)) / std::pow(sca0, 2.0)
                           * (1.0 - std::atan(sca0) / sca0);
        const double scalb = (1.0 - scala) / 2.0;
        const double scaw1 = -1.0892332e-2 + 2.0820143e-2 / smlv + 4.8553456e-6 * rhos
                           - 7.6529057e-5 / std::pow(smlv, 2.0) - 1.1783334e-9 * std::pow(rhos, 2.0)
                           - 2.9599665e-6 * rhos / smlv;
        const double scafa = 1.0 / (scala + scaw1);
        const double scafb = 1.0 / (scalb + scaw1);
        const double rads = ns * std::exp(tcwa2Gamln(afas + 7.0) - gs1 - 6.0 * std::log(lams));
        zh = std::max(1.0e-21, rads * (ora1 * std::pow(scafb, 2.0)
           + 2.0 * ora3 * scafa * scafb + ora2 * std::pow(scafa, 2.0)) / (9.0 * scakw));
        zh = zh * (std::max(0.0, -0.1 * afas + 0.5) * smlf + 1.0);
        zv = std::max(1.0e-21, rads * (ora1 * std::pow(scafa, 2.0)
           + 2.0 * ora3 * scafa * scafb + ora2 * std::pow(scafb, 2.0)) / (9.0 * scakw));
        kdp = 94247.77961 / dbzwl * ns * std::exp(gs4 - gs1 - 3.0 * std::log(lams))
            * std::fabs(scafa - scafb) * orar;
        if (zv > zh) zv = 0.9999 * zh;
      } else {
        const double las2 = std::log(afas + 2.0);
        lros = std::log(rhos);
        const double lasp = std::log(200.0 * std::max(0.01, saspr));
        const double fdbzs = 428.62235 + 18.226687 * las2 + 1.299724 * std::pow(las2, 2.0)
                           - 0.070558679 * std::pow(las2, 3.0) - 26.070123 * llms;
        const double zhrox = 1.544e-5 * std::exp(1.95714566 * lros);
        zh = ns * std::pow(10.0, (fdbzs - 250.0) / 10.0 - 18.0) * zhrox;
        const double dbzs = 10.0 * std::log10(1.0e18 * zh);
        const double zdrox = 0.14542522 * std::pow(lros, 3.0) - 1.97576772 * std::pow(lros, 2.0)
                           + 9.33469066 * lros - 14.95051206;
        const double zdasp = std::max(0.0, -8.85535e-3 * std::pow(lasp, 3.0)
                           - 0.10877652 * std::pow(lasp, 2.0) + 0.24567081 * lasp + 3.06007047);
        const double fzdrs = 0.384 * zdrox * zdasp / 0.37728;
        zv = std::pow(10.0, dbzs / 10.0 - 18.0 - fzdrs / 10.0);
        const double fkdps = std::exp(44.30884 + 2.9903578 * las2 - 0.25943237 / las2
                           - 3.0036608 * llms) / 1.0e15;
        const double kdrox = 1.974e-5 * std::exp(1.90820312 * lros);
        const double kdasp = std::max(0.0, -0.00958782 * std::pow(lasp, 3.0)
                           - 0.10043861 * std::pow(lasp, 2.0) + 0.22095032 * lasp + 3.06677159);
        kdp = ns * fkdps * kdrox * kdasp * 0.11 / dbzwl;
        if (zv > zh) zv = 0.9999 * zh;
      }
    }
  }
}

// ----------------------------------------------------------------------------
void dualpolOpGraupTcwa2(const int iband, const double tk, const double rho,
                         const double qc0, const double qr0, const double qg0,
                         const double ng0, const double gmlf,
                         double &zh, double &zv, double &kdp) {
  const double dgmin = 5.0e-5, dgmax = 2.0e-2;
  zh = zxmin;
  zv = zxmin;
  kdp = 0.0;
  double qg = qg0 / 1.0e3;
  double qc = qc0 / 1.0e3;
  double qr = qr0 / 1.0e3;
  double ng = ng0 / rho;
  const double dbzwl = (iband == 1) ? 0.11 : 0.053;

  if (qg >= 1.0e-14 && ng >= 1.0e-2) {
    int ihail = 0;
    double rhog = 400.0;
    const double tc = tk - 273.15;
    const double mvg0 = std::exp(8.8283 + 0.25 * std::log(qg)) * 1.0e-6;
    const double dsll = std::max(1.0e-3, std::min(1.0,
                        2.0 * 1.0e-2 * (std::exp(std::min(20.0,
                        -tc / std::max(0.1, (1.1e4 * (qc + qr) - 1.3e3 * qg + 1.0)))) - 1.0)));
    if (mvg0 > dsll) {
      ihail = 1;
      rhog = 900.0;
    } else {
      ihail = 0;
      const double lqg = -1.0 * std::log(qg);
      rhog = std::min(900.0, std::max(400.0, std::exp(-9.57e-5 * std::pow(lqg, 3.0)
             + 3.077e-3 * std::pow(lqg, 2.0) - 6.800923e-2 * lqg + 6.8175231)));
    }
    double mvrg = std::min(dgmax / 2.0,
                  std::max(dgmin / 2.0, std::pow(1.90986 * qg / ng / rhog, thrd)));
    double afag;
    if (qc >= 1.0e-14 || ng < 1.0) {
      afag = std::min(afamax, std::max(afamin,
             (afa1 + afa2 * std::log(ng) + afa3 * std::pow(std::log(ng), 2.0)
              + afa4 * std::log(mvrg))
             / (1.0 + afa5 * std::log(ng) + afa6 * std::log(mvrg))));
    } else {
      afag = std::min(afamax, std::max(afamin,
             (afb1 + afb2 * std::log(ng) + afb3 * std::pow(std::log(ng), 2.0)
              + afb4 * std::log(mvrg))
             / (1.0 + afb5 * std::log(ng) + afb6 * std::log(mvrg))));
    }
    if (tk > 273.15) afag = std::min(afamax, 0.01 * std::exp(-1.0 * std::log(mvrg)));
    double gg1 = tcwa2Gamln(afag + 1.0);
    double gg4 = tcwa2Gamln(afag + 4.0);
    double g3m = 1.90986 * qg / rhog;
    double lamg = std::pow(std::exp(gg4 - gg1) * ng / g3m, thrd);
    double lamgmin = std::pow(std::exp(gg4 - gg1), thrd) / dgmax;
    double lamgmax = std::pow(std::exp(gg4 - gg1), thrd) / dgmin;
    if (lamg < lamgmin) {
      lamg = lamgmin;
      ng = g3m * std::exp(gg1 - gg4 + 3.0 * std::log(lamg));
    } else if (lamg > lamgmax) {
      lamg = lamgmax;
      ng = g3m * std::exp(gg1 - gg4 + 3.0 * std::log(lamg));
    }
    double mvdg = std::pow(std::exp(gg4 - gg1), thrd) / lamg;
    if (qg >= 1.0e-6) {
      const double llmg = std::log(lamg);
      const double llmg2 = llmg * llmg;
      double gaspr = 0.0;
      if (ihail == 0) {
        const double lbdg = std::log(mvdg * 1.0e6);
        const double lrog = std::log(rhog);
        gaspr = std::max(0.2, std::min(1.0, 1.0207483 - 0.40684481 * lrog
              + 0.043489251 * std::pow(lrog, 2.0) + 0.020462134 * lbdg
              + 0.014926001 * std::pow(lbdg, 2.0) - 1.1252149e-3 * std::pow(lbdg, 3.0)));
      } else if (ihail == 1) {
        const double lag5 = std::log(afag + 5.0);
        const double lag52 = lag5 * lag5;
        gaspr = std::min(1.0, ((5.6387098 + 1.5626893 * lag5 - 1.5481826 * llmg
              + 0.12481435 * lag52 + 0.11244742 * llmg2 - 0.2336531 * lag5 * llmg)
              / (1.0 + 0.28272803 * lag5 - 0.2780308 * llmg + 2.1750225e-2 * lag52
              + 1.9960971e-2 * llmg2 - 4.1317594e-2 * lag5 * llmg)) / 10.0);
      }
      if (gmlf >= 0.01 && tk >= 273.16) {
        const double epsin = 5.27137 + 2.16474e-2 * tc - 1.31198e-3 * std::pow(tc, 2.0);
        const double epss = 78.54 * (1.0 - 4.579e-3 * (tc - 25.0)
                          + 1.19e-5 * std::pow(tc - 25.0, 2.0) - 2.8e-8 * std::pow(tc - 25.0, 3.0));
        const double alphe = -16.8129 / tk + 6.09265e-2;
        const double lamdz = 3.3836e-4 * std::exp(2513.98 / tk) * 1.0e-2;
        const double nerel = 1.0 + 2.0 * std::pow(lamdz / dbzwl, 1.0 - alphe)
                           * std::sin(1.570796327 * alphe)
                           + std::pow(lamdz / dbzwl, 2.0 - 2.0 * alphe);
        const double scawa = epsin + ((epss - epsin) * (std::pow(lamdz / dbzwl, 1.0 - alphe)
                           * std::sin(1.570796327 * alphe) + 1.0)) / nerel;
        const double scawb = ((epss - epsin) * (std::pow(lamdz / dbzwl, 1.0 - alphe)
                           * std::cos(1.570796327 * alphe))) / nerel + 1.25664 * dbzwl / 1.88496;
        const double scaw0 = std::pow(scawa, 2.0) + 4.0 * scawa + 4.0 + std::pow(scawb, 2.0);
        const double scakw =
            std::pow((std::pow(scawa, 2.0) + scawa - 2.0 + std::pow(scawb, 2.0)) / scaw0, 2.0)
            + std::pow(3.0 * scawb / scaw0, 2.0);
        g3m = 1.90986 * qg / rhog;
        mvrg = std::min(dgmax / 2.0,
               std::max(dgmin / 2.0, std::pow(1.90986 * qg / ng / rhog, thrd)));
        afag = std::min(afamax, std::max(afamin,
               std::exp(-18.420681 - 0.5 * std::log(ng) - 2.94 * std::log(mvrg))));
        gg1 = tcwa2Gamln(afag + 1.0);
        gg4 = tcwa2Gamln(afag + 4.0);
        lamg = std::pow(std::exp(gg4 - gg1) * ng / g3m, thrd);
        lamgmin = std::pow(std::exp(gg4 - gg1), thrd) / dgmax;
        lamgmax = std::pow(std::exp(gg4 - gg1), thrd) / dgmin;
        if (lamg < lamgmin) {
          lamg = lamgmin;
          ng = g3m * std::exp(gg1 - gg4 + 3.0 * std::log(lamg));
        } else if (lamg > lamgmax) {
          lamg = lamgmax;
          ng = g3m * std::exp(gg1 - gg4 + 3.0 * std::log(lamg));
        }
        mvdg = std::pow(std::exp(gg4 - gg1), thrd) / lamg;
        const double gmlv = std::pow(rhog / (1.0e3 / gmlf - 1.0e3 + rhog), 0.3);
        const double dgmm = mvdg * 1.0e3;
        gaspr = gmlf * std::max(0.4, 0.9951 + 2.51e-2 * dgmm - 3.644e-2 * std::pow(dgmm, 2.0)
              + 5.303e-3 * std::pow(dgmm, 3.0) - 2.492e-4 * std::pow(dgmm, 4.0))
              + (1.0 - gmlf) * gaspr;
        rhog = std::min(900.0, std::max(100.0, gmlf * 1.0e3 + (1.0 - gmlf) * rhog));
        const double ora0 = 10.0 * gmlf + (1.0 - gmlf) * 40.0;
        const double orar = std::exp(-2.0 * std::pow(ora0 * 0.017453293, 2.0));
        const double ora1 = (3.0 + 4.0 * orar + std::pow(orar, 4.0)) / 8.0;
        const double ora2 = (3.0 - 4.0 * orar + std::pow(orar, 4.0)) / 8.0;
        const double ora3 = (1.0 - std::pow(orar, 4.0)) / 8.0;
        const double sca0 = std::max(1.0e-6, std::sqrt(std::pow(gaspr, -2.0) - 1.0));
        const double scala = (1.0 + std::pow(sca0, 2.0)) / std::pow(sca0, 2.0)
                           * (1.0 - std::atan(sca0) / sca0);
        const double scalb = (1.0 - scala) / 2.0;
        const double scaw1 = -1.0892332e-2 + 2.0820143e-2 / gmlv + 4.8553456e-6 * rhog
                           - 7.6529057e-5 / std::pow(gmlv, 2.0) - 1.1783334e-9 * std::pow(rhog, 2.0)
                           - 2.9599665e-6 * rhog / gmlv;
        const double scafa = 1.0 / (scala + scaw1);
        const double scafb = 1.0 / (scalb + scaw1);
        const double radg = ng * std::exp(tcwa2Gamln(afag + 7.0) - gg1 - 6.0 * std::log(lamg));
        zh = std::max(1.0e-21, radg * (ora1 * std::pow(scafb, 2.0)
           + 2.0 * ora3 * scafa * scafb + ora2 * std::pow(scafa, 2.0)) / (9.0 * scakw));
        zh = zh * (std::max(0.0, -0.1 * afag + 0.5) * gmlf + 1.0);
        zv = std::max(1.0e-21, radg * (ora1 * std::pow(scafa, 2.0)
           + 2.0 * ora3 * scafa * scafb + ora2 * std::pow(scafb, 2.0)) / (9.0 * scakw));
        kdp = 94247.77961 / dbzwl * ng * std::exp(gg4 - gg1 - 3.0 * std::log(lamg))
            * std::fabs(scafa - scafb) * orar;
        if (zv > zh) zv = 0.9999 * zh;
      } else {
        const double lag2 = std::log(afag + 2.0);
        const double lrog = std::log(rhog);
        const double lasp = std::log(200.0 * std::max(0.01, gaspr));
        const double zhrox = 1.544e-5 * std::exp(1.95714566 * lrog);
        const double fdbzg = 428.62235 + 18.226687 * lag2 + 1.299724 * std::pow(lag2, 2.0)
                           - 0.070558679 * std::pow(lag2, 3.0) - 26.070123 * llmg;
        zh = ng * std::pow(10.0, (fdbzg - 250.0) / 10.0 - 18.0) * zhrox;
        const double dbzg = 10.0 * std::log10(1.0e18 * zh);
        const double zdrox = 0.14542522 * std::pow(lrog, 3.0) - 1.97576772 * std::pow(lrog, 2.0)
                           + 9.33469066 * lrog - 14.95051206;
        const double zdasp = std::max(0.0, -8.85535e-3 * std::pow(lasp, 3.0)
                           - 0.10877652 * std::pow(lasp, 2.0) + 0.24567081 * lasp + 3.06007047);
        const double fzdrg = 0.384 * zdrox * zdasp;
        zv = std::pow(10.0, dbzg / 10.0 - 18.0 - fzdrg / 10.0);
        const double fkdpg = std::exp(44.30884 + 2.9903578 * lag2 - 0.25943237 / lag2
                           - 3.0036608 * llmg) / 1.0e15;
        const double kdrox = 1.974e-5 * std::exp(1.90820312 * lrog);
        const double kdasp = std::max(0.0, -0.00958782 * std::pow(lasp, 3.0)
                           - 0.10043861 * std::pow(lasp, 2.0) + 0.22095032 * lasp + 3.06677159);
        kdp = ng * fkdpg * kdrox * kdasp * 0.11 / dbzwl;
        if (zv > zh) zv = 0.9999 * zh;
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Aggregation wrapper (mirrors tcwa2_compute_point in dualpol_op_mod.f90).
void tcwa2ComputePoint(const int iband, const double density_air,
                       const double temp_air, const double qr, const double qs,
                       const double qg, const double qc, const double qi,
                       const double nr, const double ns, const double ng,
                       const double ni, const double smlf, const double gmlf,
                       double &zh, double &zdr, double &kdp, double &phv) {
  double zh_rain, zv_rain, kdp_rain;
  double zh_ice, zv_ice, kdp_ice;
  double zh_snow, zv_snow, kdp_snow;
  double zh_graup, zv_graup, kdp_graup;

  dualpolOpRainTcwa2(iband, density_air, qc, qr, nr, zh_rain, zv_rain, kdp_rain);
  dualpolOpIceTcwa2(iband, temp_air, density_air, qc, qi, ni, zh_ice, zv_ice, kdp_ice);
  dualpolOpSnowTcwa2(iband, temp_air, density_air, qc, qr, qs, ns, smlf,
                     zh_snow, zv_snow, kdp_snow);
  dualpolOpGraupTcwa2(iband, temp_air, density_air, qc, qr, qg, ng, gmlf,
                      zh_graup, zv_graup, kdp_graup);

  zh = zh_rain + zh_ice + zh_snow + zh_graup;
  kdp = kdp_rain + kdp_ice + kdp_snow + kdp_graup;

  // ZDR = Zh / Zv (linear).
  const double zv_sum = zv_rain + zv_ice + zv_snow + zv_graup;
  zdr = zh / std::max(zv_sum, 1.0e-20);

  // TCWA2 does not provide phv explicitly.
  phv = 0.95;

  // Convert zh from m^6/m^3 (SI) to mm^6/m^3 to match the Zhang21 convention.
  zh = zh * 1.0e18;
}

}  // namespace ppro
}  // namespace ufo
