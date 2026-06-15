/*
 * (C) Copyright 2024-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_TCWA2FORWARD_H_
#define UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_TCWA2FORWARD_H_

// ----------------------------------------------------------------------------
// TCWA2 polarimetric radar operator (C++ port of tcwa2_forward_mod.f90 and the
// aggregation wrapper from dualpol_op_mod.f90).
//
// Empirically fitted analytic parameterizations assuming a gamma PSD, derived
// from offline bin scattering integrated under the Rayleigh approximation.
//
// Author: Rong Kong (NCAR/MMM) - C++ rewrite of the operator from Fortran;
//         integration into the PPRO library (2025).
//
// Original Fortran author: Tzu-Chin Tsai (CWA, original implementation).
//
// Reference: Tsai, T.-C., J.-P. Chen, Z. Liu, S.-Y. Jiang, R. Kong, et al.,
//            2026 (J. Adv. Model. Earth Syst., submitted).
// ----------------------------------------------------------------------------

namespace ufo {
namespace ppro {

/// ln(Gamma(xx)) for xx > 0 (wraps the C++ standard-library std::lgamma).
double tcwa2Gamln(double xx);

/// Radar-band convention used throughout the TCWA2 routines below:
///   iband = 1 -> S-band (wavelength 11 cm), iband = 2 -> C-band (5.3 cm).

/// Individual hydrometeor contributions.  Outputs zh, zv [mm^6/m^3] and
/// kdp [deg/km].
void dualpolOpRainTcwa2(const int iband, const double rho, const double qc0,
                        const double qr0, const double nr0,
                        double &zh, double &zv, double &kdp);
void dualpolOpIceTcwa2(const int iband, const double tk, const double rho,
                       const double qc0, const double qi0, const double ni0,
                       double &zh, double &zv, double &kdp);
void dualpolOpSnowTcwa2(const int iband, const double tk, const double rho,
                        const double qc0, const double qr0, const double qs0,
                        const double ns0, const double smlf,
                        double &zh, double &zv, double &kdp);
void dualpolOpGraupTcwa2(const int iband, const double tk, const double rho,
                         const double qc0, const double qr0, const double qg0,
                         const double ng0, const double gmlf,
                         double &zh, double &zv, double &kdp);

/// Aggregation wrapper: sums the per-species contributions and converts to the
/// Zhang21 output convention (zh in mm^6/m^3, zdr linear, phv placeholder).
///   density_air [kg/m^3], temp_air [K], q* [g/kg], n* [1/m^3].
void tcwa2ComputePoint(const int iband, const double density_air,
                       const double temp_air, const double qr, const double qs,
                       const double qg, const double qc, const double qi,
                       const double nr, const double ns, const double ng,
                       const double ni, const double smlf, const double gmlf,
                       double &zh, double &zdr, double &kdp, double &phv);

}  // namespace ppro
}  // namespace ufo

#endif  // UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_TCWA2FORWARD_H_
