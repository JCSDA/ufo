/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_SYMMCLDIMPACTUTILS_H_
#define UFO_UTILS_SYMMCLDIMPACTUTILS_H_

#include <cstddef>
#include <vector>

namespace ufo {

/// \brief Core computation of Symmetric Cloud Impact using Harnish et al. (2016)
///
/// Operates on raw vectors so it can be called from both obsfunction and predictor.
/// Result is written into out[jloc * nvars + jvar].
/// All brightness temperatures must be raw (pre-bias-correction) values.
///
/// \param bak    Raw model brightness temperature [nlocs]
/// \param obs    Observed brightness temperature [nlocs]
/// \param btlim  Per-channel BT limit threshold for this channel
/// \param nlocs  Number of locations
/// \param jvar   Channel index
/// \param nvars  Number of channels/variables
/// \param order  Power to raise SCI to
/// \param out    Output vector [nlocs * nvars]
///
/// Harnisch, F., M. Weissmann, and A. Perianez, 2016: Error model for the
/// assimilation of cloud-affected infrared satellite observations in an
/// ensemble data assimilation system. Q.J.R. Meteorol. Soc., 142, 1797-1808.
/// doi:10.1002/qj.2776
void computeSCIHarnish(const std::vector<float> & bak,
                       const std::vector<float> & obs,
                       float btlim,
                       size_t nlocs,
                       size_t jvar,
                       size_t nvars,
                       int order,
                       std::vector<float> & out);

/// \brief Core computation of Symmetric Cloud Impact using Okamoto et al. (2014)
///
/// Operates on raw vectors so it can be called from both obsfunction and predictor.
/// Result is written into out[jloc * nvars + jvar].
/// All brightness temperatures must be raw (pre-bias-correction) values.
///
/// \param clr          Raw CRTM clear-sky model brightness temperature [nlocs]
/// \param bak          Raw model brightness temperature [nlocs]
/// \param obs          Observed brightness temperature [nlocs]
/// \param scale_by_omb If true, scale SCI by sigmoid function of OmB
/// \param sigmoid_c1   Slope of sigmoid function
/// \param sigmoid_c2   50-percent value of sigmoid function
/// \param nlocs        Number of locations
/// \param jvar         Channel index
/// \param nvars        Number of channels/variables
/// \param order        Power to raise SCI to
/// \param out          Output vector [nlocs * nvars]
///
/// Okamoto, K., McNally, A.P. and Bell, W. (2014), Progress towards the
/// assimilation of all-sky infrared radiances: an evaluation of cloud
/// effects. Q.J.R. Meteorol. Soc., 140: 1603-1614. doi:10.1002/qj.2242
void computeSCIOkamoto(const std::vector<float> & clr,
                       const std::vector<float> & bak,
                       const std::vector<float> & obs,
                       bool scale_by_omb,
                       float sigmoid_c1,
                       float sigmoid_c2,
                       size_t nlocs,
                       size_t jvar,
                       size_t nvars,
                       int order,
                       std::vector<float> & out);

}  // namespace ufo

#endif  // UFO_UTILS_SYMMCLDIMPACTUTILS_H_
