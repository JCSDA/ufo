/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_SAVEBIASCOEFFS_H_
#define UFO_UTILS_SAVEBIASCOEFFS_H_

#include <string>
#include <vector>

#include "ioda/ObsGroup.h"

namespace ufo {

/// Create ObsGroup with dimensions npredictors = size(predictors) and
/// nchannels = size(channels), variables predictors, channels and
/// bias_cooefficients (npredictors x nchannels)
ioda::ObsGroup saveBiasCoeffsWithChannels(ioda::Group &,
                                          const std::vector<std::string> &,
                                          const std::vector<int> &,
                                          const Eigen::MatrixXd &);

}  // namespace ufo

#endif  // UFO_UTILS_SAVEBIASCOEFFS_H_
