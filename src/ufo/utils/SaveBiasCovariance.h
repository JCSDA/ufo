/*
 * (C) Crown Copyright 2025, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_SAVEBIASCOVARIANCE_H_
#define UFO_UTILS_SAVEBIASCOVARIANCE_H_

#include <string>
#include <vector>

#include "ioda/ObsGroup.h"

namespace ufo {

void saveBiasCovarianceWithChannels(ioda::Group &,
                              const std::vector<std::string> &,
                              const std::vector<int> &,
                              const std::vector<int> &,
                              const Eigen::MatrixXd &);

void saveBiasCovarianceWithRecords(ioda::Group &,
                             const std::vector<std::string> &,
                             const std::vector<std::string> &,
                             const std::vector<std::string> &,
                             const std::vector<int> &,
                             const Eigen::MatrixXd &);

}  // namespace ufo

#endif  // UFO_UTILS_SAVEBIASCOVARIANCE_H_
