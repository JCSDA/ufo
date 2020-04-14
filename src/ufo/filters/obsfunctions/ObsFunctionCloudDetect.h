/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLOUDDETECT_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLOUDDETECT_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Cloud Detection Algorithm (Minimum Residual Method) for Infrared sensors
// using selected channels from 15 microns CO2 absorption band
// Output of this function:
// 0 = channel is not affected by clouds (clear channel)
// 1 = channel is affected by clouds (cloudy channel)
// 2 = channel is not affected by clouds but too sensitive to surface condition
// -----------------------------------------------------------------------------

class ObsFunctionCloudDetect : public ObsFunctionBase {
 public:
  explicit ObsFunctionCloudDetect(const eckit::LocalConfiguration);
  ~ObsFunctionCloudDetect();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::string biasgrp_;
  std::string errgrp_;
  std::string hofxgrp_;
  std::vector<int> channels_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLOUDDETECT_H_
