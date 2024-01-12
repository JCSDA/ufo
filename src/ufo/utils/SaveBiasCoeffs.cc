/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/utils/SaveBiasCoeffs.h"

#include "ioda/ObsGroup.h"
#include "oops/util/missingValues.h"

namespace ufo {

ioda::ObsGroup saveBiasCoeffsWithChannels(ioda::Group & parent,
                                          const std::vector<std::string> & predictors,
                                          const std::vector<int> & channels,
                                          const Eigen::MatrixXd & coeffs) {
  // dimensions
  ioda::NewDimensionScales_t dims {
      ioda::NewDimensionScale<int>("npredictors", predictors.size()),
      ioda::NewDimensionScale<int>("nchannels", channels.size())
  };
  // new ObsGroup
  ioda::ObsGroup ogrp = ioda::ObsGroup::generate(parent, dims);

  // save the predictors
  ioda::Variable predsVar = ogrp.vars.createWithScales<std::string>(
                            "predictors", {ogrp.vars["npredictors"]});
  predsVar.write(predictors);
  // and the variables
  ioda::Variable chansVar = ogrp.vars.createWithScales<int>("channels", {ogrp.vars["nchannels"]});
  chansVar.write(channels);

  // Set up the creation parameters for the bias coefficients variable
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;               // allow chunking
  float_params.compressWithGZIP();         // compress using gzip
  const float missing_value = util::missingValue<float>();
  float_params.setFillValue<float>(missing_value);

  // Create a variable for bias coefficients, save bias coeffs to the variable
  ioda::Variable biasVar = ogrp.vars.createWithScales<float>("bias_coefficients",
                           {ogrp.vars["npredictors"], ogrp.vars["nchannels"]}, float_params);
  biasVar.writeWithEigenRegular(coeffs);
  return ogrp;
}

}  // namespace ufo
