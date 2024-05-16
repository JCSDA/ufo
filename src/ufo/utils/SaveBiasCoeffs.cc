/*
 * (C) Copyright 2024 UCAR
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
  ioda::NewDimensionScales_t newDims {
      ioda::NewDimensionScale<int>("Record", 1),
      ioda::NewDimensionScale<int>("Channel", channels.size())
  };

  // new ObsGroup
  ioda::ObsGroup ogrp = ioda::ObsGroup::generate(parent, newDims);

  // Set up the creation parameters for the bias coefficients variable
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;               // allow chunking
  float_params.compressWithGZIP();         // compress using gzip
  const float missing_value = util::missingValue<float>();
  float_params.setFillValue<float>(missing_value);

  // Loop over predictors and create variables
  for (size_t jpred = 0; jpred < predictors.size(); ++jpred) {
    // create and write the bias coeffs
    ioda::Variable biasVar = ogrp.vars.createWithScales<float>(
                             "BiasCoefficients/"+predictors[jpred],
                             {ogrp.vars["Record"], ogrp.vars["Channel"]}, float_params);
    biasVar.writeWithEigenRegular(coeffs(jpred, Eigen::all));
  }

  // Save the Channel
  // and the variables
  ioda::Variable chansVar = ogrp.vars.open("Channel");
  chansVar.write(channels);

  return ogrp;
}

}  // namespace ufo
