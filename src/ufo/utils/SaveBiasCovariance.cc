/*
 * (C) Crown Copyright 2025, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/utils/SaveBiasCovariance.h"

#include "ioda/ObsGroup.h"
#include "oops/util/missingValues.h"

namespace ufo {

void saveBiasCovarianceWithChannels(ioda::Group & parent,
                                    const std::vector<std::string> & predictors,
                                    const std::vector<int> & channels,
                                    const std::vector<int> & obs_assimilated,
                                    const Eigen::MatrixXd & allbcerrors) {
  // dimensions
  ioda::NewDimensionScales_t dims {
        ioda::NewDimensionScale<int>("Record", 1),
        ioda::NewDimensionScale<int>("Channel", channels.size())
  };
  // new ObsGroup
  ioda::ObsGroup ogrp = ioda::ObsGroup::generate(parent, dims);

  // Set up the creation parameters for the bias covariance coefficients variable
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;               // allow chunking
  float_params.compressWithGZIP();         // compress using gzip
  const float missing_value = util::missingValue<float>();
  float_params.setFillValue<float>(missing_value);

  // Loop over predictors and create variables
  for (size_t jpred = 0; jpred < predictors.size(); ++jpred) {
    // create and write the bias covariance coeffs
    ioda::Variable anvarVar = ogrp.vars.createWithScales<float>(
                             "BiasCoefficientErrors/"+predictors[jpred],
                             {ogrp.vars["Record"], ogrp.vars["Channel"]}, float_params);
    anvarVar.writeWithEigenRegular(allbcerrors(jpred, Eigen::all));
  }

  // write number_obs_assimilated
  ioda::Variable nobsVar = ogrp.vars.createWithScales<int>(
                           "numberObservationsUsed", {ogrp.vars["Record"], ogrp.vars["Channel"]});
  nobsVar.write(obs_assimilated);

  // write out the channels
  ioda::Variable chansVar = ogrp.vars.open("Channel");
  chansVar.write(channels);
}

void saveBiasCovarianceWithRecords(ioda::Group & parent,
                                   const std::vector<std::string> & predictors,
                                   const std::vector<std::string> & records,
                                   const std::vector<std::string> & vars,
                                   const std::vector<int> & obs_assimilated,
                                   const Eigen::MatrixXd & allbcerrors) {
    // dimensions
    ioda::NewDimensionScales_t dims {
        ioda::NewDimensionScale<int>("Record", records.size()),
        ioda::NewDimensionScale<int>("Variable", vars.size())
    };

    // new ObsGroup
    ioda::ObsGroup ogrp = ioda::ObsGroup::generate(parent, dims);

    // Set up the creation parameters for the bias covariance coefficients variable
    ioda::VariableCreationParameters float_params;
    float_params.chunk = true;               // allow chunking
    float_params.compressWithGZIP();         // compress using gzip
    const float missing_value = util::missingValue<float>();
    float_params.setFillValue<float>(missing_value);

    // Loop over predictors and create variables
    for (size_t jpred = 0; jpred < predictors.size(); ++jpred) {
      // create and write the bias coeffs
      ioda::Variable biasVar = ogrp.vars.createWithScales<float>(
                               "BiasCoefficientErrors/" + predictors[jpred],
                               {ogrp.vars["Record"], ogrp.vars["Variable"]}, float_params);
      biasVar.writeWithEigenRegular(allbcerrors(jpred, Eigen::all));
    }

    // write number_obs_assimilated
    ioda::Variable nobsVar = ogrp.vars.createWithScales<int>("numberObservationsUsed",
                                            {ogrp.vars["Record"], ogrp.vars["Variable"]});
    nobsVar.write(obs_assimilated);

    // Save the stationIdentification
    // and the variables
    ioda::Variable variablesVar = ogrp.vars.createWithScales<std::string>(
                "Variables", {ogrp.vars["Variable"]});
    variablesVar.write(vars);
    ioda::Variable recordVar = ogrp.vars.createWithScales<std::string>(
                "stationIdentification", {ogrp.vars["Record"]});
    recordVar.write(records);
}

}  // namespace ufo
