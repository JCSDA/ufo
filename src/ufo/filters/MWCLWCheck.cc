/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/MWCLWCheck.h"

#include <math.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

MWCLWCheck::MWCLWCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters) {
  oops::Log::debug() << "MWCLWCheck: config = " << parameters_ << std::endl;

  const oops::ObsVariables &invars = parameters_.clwVariables.value();
  ASSERT(invars.size() == 2);
  const Variable var0("HofX/" + invars[0]);
  const Variable var1("HofX/" + invars[1]);
  allvars_ += var0;
  allvars_ += var1;
}

// -----------------------------------------------------------------------------

MWCLWCheck::~MWCLWCheck() {}

// -----------------------------------------------------------------------------

void MWCLWCheck::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "MWCLWCheck postFilter" << std::endl;

  const oops::ObsVariables observed = obsdb_.assimvariables();
  const float missing = util::missingValue<float>();
  float amsua_clw(float, float, float);

// Get config
  const std::vector<float> &clw_thresholds = parameters_.clwThresholds;
// clw_option controls how the clw is calculated:
//     1) Use observed BTs
//     2) Use calculated BTs
//     3) Symmetric calculation
  const int clw_option = parameters_.clwOption;
  const oops::ObsVariables &invars = parameters_.clwVariables.value();
  oops::Log::debug() << "MWCLWCheck: config = " << parameters_ << std::endl;

// Number of channels to be tested and number of thresholds needs to be the same
  ASSERT(clw_thresholds.size() == filtervars.nvars());
// Thresholds need to have non-missing values
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    ASSERT(clw_thresholds[jv] != missing);
  }

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");
  ioda::ObsDataVector<float> obs_for_calc(obsdb_, invars, "ObsValue");
  ioda::ObsDataVector<float> sza(obsdb_, "sensorZenithAngle", "MetaData");
  ioda::ObsDataVector<float> clw(obsdb_, "cloud_liquid_water", "Diagnostic", false);
  ioda::ObsDataVector<float> clw_guess_out(obsdb_, "clws_guess", "Diagnostic", false);
  ioda::ObsDataVector<float> clw_obs_out(obsdb_, "clw_obs", "Diagnostic", false);

//  H(x)
  const Variable var0("HofX/" + invars[0]);
  const Variable var1("HofX/" + invars[1]);
  std::vector<float> hofx0;
  data_.get(var0, hofx0);
  std::vector<float> hofx1;
  data_.get(var1, hofx1);

  auto clw_obs = [&](size_t jobs) { return amsua_clw(
    obs_for_calc[0][jobs], obs_for_calc[1][jobs], sza[0][jobs]); };

  auto clw_guess = [&](size_t jobs) { return amsua_clw(
    hofx0[jobs], hofx1[jobs], sza[0][jobs]); };

// Loop over obs locations calculating CLW from observations
  for (size_t jobs = 0; jobs < obs.nlocs(); ++jobs) {
    switch (clw_option) {
      case 1 :
        clw[0][jobs] = clw_obs(jobs);
        break;
      case 2 :
        clw[0][jobs] = clw_guess(jobs);
        break;
      case 3 :
        clw_obs_out[0][jobs] = clw_obs(jobs);
        clw_guess_out[0][jobs] = clw_guess(jobs);
        if (clw_obs_out[0][jobs] == missing || clw_guess_out[0][jobs] == missing) {
          clw[0][jobs] = missing;
        } else {
          clw[0][jobs] = (clw_obs_out[0][jobs] + clw_guess_out[0][jobs])/2.0;
        }
        break;
      default:
        oops::Log::error() << "Invalid value for clw_option:" << clw_option << std::endl;
        ABORT("Invalid value for clw_option");
    }

// Apply CLW threshold to observations
    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      if (apply[jobs]) {
        if (clw[0][jobs] == missing ||
            clw[0][jobs] > clw_thresholds[jv]) flagged[jv][jobs] = true;
      }
    }
  }
  clw.save("Derived");
  clw_obs_out.save("Derived");
  clw_guess_out.save("Derived");
}

// -----------------------------------------------------------------------------

float amsua_clw(float tobs1, float tobs2, float sza) {
    const float d1 = 0.754;
    const float d2 = -2.265;
    const float missing = util::missingValue<float>();

    float clw;

    if (tobs1 != missing && tobs2 != missing && sza != missing &&
        tobs1 <= 284.0 && tobs2 <= 284.0 && tobs1 > 0.0 && tobs2 > 0.0) {
      float cossza = cos(M_PI * sza/180.0);
      float d0 = 8.240 - (2.622 - 1.846*cossza)*cossza;
      clw = cossza*(d0 + d1*std::log(285.0-tobs1)) + d2*std::log(285.0-tobs2);
      clw = std::max(0.0f, clw);
    } else {
      clw = missing;
    }

    return clw;
}

// -----------------------------------------------------------------------------

void MWCLWCheck::print(std::ostream & os) const {
  os << "MWCLWCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
