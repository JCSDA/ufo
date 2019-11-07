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

namespace ufo {

// -----------------------------------------------------------------------------

MWCLWCheck::MWCLWCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                       boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                       boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr),
    invars_(eckit::LocalConfiguration(config_, "inputs"))
{
  oops::Log::debug() << "MWCLWCheck: config = " << config_ << std::endl;
  const Variable var0(invars_[0] + "@HofX");
  const Variable var1(invars_[1] + "@HofX");
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

  const float missing = util::missingValue(missing);
  float amsua_clw(float, float, float);

// Get config
  std::vector<float> clw_thresholds = config_.getFloatVector("clw_thresholds");
// clw_option controls how the clw is calculated:
//     1) Use observed BTs
//     2) Use calculated BTs
//     3) Symmetric calculation
  const int clw_option = config_.getInt("clw_option");
  oops::Log::debug() << "MWCLWCheck: config = " << config_ << std::endl;

// Number of channels to be tested and number of thresholds needs to be the same
  ASSERT(clw_thresholds.size() == filtervars.nvars());
// Check we have the correct number of channels to do the CLW calculation
  ASSERT(invars_.size() == 2);
// Check clw_option is in range
  ASSERT(clw_option >= 1 && clw_option <=3);

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsVariables(), "ObsValue");
  ioda::ObsDataVector<float> obs_for_calc(obsdb_, invars_, "ObsValue");
  ioda::ObsDataVector<float> sza(obsdb_, "sensor_zenith_angle", "MetaData");
  ioda::ObsDataVector<float> clw(obsdb_, "cloud_liquid_water", "Diagnostic", false);
  ioda::ObsDataVector<float> clw_guess_out(obsdb_, "clws_guess", "Diagnostic", false);
  ioda::ObsDataVector<float> clw_obs_out(obsdb_, "clw_obs", "Diagnostic", false);

//  H(x)
    const Variable var0(invars_[0] + "@HofX");
    const Variable var1(invars_[1] + "@HofX");
    std::vector<float> hofx0;
    data_.get(var0, hofx0);
    std::vector<float> hofx1;
    data_.get(var1, hofx1);

// Loop over obs locations calculating CLW from observations
  float clw_obs = missing, clw_guess = missing;
  for (size_t jobs = 0; jobs < obs.nlocs(); ++jobs) {
      switch (clw_option) {
           case 1 :
               clw_obs = amsua_clw(obs_for_calc[0][jobs], obs_for_calc[1][jobs],
                               sza[0][jobs]);
               clw[0][jobs] = clw_obs;
               break;
           case 2 :
               clw_guess = amsua_clw(hofx0[jobs], hofx1[jobs], sza[0][jobs]);
               clw[0][jobs] = clw_guess;
               break;
           case 3 :
               clw_obs = amsua_clw(obs_for_calc[0][jobs], obs_for_calc[1][jobs], sza[0][jobs]);
               clw_guess = amsua_clw(hofx0[jobs], hofx1[jobs], sza[0][jobs]);
               clw_obs_out[0][jobs] = clw_obs;
               clw_guess_out[0][jobs] = clw_guess;
               if (clw_obs != missing && clw_guess != missing)
               {
                 clw[0][jobs] = (clw_obs + clw_guess)/2.0;
               } else {
                 clw[0][jobs] = missing;
               }
               break;
      }

// Apply CLW threshold to observations
     for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
        if (apply[jobs]) {
           if (clw_thresholds[jv] != missing && (clw[0][jobs] == missing ||
               clw[0][jobs] > clw_thresholds[jv])) flagged[jv][jobs] = true;
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
    const float missing = util::missingValue(missing);

    float clw;
    float cossza = cos(M_PI * sza/180.0);
    float d0 = 8.240 - (2.622 - 1.846*cossza)*cossza;
    if (tobs1 <= 284.0 && tobs2 <= 284.0 && tobs1 > 0.0 && tobs2 > 0.0)
    {
        clw = cossza*(d0 + d1*log(285.0-tobs1)) + d2*log(285.0-tobs2);
        clw = std::max(0.0f, clw);
    } else {
        clw = missing;
    }

    return clw;
}

// -----------------------------------------------------------------------------

void MWCLWCheck::print(std::ostream & os) const {
  os << "MWCLWCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
