/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <algorithm>
#include <cmath>
#include <valarray>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/WindUVfromObsSpdHofxDir.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<WindUVfromObsSpdHofxDir>
makerObsFuncWindUVfromObsSpdHofxDir_("WindUVfromObsSpdHofxDir");

// -----------------------------------------------------------------------------

WindUVfromObsSpdHofxDir::WindUVfromObsSpdHofxDir(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Retrieve windspeed found through a retrieval method (direction not included)
  invars_ += Variable("ObsValue/windSpeed");

  // Typical use would be HofX group, but during testing, we include option for other
  // reference H(x) groups
  const std::string &hofxgrp = options_.hofxGroup.value();
  invars_ += Variable(hofxgrp + "/windEastward");
  invars_ += Variable(hofxgrp + "/windNorthward");
}

// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------

void WindUVfromObsSpdHofxDir::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();
  const float rad2deg = 1.0f / Constants::deg2rad;

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  ioda::ObsSpace &obsdb_ = in.obsspace();


  // Retrieve wind components from hofx results. Used to calculate direction
  std::vector<float> wspd;
  in.get(Variable("ObsValue/windSpeed"), wspd);
  // Retrieve Model HofX wind components
  const std::string &hofxgrp = options_.hofxGroup.value();
  std::vector<float> um, vm;
  in.get(Variable(hofxgrp + "/windEastward"), um);
  in.get(Variable(hofxgrp + "/windNorthward"), vm);

  float wdir_model;
  std::vector<float> w_east(nlocs);
  std::vector<float> w_north(nlocs);


  // Wind direction is not provided in the obs file. As a proxy, we take the
  // modeled direction to represent an observed direction.
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (wspd[jj] != missing) {
      // compute wind direction in radian
      wdir_model = std::atan2(-um[jj], -vm[jj]);
      // compute wind components
      w_east[jj] = -wspd[jj] * std::sin(wdir_model);
      w_north[jj] = -wspd[jj] * std::cos(wdir_model);
      // save wind direction in degree
      out[0][jj] = wdir_model * rad2deg;
    } else {
      out[0][jj] = missing;
      w_east[jj] = missing;
      w_north[jj] = missing;
    }
  }
  out.save("DerivedValue");

  // Store wind components to the DerivedObsValue group
  obsdb_.put_db("DerivedObsValue", "windEastward", w_east);
  obsdb_.put_db("DerivedObsValue", "windNorthward", w_north);
}

// -----------------------------------------------------------------------------

const ufo::Variables & WindUVfromObsSpdHofxDir::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
