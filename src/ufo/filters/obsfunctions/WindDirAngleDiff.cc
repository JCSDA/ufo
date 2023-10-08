/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/WindDirAngleDiff.h"

#include <algorithm>
#include <cmath>
#include <valarray>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<WindDirAngleDiff> makerObsFuncWindDirAngleDiff_("WindDirAngleDiff");

// -----------------------------------------------------------------------------

WindDirAngleDiff::WindDirAngleDiff(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "WindDirAngleDiff: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // We need to retrieve the observed wind components.
  invars_ += Variable("ObsValue/windEastward");
  invars_ += Variable("ObsValue/windNorthward");

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string test_hofx = options_.test_hofx.value();
  invars_ += Variable(test_hofx + "/windEastward");
  invars_ += Variable(test_hofx + "/windNorthward");

  // TODO(gthompsn): Need to include a check that whatever HofX group name used actually exists.
}

// -----------------------------------------------------------------------------

WindDirAngleDiff::~WindDirAngleDiff() {}

// -----------------------------------------------------------------------------

void WindDirAngleDiff::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();
  const double deg = Constants::rad2deg;

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  // Retrieve minimum_uv value and assure it is sensible.
  const float min_uv = std::max(0.0001f, options_.minimum_uv.value());

  // Retrieve observations of wind components
  std::vector<float> u, v;
  in.get(Variable("ObsValue/windEastward"), u);
  in.get(Variable("ObsValue/windNorthward"), v);
  // Retrieve Model HofX wind components
  std::string test_hofx = options_.test_hofx.value();
  std::vector<float> um, vm;
  in.get(Variable(test_hofx + "/windEastward"), um);
  in.get(Variable(test_hofx + "/windNorthward"), vm);

  double wdir_obs, wdir_model;

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (u[jj] != missing && v[jj] != missing) {
      if (std::abs(u[jj]) < min_uv && std::abs(v[jj]) < min_uv) {
        out[0][jj] = 0.0;
      } else {
        wdir_obs = std::atan2(-u[jj], -v[jj])*deg;
        wdir_model = std::atan2(-um[jj], -vm[jj])*deg;
        out[0][jj] = std::min({std::abs(wdir_obs-wdir_model),
                              std::abs(wdir_obs-wdir_model+360.0),
                              std::abs(wdir_obs-wdir_model-360.0)});
      }
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & WindDirAngleDiff::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
