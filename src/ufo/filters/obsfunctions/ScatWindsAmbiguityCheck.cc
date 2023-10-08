/*
 * (C) Copyright 2022 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "ufo/filters/obsfunctions/ScatWindsAmbiguityCheck.h"

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

static ObsFunctionMaker<ScatWindsAmbiguityCheck> makerObsFuncScatWindsAmbiguityCheck\
_("ScatWindsAmbiguityCheck");

// -----------------------------------------------------------------------------

ScatWindsAmbiguityCheck::ScatWindsAmbiguityCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "ScatWindsAmbiguityCheck: config = " << conf << std::endl;
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

ScatWindsAmbiguityCheck::~ScatWindsAmbiguityCheck() {}

// -----------------------------------------------------------------------------

void ScatWindsAmbiguityCheck::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

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

  double vecdiff_obs, vecdiff_opp;

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (u[jj] != missing && v[jj] != missing) {
      if (std::abs(u[jj]) < min_uv && std::abs(v[jj]) < min_uv) {
        out[0][jj] = 0.0;
      } else {
        vecdiff_obs = std::sqrt(std::pow(u[jj]-um[jj], 2.0) + std::pow(v[jj]-vm[jj], 2.0));
        vecdiff_opp = std::sqrt(std::pow(-u[jj]-um[jj], 2.0) + std::pow(-v[jj]-vm[jj], 2.0));
        out[0][jj] = vecdiff_obs-vecdiff_opp;
      }
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ScatWindsAmbiguityCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
