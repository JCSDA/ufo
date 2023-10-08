/*
 * (C) Copyright 2023 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "ufo/filters/obsfunctions/SatWindsErrnormCheck.h"

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

static ObsFunctionMaker<SatWindsErrnormCheck> makerObsFuncSatWindsErrnormCheck\
_("SatWindsErrnormCheck");

// -----------------------------------------------------------------------------

SatWindsErrnormCheck::SatWindsErrnormCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "SatWindsErrnormCheck: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // We need to retrieve the observed wind components.
  invars_ += Variable("ObsValue/windEastward");
  invars_ += Variable("ObsValue/windNorthward");
  // We need to retrieve the expected error
  invars_ += Variable("MetaData/expectedError");
}

// -----------------------------------------------------------------------------

SatWindsErrnormCheck::~SatWindsErrnormCheck() {}

// -----------------------------------------------------------------------------

void SatWindsErrnormCheck::compute(const ObsFilterData & in,
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
  // Retrieve observation expected error:
  // This value is packed in the obs file in the same format as
  // a quality indicator value, with discrete values between
  // 0. and 100. A value > 100. indicates a missing value. These
  // packed values can be unpacked into an expected observation
  // error value in m/s via: 10. - 0.1*ee. A missing value for
  // ee assigned to > 100. will generate a negative unpacked
  // expected error value, which needs to be accounted for
  // in the errnorm test.
  std::vector<float> ee;
  in.get(Variable("MetaData/expectedError"), ee);
  // Variables to store exError, windSpeed, and errnorm
  double exError;
  double windSpeed;
  double errnorm;

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (u[jj] != missing && v[jj] != missing) {
      if (std::abs(u[jj]) < min_uv && std::abs(v[jj]) < min_uv) {
        out[0][jj] = 0.0;
      } else {
        // Define windSpeed from wind components
        windSpeed = std::sqrt(std::pow(u[jj], 2.0) + std::pow(v[jj], 2.0));
        // Define exError from ee
        exError = 10.0 - 0.1*ee[jj];
        // Define errnorm as ratio of exError to windSpeed:
        // Automatically assign a value of missing if:
        //    (1) windSpeed <= 0.1, or:
        //    (2) exError < 0.
        // The first indicates a wind too slow to generate a comparison,
        // the second indicates a missing-value for ee generated a negative
        // exError value, and should be rejected for missing data
        if (windSpeed > 0.1 || exError >= 0.) {
          errnorm = exError/windSpeed;
        } else {
          errnorm = missing;
        }
        // Output errnrom
        out[0][jj] = errnorm;
      }
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SatWindsErrnormCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
