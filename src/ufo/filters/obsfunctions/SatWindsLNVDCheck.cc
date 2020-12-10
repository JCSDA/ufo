/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SatWindsLNVDCheck.h"

#include <math.h>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<SatWindsLNVDCheck> makerObsFuncSatWindsLNVDCheck_("SatWindsLNVDCheck");

// -----------------------------------------------------------------------------

SatWindsLNVDCheck::SatWindsLNVDCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "SatWindsLNVDCheck: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // We need to retrieve the observed wind components.
  invars_ += Variable("eastward_wind@ObsValue");
  invars_ += Variable("northward_wind@ObsValue");

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string testHofX = options_.testHofX.value();
  invars_ += Variable("eastward_wind@" + testHofX);
  invars_ += Variable("northward_wind@" + testHofX);

  // TODO(gthompsn): Need to include a check that whatever HofX group name used actually exists.
}

// -----------------------------------------------------------------------------

SatWindsLNVDCheck::~SatWindsLNVDCheck() {}

// -----------------------------------------------------------------------------

void SatWindsLNVDCheck::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue(missing);

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  // Retrieve SatWinds observations of wind components
  std::vector<float> u, v;
  in.get(Variable("eastward_wind@ObsValue"), u);
  in.get(Variable("northward_wind@ObsValue"), v);
  // Retrieve Model HofX wind components
  std::string testHofX = options_.testHofX.value();
  std::vector<float> um, vm;
  in.get(Variable("eastward_wind@" + testHofX), um);
  in.get(Variable("northward_wind@" + testHofX), vm);

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (u[jj] != missing && v[jj] != missing) {
      out[0][jj] = sqrt((u[jj]-um[jj])*(u[jj]-um[jj]) + (v[jj]-vm[jj])*(v[jj]-vm[jj]))
                    / log(sqrt(u[jj]*u[jj] + v[jj]*v[jj]));
      oops::Log::debug() << "u, v: " << u[jj] << ", " << v[jj]
                         << " um, vm: " << um[jj] << ", " << vm[jj]
                         << " LNVD: " << out[0][jj] << std::endl;
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SatWindsLNVDCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
