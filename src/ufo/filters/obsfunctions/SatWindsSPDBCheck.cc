/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SatWindsSPDBCheck.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<SatWindsSPDBCheck> makerObsFuncSatWindsSPDBCheck_("SatWindsSPDBCheck");

// -----------------------------------------------------------------------------

SatWindsSPDBCheck::SatWindsSPDBCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "SatWindsSPDBCheck: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // Initialize error_min, and max from options. Make sure they are sane.
  const float error_min = options_.error_min.value();
  const float error_max = options_.error_max.value();
  ASSERT(error_min < error_max);

  // We need to retrieve the observed wind components.
  invars_ += Variable("eastward_wind@ObsValue");
  invars_ += Variable("northward_wind@ObsValue");

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string test_hofx = options_.test_hofx.value();
  invars_ += Variable("eastward_wind@" + test_hofx);
  invars_ += Variable("northward_wind@" + test_hofx);

  // The starting (un-inflated) value of obserror. If running in sequence of filters,
  // then it is probably found in ObsErrorData, otherwise, it is probably ObsError.
  const std::string errgrp = options_.original_obserr.value();
  invars_ += Variable("eastward_wind@"+errgrp);

  // TODO(gthompsn): Need to include a check that whatever HofX/obserr group name used exists.
}

// -----------------------------------------------------------------------------

SatWindsSPDBCheck::~SatWindsSPDBCheck() {}

// -----------------------------------------------------------------------------

void SatWindsSPDBCheck::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue(missing);
  float spdb = 0.0f;
  float residual = 0.0f;
  float obserr = 1.0f;

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  // Get min, max, gross error values
  const float error_min = options_.error_min.value();
  const float error_max = options_.error_max.value();

  // Retrieve SatWinds observations of wind components
  std::vector<float> u, v;
  in.get(Variable("eastward_wind@ObsValue"), u);
  in.get(Variable("northward_wind@ObsValue"), v);
  // Retrieve Model HofX wind components
  std::string test_hofx = options_.test_hofx.value();
  std::vector<float> um, vm;
  in.get(Variable("eastward_wind@" + test_hofx), um);
  in.get(Variable("northward_wind@" + test_hofx), vm);

  // Get original ObsError of eastward_wind (would make little sense if diff from northward)
  std::vector<float> currentObserr(nlocs);
  const std::string errgrp = options_.original_obserr.value();
  in.get(Variable("eastward_wind@"+errgrp), currentObserr);

  for (size_t jj = 0; jj < nlocs; ++jj) {
    out[0][jj] = 0.0f;
    if (u[jj] != missing && v[jj] != missing) {
      spdb = sqrt(u[jj]*u[jj]+v[jj]*v[jj]) - sqrt(um[jj]*um[jj]+vm[jj]*vm[jj]);
      if (spdb < 0.0f) {
        obserr = currentObserr[jj];
        obserr = std::max(error_min, std::min(obserr, error_max));
        residual = sqrt((u[jj]-um[jj])*(u[jj]-um[jj]) + (v[jj]-vm[jj])*(v[jj]-vm[jj]));
        out[0][jj] = residual/obserr;
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SatWindsSPDBCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
