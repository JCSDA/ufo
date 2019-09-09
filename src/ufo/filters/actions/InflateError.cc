/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/InflateError.h"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<InflateError> makerInflateErr_("inflate error");

// -----------------------------------------------------------------------------

InflateError::InflateError(const eckit::Configuration & conf)
  : strfactor_(conf.getString("inflation")) {
}

// -----------------------------------------------------------------------------

void InflateError::apply(const oops::Variables & vars,
                         const std::vector<std::vector<bool>> & flagged,
                         const ObsFilterData & data,
                         ioda::ObsDataVector<int> &,
                         ioda::ObsDataVector<float> & obserr) const {
  std::vector<float> factors = getScalarOrFilterData(strfactor_, data);
  oops::Log::debug() << " input obserr: " << obserr << std::endl;
  ASSERT(factors.size() == obserr.nlocs());
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = obserr.varnames().find(vars[jv]);
    for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
      if (flagged[iv][jobs]) obserr[iv][jobs] *= factors[jobs];
    }
  }
  oops::Log::debug() << " output obserr: " << obserr << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
