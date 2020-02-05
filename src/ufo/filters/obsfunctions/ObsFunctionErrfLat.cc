/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionErrfLat.h"

#include <math.h>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionErrfLat> makerObsFuncErrfLat_("ErrfLat");

// -----------------------------------------------------------------------------

ObsFunctionErrfLat::ObsFunctionErrfLat(const eckit::LocalConfiguration conf)
  : invars_(), conf_(conf) {
  // Check options
  ASSERT(conf_.has("latitude_parameters"));
  invars_ += Variable("latitude@MetaData");
}

// -----------------------------------------------------------------------------

ObsFunctionErrfLat::~ObsFunctionErrfLat() {}

// -----------------------------------------------------------------------------

void ObsFunctionErrfLat::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get parameters from options
  // Parameters for reducing observation error bounds within latitude band defined by params[0]
  // params[1-3] give the maximum reduction at equator and decreasing towards params[0]
  std::vector<float> params = conf_.getFloatVector("latitude_parameters");

  const size_t nlocs = in.nlocs();
  std::vector<float> lats;
  in.get(Variable("latitude@MetaData"), lats);
  for (size_t jj = 0; jj < nlocs; ++jj) {
    out[0][jj] = 1.0;
    if ( std::abs(lats[jj]) < params[0] ) {
      out[0][jj] = params[1] *(std::abs(lats[jj]) * params[2] + params[3]);
    }
    out[0][jj] = sqrt(1.0 / out[0][jj]);
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionErrfLat::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
