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

namespace ufo {

static ObsFunctionMaker<ObsFunctionErrfLat> makerObsFuncErrfLat_("ErrfLat");

// -----------------------------------------------------------------------------

ObsFunctionErrfLat::ObsFunctionErrfLat()
  : invars_() {
  invars_ += "latitude@MetaData";
}

// -----------------------------------------------------------------------------

ObsFunctionErrfLat::~ObsFunctionErrfLat() {}

// -----------------------------------------------------------------------------

void ObsFunctionErrfLat::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // TODO(AS): should use constants for variable names
  const size_t nlocs = in.nlocs();
  std::vector<float> lats;
  in.get("latitude@MetaData", lats);
  for (size_t jj = 0; jj < nlocs; ++jj) {
    out[0][jj] = 1.0;
    if ( std::abs(lats[jj]) < 25.0 ) {
      out[0][jj] = 0.5 *(std::abs(lats[jj]) * 0.04 + 1.0);
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionErrfLat::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
