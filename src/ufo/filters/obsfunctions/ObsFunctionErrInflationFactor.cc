/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionErrInflationFactor.h"

#include <math.h>
#include <algorithm>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionErrInflationFactor>
       makerObsFuncErrInflationFactor_("ErrInflationFactor");

// -----------------------------------------------------------------------------

ObsFunctionErrInflationFactor::ObsFunctionErrInflationFactor()
  : geovars_() {
}

// -----------------------------------------------------------------------------

ObsFunctionErrInflationFactor::~ObsFunctionErrInflationFactor() {}

// -----------------------------------------------------------------------------

void ObsFunctionErrInflationFactor::compute(const ObsFilterData & input,
                                            ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = input.nlocs();
  std::vector<float> lats = input.get("latitude@MetaData");
  for (size_t jj = 0; jj < nlocs; ++jj) {
      out[0][jj] = 1.0;
     if ( std::abs(lats[jj]) > 25.0 ) {
        out[0][jj] = 0.5 *(std::abs(lats[jj]) * 0.04 + 1.0);
        oops::Log::debug() << "latitude: " << lats[jj]
                           << ", Inflation Factor = " << out[0][jj] << std::endl;
     } else {
        oops::Log::debug() << "latitude: " << lats[jj]
                           << ", Inflation Factor = " << out[0][jj] << std::endl;
     }
  }
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsFunctionErrInflationFactor::requiredGeoVaLs() const {
  return geovars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
