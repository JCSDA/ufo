/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorLatRad.h"

#include <math.h>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorLatRad> makerObsErrorFactorLatRad_("ObsErrorFactorLatRad");

// -----------------------------------------------------------------------------

ObsErrorFactorLatRad::ObsErrorFactorLatRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  ASSERT((options_.latitudeParameters.value()).size() == 4);

  invars_ += Variable("latitude@MetaData");
}

// -----------------------------------------------------------------------------

ObsErrorFactorLatRad::~ObsErrorFactorLatRad() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorLatRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get parameters from options
  const std::vector<float> &params = options_.latitudeParameters.value();

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

const ufo::Variables & ObsErrorFactorLatRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
