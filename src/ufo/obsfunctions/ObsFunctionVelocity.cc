/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/obsfunctions/ObsFunctionVelocity.h"

#include <math.h>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionVelocity> makerObsFuncVelocity_("Velocity");

// -----------------------------------------------------------------------------

ObsFunctionVelocity::ObsFunctionVelocity()
  : obsvars_(), metadatavars_(), geovars_() {
  // needs two wind components to calculate speed
  obsvars_.push_back("eastward_wind");
  obsvars_.push_back("northward_wind");
}

// -----------------------------------------------------------------------------

ObsFunctionVelocity::~ObsFunctionVelocity() {}

// -----------------------------------------------------------------------------

void ObsFunctionVelocity::compute(const ioda::ObsDataVector<float> & metadata,
                                  const ioda::ObsDataVector<float> & obs,
                                  const GeoVaLs & geovals,
                                  ioda::ObsDataVector<float> & out) const {
  // TODO(AS): should use constants for variable names
  const size_t nlocs = obs.nlocs();
  for (size_t jj = 0; jj < nlocs; ++jj) {
    out[0][jj] = sqrt(pow(obs["eastward_wind"][jj], 2) + pow(obs["northward_wind"][jj], 2));
    oops::Log::debug() << "u, v: " << obs["eastward_wind"][jj] << ", "
                       << obs["northward_wind"][jj] << ", speed=" << out[0][jj] << std::endl;
  }
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsFunctionVelocity::requiredObsData() const {
  return obsvars_;
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsFunctionVelocity::requiredMetaData() const {
  return metadatavars_;
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsFunctionVelocity::requiredGeoVaLs() const {
  return geovars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
