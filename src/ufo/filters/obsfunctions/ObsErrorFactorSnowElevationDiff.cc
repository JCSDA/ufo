/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorSnowElevationDiff.h"

#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorSnowElevationDiff> 
    makerSnowElevDiff_("ObsErrorFactorSnowElevationDiff");

// -----------------------------------------------------------------------------

ObsErrorFactorSnowElevationDiff::ObsErrorFactorSnowElevationDiff(
    const eckit::Configuration &config)
  : invars_() {
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff constructor" << std::endl;
  oops::Log::debug() << "ObsErrorFactorSnowElevationDiff: config = " << config << std::endl;
  
  // Initialize options
  options_.reset(new ObsErrorFactorSnowElevationDiffParameters());
  options_->deserialize(config);

  // Include required observation elevation variable
  const std::string obs_elev_var = options_->obs_elevation_var.value();
  invars_ += Variable(obs_elev_var);

  // Include required model elevation variable
  const std::string model_elev_var = options_->model_elevation_var.value();
  invars_ += Variable(model_elev_var);
}

// -----------------------------------------------------------------------------

ObsErrorFactorSnowElevationDiff::~ObsErrorFactorSnowElevationDiff() {
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorFactorSnowElevationDiff::compute(
    const ObsFilterData & data,
    ioda::ObsDataVector<float> & obserr) const {
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff compute start" << std::endl;
  
  const float missing = util::missingValue<float>();
  const float h_scale = options_->elevation_scale_h.value();

  // If no observations on this processor then nothing to do
  if (data.nlocs() == 0) return;

  // Ensure that only one output variable is expected
  ASSERT(obserr.nvars() == 1);

  // Get dimensions
  size_t nlocs = data.nlocs();

  // Get observation elevation variable
  std::vector<float> ob_elevation(nlocs);
  const std::string obs_elev_var = options_->obs_elevation_var.value();
  data.get(Variable(obs_elev_var), ob_elevation);

  // Get model surface elevation variable
  std::vector<float> model_elevation(nlocs);
  const std::string model_elev_var = options_->model_elevation_var.value();
  data.get(Variable(model_elev_var), model_elevation);

  // Compute inflation factor for each observation
  float dz, inflation_factor_1, inflation_factor;
  int iv = 0;

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // If missing observation or model elevation, set factor to 1.0 (no inflation)
    if (ob_elevation[iloc] == missing || model_elevation[iloc] == missing) {
      obserr[iv][iloc] = 1.0f;
    } else {
      // Compute elevation difference
      dz = std::abs(model_elevation[iloc] - ob_elevation[iloc]);
      
      // Compute inflation_factor_1 = exp(-1 * dz^2 / (h^2))
      inflation_factor_1 = std::exp(-1.0f * dz * dz / (h_scale * h_scale));
      
      // Compute output inflation_factor = 1 / inflation_factor_1
      inflation_factor = 1.0f / inflation_factor_1;
      
      // Return the inflation factor (ratio)
      obserr[iv][iloc] = inflation_factor;
    }
  }
  
  oops::Log::trace() << "ObsErrorFactorSnowElevationDiff compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorSnowElevationDiff::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
