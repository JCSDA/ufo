/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SIRetMW.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<SIRetMW> makerSIRetMW_("SIRetMW");

SIRetMW::SIRetMW(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Check required parameters
  // Get variable group types for SI retrieval from option
  ASSERT(options_.varGroup.value().size() == 1 || options_.varGroup.value().size() == 2);

  // Get channels for SI retrieval from options
  const std::vector<int> channels_ = {options_.ch90.value(), options_.ch150.value()};
  ASSERT(options_.ch90 != 0 && options_.ch150 != 0 && channels_.size() == 2);

  // Include list of required data from ObsSpace
  for (size_t igrp = 0; igrp < options_.varGroup.value().size(); ++igrp) {
    invars_ += Variable(options_.varGroup.value()[igrp] + "/brightnessTemperature", channels_);
  }
  invars_ += Variable(options_.testBias.value() + "/brightnessTemperature", channels_);

  // Include list of required data from ObsDiag
  invars_ += Variable("ObsDiag/brightness_temperature_assuming_clear_sky" , channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");
}

// -----------------------------------------------------------------------------

SIRetMW::~SIRetMW() {}

// -----------------------------------------------------------------------------

void SIRetMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get required parameters
  const std::vector<std::string> &vargrp = options_.varGroup.value();
  const std::vector<int> channels_ = {options_.ch90.value(), options_.ch150.value()};

  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t ngrps = vargrp.size();

  // Get variables from GeoVaLs

  // Get brightness temperature assuming all pixels are in clear skies
  std::vector<float> clr90(nlocs), clr150(nlocs);
  in.get(Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_)[0], clr90);
  in.get(Variable("ObsDiag/brightness_temperature_assuming_clear_sky", channels_)[1], clr150);

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);

  // Calculate scattering index
  std::vector<float> bt90(nlocs), bt150(nlocs);
  for (size_t igrp = 0; igrp < ngrps; ++igrp) {
    // Get data based on group type
    in.get(Variable(vargrp[igrp] + "/brightnessTemperature", channels_)[0], bt90);
    in.get(Variable(vargrp[igrp] + "/brightnessTemperature", channels_)[1], bt150);
    // Get bias based on group type
    if (options_.addBias.value() == vargrp[igrp]) {
      std::vector<float> bias90(nlocs), bias150(nlocs);
      if (in.has(Variable(options_.testBias.value() + "/brightnessTemperature", channels_)[0])) {
      in.get(Variable(options_.testBias.value() + "/brightnessTemperature", channels_)
                      [0], bias90);
      in.get(Variable(options_.testBias.value() + "/brightnessTemperature", channels_)
                      [1], bias150);
      } else {
      bias90.assign(nlocs, 0.0f);
      bias150.assign(nlocs, 0.0f);
      }
      // Add bias correction to the assigned group (only need to do it for ObsValue, since
      // HofX includes bias correction
      if (options_.addBias.value() == "ObsValue") {
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          bt90[iloc] = bt90[iloc] - bias90[iloc];
          bt150[iloc] = bt150[iloc] - bias150[iloc];
        }
      } else {
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          // Temporarily account for ZERO clear-sky BT output from CRTM
          if (clr90[iloc] > -1.0f && clr90[iloc] < 1.0f) {
            clr90[iloc] = bt90[iloc];
          } else {
            clr90[iloc] = clr90[iloc] + bias90[iloc];
          }
          if (clr150[iloc] > -1.0f && clr150[iloc] < 1.0f) {
            clr150[iloc] = bt150[iloc];
          } else {
            clr150[iloc] = clr150[iloc] + bias150[iloc];
          }
        }
      }
    }

    // Retrieve scattering index
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (water_frac[iloc] >= 0.99) {
        float siclr = clr90[iloc] - clr150[iloc];
        out[igrp][iloc] = bt90[iloc] - bt150[iloc] - siclr;
      } else {
        out[igrp][iloc] = bt90[iloc] - bt150[iloc];
      }
      }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SIRetMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
