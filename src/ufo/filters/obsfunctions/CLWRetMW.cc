/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/CLWRetMW.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<CLWRetMW> makerCLWRetMW_("CLWRetMW");

CLWRetMW::CLWRetMW(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Check required parameters
  // Get variable group types for CLW retrieval from option
  ASSERT(options_.varGroup.value().size() == 1 || options_.varGroup.value().size() == 2);

  // Get channels for CLW retrieval from options
  const std::vector<int> channels_ = {options_.ch238.value(), options_.ch314.value()};
  ASSERT(options_.ch238 != 0 && options_.ch314 != 0 && channels_.size() == 2);

  // Include list of required data from ObsSpace
  for (size_t igrp = 0; igrp < options_.varGroup.value().size(); ++igrp) {
    invars_ += Variable("brightness_temperature@" + options_.varGroup.value()[igrp], channels_);
  }
  invars_ += Variable("brightness_temperature@" + options_.testBias.value(), channels_);
  invars_ += Variable("sensor_zenith_angle@MetaData");

  // Include list of required data from GeoVaLs
  invars_ += Variable("average_surface_temperature_within_field_of_view@GeoVaLs");
  invars_ += Variable("water_area_fraction@GeoVaLs");
  invars_ += Variable("surface_temperature_where_sea@GeoVaLs");
}

// -----------------------------------------------------------------------------

CLWRetMW::~CLWRetMW() {}

// -----------------------------------------------------------------------------

void CLWRetMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get required parameters
  const std::vector<std::string> &vargrp_ = options_.varGroup.value();
  const std::vector<int> channels_ = {options_.ch238.value(), options_.ch314.value()};

  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t ngrps = vargrp_.size();

  // Get variables from ObsSpace
  // Get sensor zenith angle
  std::vector<float> szas(nlocs);
  in.get(Variable("sensor_zenith_angle@MetaData"), szas);

  // Get variables from GeoVaLs
  // Get average surface temperature in FOV
  std::vector<float> tsavg(nlocs);
  in.get(Variable("average_surface_temperature_within_field_of_view@GeoVaLs"), tsavg);

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // Calculate retrieved cloud liquid water
  std::vector<float> bt238(nlocs), bt314(nlocs);
  for (size_t igrp = 0; igrp < ngrps; ++igrp) {
    // Get data based on group type
    in.get(Variable("brightness_temperature@" + vargrp_[igrp], channels_)[channels_[0]-1], bt238);
    in.get(Variable("brightness_temperature@" + vargrp_[igrp], channels_)[channels_[1]-1], bt314);
    // Get bias based on group type
    if (options_.addBias.value() == vargrp_[igrp]) {
      std::vector<float> bias238(nlocs), bias314(nlocs);
      in.get(Variable("brightness_temperature@" + options_.testBias.value(), channels_)
                      [channels_[0]-1], bias238);
      in.get(Variable("brightness_temperature@" + options_.testBias.value(), channels_)
                      [channels_[1]-1], bias314);
      // Add bias correction to the assigned group
      if (options_.addBias.value() == "ObsValue") {
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          bt238[iloc] = bt238[iloc] - bias238[iloc];
          bt314[iloc] = bt314[iloc] - bias314[iloc];
        }
      } else {
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          bt238[iloc] = bt238[iloc] + bias238[iloc];
          bt314[iloc] = bt314[iloc] + bias314[iloc];
        }
      }
    }
    const float t0c = Constants::t0c;
    const float d1 = 0.754, d2 = -2.265;
    const float c1 = 8.240, c2 = 2.622, c3 = 1.846;
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (water_frac[iloc] >= 0.99) {
        float cossza = cos(Constants::deg2rad * szas[iloc]);
        float d0 = c1 - (c2 - c3 * cossza) * cossza;
        if (tsavg[iloc] > t0c - 1.0 && bt238[iloc] <= 284.0 && bt314[iloc] <= 284.0
                                    && bt238[iloc] > 0.0 && bt314[iloc] > 0.0) {
          out[igrp][iloc] = cossza * (d0 + d1 * std::log(285.0 - bt238[iloc])
                                          + d2 * std::log(285.0 - bt314[iloc]));
          out[igrp][iloc] = std::max(0.f, out[igrp][iloc]);
        } else {
          out[igrp][iloc] = getBadValue();
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & CLWRetMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
