/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionCLWRet.h"

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

static ObsFunctionMaker<ObsFunctionCLWRet> makerObsFuncCLWRet_("CLWRet");

ObsFunctionCLWRet::ObsFunctionCLWRet(const eckit::LocalConfiguration conf)
  : invars_(), group_() {
  // Check options
  ASSERT(conf.has("clwret_type"));

  // Get group type from option
  group_ = conf.getStringVector("clwret_type");

  // Set channels
  std::vector<int> channels{1, 2};

  // Include list of required data from ObsSpace and HofX (or GsiHofX for testing)
  for (size_t igrp = 0; igrp < group_.size(); ++igrp) {
    invars_ += Variable("brightness_temperature@" + group_[igrp], channels);
    if (group_[igrp] == "HofX" || group_[igrp] == "GsiHofX") {
      invars_ += Variable("brightness_temperature@ObsBias", channels);
    }
  }
  invars_ += Variable("sensor_zenith_angle@MetaData");

  // Include list of required data from GeoVaLs
  invars_ += Variable("average_surface_temperature_within_field_of_view@GeoVaLs");
  invars_ += Variable("water_area_fraction@GeoVaLs");
  invars_ += Variable("surface_temperature_where_sea@GeoVaLs");
}

// -----------------------------------------------------------------------------

ObsFunctionCLWRet::~ObsFunctionCLWRet() {}

// -----------------------------------------------------------------------------

void ObsFunctionCLWRet::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t ngrps = group_.size();

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
  std::vector<float> bt1(nlocs), bt2(nlocs);
  oops::Log::debug() << "ObsFunctionCLWRet: ngrps = " << ngrps << std::endl;
  for (size_t igrp = 0; igrp < ngrps; ++igrp) {
    // Get data based on group type
    in.get(Variable("brightness_temperature_1@"+group_[igrp]), bt1);
    in.get(Variable("brightness_temperature_2@"+group_[igrp]), bt2);

    // Get bias and bias corrected obs
    if (group_[igrp] == "HofX" || group_[igrp] == "GsiHofX") {
      std::vector<float> bias1(nlocs), bias2(nlocs);
      in.get(Variable("brightness_temperature_1@ObsBias"), bias1);
      in.get(Variable("brightness_temperature_2@ObsBias"), bias2);

      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        bt1[iloc] = bt1[iloc] + bias1[iloc];
        bt2[iloc] = bt2[iloc] + bias2[iloc];
      }
    }
    const float t0c = Constants::t0c;
    const float d1 = 0.754, d2 = -2.265;
    const float c1 = 8.240, c2 = 2.622, c3 = 1.846;
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (water_frac[iloc] >= 0.99) {
        float cossza = cos(Constants::deg2rad * szas[iloc]);
        float d0 = c1 - (c2 - c3 * cossza) * cossza;
        if (tsavg[iloc] > t0c - 1.f && bt1[iloc] <= 284.f && bt2[iloc] <= 284.f
                                    && bt1[iloc] > 0.f && bt2[iloc] > 0.f) {
          out[igrp][iloc] = cossza * (d0 + d1 * std::log(285.f - bt1[iloc])
                                          + d2 * std::log(285.f - bt2[iloc]));
          out[igrp][iloc] = std::fmax(0.f, out[igrp][iloc]);
        } else {
          out[igrp][iloc] = getBadValue();
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionCLWRet::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
