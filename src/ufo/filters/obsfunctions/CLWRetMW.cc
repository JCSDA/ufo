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
  ASSERT((options_.ch238.value() != boost::none && options_.ch314.value() != boost::none) ||
         (options_.ch37v.value() != boost::none && options_.ch37h.value() != boost::none));

  if (options_.ch238.value() != boost::none && options_.ch314.value() != boost::none) {
    // For AMSUA and ATMS retrievals
    // Get channels for CLW retrieval from options
    const std::vector<int> channels = {options_.ch238.value().get(), options_.ch314.value().get()};
    ASSERT(options_.ch238.value().get() != 0 && options_.ch314.value().get() != 0
           && channels.size() == 2);
    // Include list of required data from ObsSpace
    for (size_t igrp = 0; igrp < options_.varGroup.value().size(); ++igrp) {
      invars_ += Variable("brightness_temperature@" + options_.varGroup.value()[igrp], channels);
    }
    invars_ += Variable("brightness_temperature@" + options_.testBias.value(), channels);
    invars_ += Variable("sensor_zenith_angle@MetaData");

    // Include list of required data from GeoVaLs
    invars_ += Variable("average_surface_temperature_within_field_of_view@GeoVaLs");
    invars_ += Variable("water_area_fraction@GeoVaLs");
    invars_ += Variable("surface_temperature_where_sea@GeoVaLs");

  } else if (options_.ch37v.value() != boost::none && options_.ch37h.value() != boost::none) {
    // For cloud index like GMI's.
    // Get channels for CLW retrieval from options
    const std::vector<int> channels = {options_.ch37v.value().get(), options_.ch37h.value().get()};

    ASSERT(options_.ch37v.value().get() != 0 && options_.ch37h.value().get() != 0 &&
           channels.size() == 2);
    // Include list of required data from ObsSpace
    for (size_t igrp = 0; igrp < options_.varGroup.value().size(); ++igrp) {
      invars_ += Variable("brightness_temperature@" + options_.varGroup.value()[igrp], channels);
    }
    invars_ += Variable("brightness_temperature@" + options_.testBias.value(), channels);
    // Include list of required data from ObsDiag
    invars_ += Variable("brightness_temperature_assuming_clear_sky@ObsDiag" , channels);

    // Include list of required data from GeoVaLs
    invars_ += Variable("water_area_fraction@GeoVaLs");
  }
}

// -----------------------------------------------------------------------------

CLWRetMW::~CLWRetMW() {}

// -----------------------------------------------------------------------------

void CLWRetMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get required parameters
  const std::vector<std::string> &vargrp = options_.varGroup.value();

  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t ngrps = vargrp.size();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // --------------- amsua or atms --------------------------
  if (options_.ch238.value() != boost::none && options_.ch314.value() != boost::none) {
    const std::vector<int> channels = {options_.ch238.value().get(), options_.ch314.value().get()};

    // Get variables from ObsSpace
    // Get sensor zenith angle
    std::vector<float> szas(nlocs);
    in.get(Variable("sensor_zenith_angle@MetaData"), szas);

    // Get variables from GeoVaLs
    // Get average surface temperature in FOV
    std::vector<float> tsavg(nlocs);
    in.get(Variable("average_surface_temperature_within_field_of_view@GeoVaLs"), tsavg);

    // Calculate retrieved cloud liquid water
    std::vector<float> bt238(nlocs), bt314(nlocs);
    for (size_t igrp = 0; igrp < ngrps; ++igrp) {
      // Get data based on group type
      in.get(Variable("brightness_temperature@" + vargrp[igrp], channels)[0], bt238);
      in.get(Variable("brightness_temperature@" + vargrp[igrp], channels)[1], bt314);
      // Get bias based on group type
      if (options_.addBias.value() == vargrp[igrp]) {
        std::vector<float> bias238(nlocs), bias314(nlocs);
        if (in.has(Variable("brightness_temperature@" + options_.testBias.value(), channels)[0])) {
        in.get(Variable("brightness_temperature@" + options_.testBias.value(), channels)[0],
                        bias238);
        in.get(Variable("brightness_temperature@" + options_.testBias.value(), channels)[1],
                        bias314);
        } else {
        bias238.assign(nlocs, 0.0f);
        bias314.assign(nlocs, 0.0f);
        }
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

      // Compute the cloud liquid water
      cloudLiquidWater(szas, tsavg, water_frac, bt238, bt314, out[igrp]);
    }
  // -------------------- GMI ---------------------------
  } else if (options_.ch37v.value() != boost::none && options_.ch37h.value() != boost::none) {
    const std::vector<int> channels = {options_.ch37v.value().get(), options_.ch37h.value().get()};
    // Indices of data at channeles 37v and 37h in the above array "channels"
    const int jch37v = 0;
    const int jch37h = 1;
    std::vector<float> bt_clr_37v(nlocs), bt_clr_37h(nlocs);
    in.get(Variable("brightness_temperature_assuming_clear_sky@ObsDiag" , channels)
           [jch37v], bt_clr_37v);
    in.get(Variable("brightness_temperature_assuming_clear_sky@ObsDiag" , channels)
           [jch37h], bt_clr_37h);
    // Calculate retrieved cloud liquid water
    std::vector<float> bt37v(nlocs), bt37h(nlocs);
    for (size_t igrp = 0; igrp < ngrps; ++igrp) {
      // Get data based on group type
      in.get(Variable("brightness_temperature@" + vargrp[igrp], channels) [jch37v], bt37v);
      in.get(Variable("brightness_temperature@" + vargrp[igrp], channels) [jch37h], bt37h);
      // Get bias based on group type
      if (options_.addBias.value() == vargrp[igrp]) {
        std::vector<float> bias37v(nlocs), bias37h(nlocs);
        if (in.has(Variable("brightness_temperature@" + options_.testBias.value(), channels)
            [jch37v])) {
          in.get(Variable("brightness_temperature@" + options_.testBias.value(), channels)
                 [jch37v], bias37v);
          in.get(Variable("brightness_temperature@" + options_.testBias.value(), channels)
                 [jch37h], bias37h);
        } else {
        bias37v.assign(nlocs, 0.0f);
        bias37h.assign(nlocs, 0.0f);
        }
        // Add bias correction to the assigned group
        if (options_.addBias.value() == "ObsValue") {
          for (size_t iloc = 0; iloc < nlocs; ++iloc) {
            bt37v[iloc] = bt37v[iloc] - bias37v[iloc];
            bt37h[iloc] = bt37h[iloc] - bias37h[iloc];
          }
        } else {
          for (size_t iloc = 0; iloc < nlocs; ++iloc) {
            bt37v[iloc] = bt37v[iloc] + bias37v[iloc];
            bt37h[iloc] = bt37h[iloc] + bias37h[iloc];
          }
        }
      }


      // Compute cloud index
      CIret_37v37h_diff(bt_clr_37v, bt_clr_37h, water_frac, bt37v, bt37h, out[igrp]);
    }
  }
}

// -----------------------------------------------------------------------------

void CLWRetMW::cloudLiquidWater(const std::vector<float> & szas,
                                         const std::vector<float> & tsavg,
                                         const std::vector<float> & water_frac,
                                         const std::vector<float> & bt238,
                                         const std::vector<float> & bt314,
                                         std::vector<float> & out) {
  ///
  /// \brief Retrieve cloud liquid water from AMSU-A 23.8 GHz and 31.4 GHz channels.
  ///
  /// Reference: Grody et al. (2001)
  /// Determination of precipitable water and cloud liquid water over oceans from
  /// the NOAA 15 advanced microwave sounding unit
  ///
  const float t0c = Constants::t0c;
  const float d1 = 0.754, d2 = -2.265;
  const float c1 = 8.240, c2 = 2.622, c3 = 1.846;
  for (size_t iloc = 0; iloc < water_frac.size(); ++iloc) {
    if (water_frac[iloc] >= 0.99) {
      float cossza = cos(Constants::deg2rad * szas[iloc]);
      float d0 = c1 - (c2 - c3 * cossza) * cossza;
      if (tsavg[iloc] > t0c - 1.0 && bt238[iloc] <= 284.0 && bt314[iloc] <= 284.0
                                  && bt238[iloc] > 0.0 && bt314[iloc] > 0.0) {
        out[iloc] = cossza * (d0 + d1 * std::log(285.0 - bt238[iloc])
                                 + d2 * std::log(285.0 - bt314[iloc]));
        out[iloc] = std::max(0.f, out[iloc]);
      } else {
        out[iloc] = getBadValue();
      }
    }
  }
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

void CLWRetMW::CIret_37v37h_diff(const std::vector<float> & bt_clr_37v,
                                         const std::vector<float> & bt_clr_37h,
                                         const std::vector<float> & water_frac,
                                         const std::vector<float> & bt37v,
                                         const std::vector<float> & bt37h,
                                         std::vector<float> & out) {
  ///
  /// \brief Retrieve cloud index from GMI 37V and 37H channels.
  ///
  /// GMI cloud index: 1.0 - (Tb_37v - Tb_37h)/(Tb_37v_clr - Tb_37h_clr), in which
  /// Tb_37v_clr and Tb_37h_clr for calculated Tb at 37 V and 37H GHz from module values
  /// assuming in clear-sky condition. Tb_37v and Tb_37h are Tb observations at 37 V and 37H GHz.
  ///
  for (size_t iloc = 0; iloc < water_frac.size(); ++iloc) {
    if (water_frac[iloc] >= 0.99) {
      if (bt37h[iloc] <= bt37v[iloc]) {
        out[iloc] = 1.0 - (bt37v[iloc] - bt37h[iloc])/(bt_clr_37v[iloc] - bt_clr_37h[iloc]);
        out[iloc] = std::max(0.f, out[iloc]);
      } else {
        out[iloc] = getBadValue();
      }
    } else {
      out[iloc] = getBadValue();
    }
  }
}

// -----------------------------------------------------------------------------
const ufo::Variables & CLWRetMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
