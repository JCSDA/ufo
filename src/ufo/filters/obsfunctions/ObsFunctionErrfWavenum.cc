/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionErrfWavenum.h"

#include <math.h>

#include <algorithm>
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

static ObsFunctionMaker<ObsFunctionErrfWavenum> makerObsFuncErrfWavenum_("ErrfWavenum");

// -----------------------------------------------------------------------------

ObsFunctionErrfWavenum::ObsFunctionErrfWavenum(const eckit::LocalConfiguration conf)
  : invars_(), channels_(), conf_(conf) {
  // Check options
  ASSERT(conf_.has("channels"));

  // Get channels from options
  const std::string chlist = conf_.getString("channels");
  std::set<int> channelset = oops::parseIntSet(chlist);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));

  // Include required variables from ObsDiag
  invars_ += Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("water_area_fraction@GeoVaLs");

  // Include list of required data from ObsSpace
  invars_ += Variable("solar_zenith_angle@MetaData");
  invars_ += Variable("sensor_band_central_radiation_wavenumber@VarMetaData");
}

// -----------------------------------------------------------------------------

ObsFunctionErrfWavenum::~ObsFunctionErrfWavenum() {}

// -----------------------------------------------------------------------------

void ObsFunctionErrfWavenum::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_)[0]);

  // Get surface geopotential height
  std::vector<float> water_frac(nlocs);
  in.get(Variable("water_area_fraction@GeoVaLs"), water_frac);

  // Get sensor zenith angle
  std::vector<float> solza(nlocs);
  in.get(Variable("solar_zenith_angle@MetaData"), solza);

  // Get sensor band central radiation wavenumber
  std::vector<float> wavenumber(nchans);
  in.get(Variable("sensor_band_central_radiation_wavenumber@VarMetaData"), wavenumber);

  // Inflate obs error for wavenumber in the range of (2000, 2400] during daytime over water surface
  // as a function of wavenumber number, surface-to-space transmittance, solar zenith angle, and
  // surface type
  std::vector<float> tao_sfc(nlocs);
  for (size_t ich = 0; ich < nchans; ++ich) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) out[ich][iloc] = 1.0; 
    if (wavenumber[ich] > 2000. && wavenumber[ich] <= 2400.0) {
      in.get(Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_)[ich], nlevs, tao_sfc);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        if (water_frac[iloc] > 0.0 && solza[iloc] <= 89.0) {
          float factor = std::max(0.0, cos(Constants::deg2rad * solza[iloc]));
          factor = tao_sfc[iloc] * factor *(1.0 / 400.0);
          out[ich][iloc] = sqrt(1.0 / (1.0 - (wavenumber[ich] - 2000.0) * factor));
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionErrfWavenum::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
