/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorTopoRad.h"

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

static ObsFunctionMaker<ObsErrorFactorTopoRad>
       makerObsErrorFactorTopoRad_("ObsErrorFactorTopoRad");

// -----------------------------------------------------------------------------

ObsErrorFactorTopoRad::ObsErrorFactorTopoRad(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Include required variables from ObsDiag
  invars_ += Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("surface_geopotential_height@GeoVaLs");
}

// -----------------------------------------------------------------------------

ObsErrorFactorTopoRad::~ObsErrorFactorTopoRad() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorTopoRad::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();
  size_t nlevs = in.nlevs(Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_)[0]);

  // Get surface geopotential height
  std::vector<float> zsges(nlocs);
  in.get(Variable("surface_geopotential_height@GeoVaLs"), zsges);

  // Inflate obs error as a function of terrian height (>2000) and surface-to-space transmittance
  std::vector<float> tao_sfc(nlocs);
  for (size_t ich = 0; ich < nchans; ++ich) {
    in.get(Variable("transmittances_of_atmosphere_layer@ObsDiag", channels_)[ich], nlevs, tao_sfc);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ich][iloc] = 1.0;
      if (zsges[iloc] > 2000.0) {
        float factor = pow((2000.0/zsges[iloc]), 4);
        out[ich][iloc] = sqrt(1.0 / (1.0 - (1.0 - factor) * tao_sfc[iloc]));
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorTopoRad::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
