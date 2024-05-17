/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/CLWMatchIndexMW.h"

#include <cmath>

#include <algorithm>
#include <set>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<CLWMatchIndexMW> makerCLWMatchIndexMW_("CLWMatchIndexMW");

// -----------------------------------------------------------------------------

CLWMatchIndexMW::CLWMatchIndexMW(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");

  const Variable &clwobs = options_.clwobsFunction.value();
  invars_ += clwobs;

  const Variable &clwbkg = options_.clwbkgFunction.value();
  invars_ += clwbkg;
}

// -----------------------------------------------------------------------------

CLWMatchIndexMW::~CLWMatchIndexMW() {}

// -----------------------------------------------------------------------------

void CLWMatchIndexMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t nchans = channels_.size();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);

  // Get CLW retrieval based on observation from ObsFunction
  const Variable &clwobsvar = options_.clwobsFunction.value();
  ioda::ObsDataVector<float> clwobs(in.obsspace(), clwobsvar.toOopsObsVariables());
  in.get(clwobsvar, clwobs);

  // Get CLW retrieval based on simulated observation from ObsFunction
  const Variable &clwbkgvar = options_.clwbkgFunction.value();
  ioda::ObsDataVector<float> clwbkg(in.obsspace(), clwbkgvar.toOopsObsVariables());
  in.get(clwbkgvar, clwbkg);

  // Get parameters for observation errors from options
  const std::vector<float> &clw_clr = options_.clwretClearSky.value();
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ichan][iloc] = 1.0;
      if (water_frac[iloc] > 0.99) {
        float condition1 = (clwobs[0][iloc] - clw_clr[ichan]) * (clwbkg[0][iloc] - clw_clr[ichan]);
        float condition2 = std::abs(clwobs[0][iloc] - clwbkg[0][iloc]);
        if ( condition1 < 0 && condition2 >= 0.005) out[ichan][iloc] = 0.0;
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & CLWMatchIndexMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
