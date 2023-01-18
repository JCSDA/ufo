/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/TotalColumnVaporGuess.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<TotalColumnVaporGuess> makerTotalColumnVaporGuess_("TotalColumnVaporGuess");

TotalColumnVaporGuess::TotalColumnVaporGuess(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/humidity_mixing_ratio");
  invars_ += Variable("GeoVaLs/air_pressure_levels");
}

// -----------------------------------------------------------------------------

TotalColumnVaporGuess::~TotalColumnVaporGuess() {}

// -----------------------------------------------------------------------------

void TotalColumnVaporGuess::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/air_pressure_levels"));
  const float GK = 1.0/Constants::grav;

  // column q (kg/m^2) = sum( pressure_thickness * (q_mixrati/(1 + q_mixrati)) / grav)
  std::vector<float> q_mixrati(nlocs);
  std::vector<float> tcwv(nlocs, 0.0);
  std::vector<float> pre_lev0(nlocs), pre_levl(nlocs);

  in.get(Variable("GeoVaLs/air_pressure_levels"), 0, pre_lev0);
  for (size_t ilev = 1; ilev < nlevs; ++ilev) {
    in.get(Variable("GeoVaLs/air_pressure_levels"), ilev, pre_levl);
    in.get(Variable("GeoVaLs/humidity_mixing_ratio"), ilev - 1, q_mixrati);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      // Change the unit of q_mixing g/kg => kg/kg.
      q_mixrati[iloc] *= 0.001;
      tcwv[iloc] = tcwv[iloc] + q_mixrati[iloc]/(q_mixrati[iloc]+1)
                   * fabs(pre_levl[iloc] - pre_lev0[iloc]);
    }
    pre_lev0 = pre_levl;
  }
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    tcwv[iloc] *= GK;
  }
  out[0] = tcwv;
}

// -----------------------------------------------------------------------------

const ufo::Variables & TotalColumnVaporGuess::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
