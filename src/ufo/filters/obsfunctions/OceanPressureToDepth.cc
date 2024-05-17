/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/OceanPressureToDepth.h"

#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<OceanPressureToDepth>
         makerOceanPressureToDepth_("OceanPressureToDepth");

// -----------------------------------------------------------------------------

OceanPressureToDepth::OceanPressureToDepth(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  invars_ += options_.pressure.value();
  invars_ += Variable("MetaData/latitude");
}

// -----------------------------------------------------------------------------

OceanPressureToDepth::~OceanPressureToDepth() {}

// -----------------------------------------------------------------------------

void OceanPressureToDepth::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // parameters (see eqn 3 of https://archimer.ifremer.fr/doc/00447/55889/57949.pdf)
  const float missingFloat = util::missingValue<float>();
  const double param1 = 9.780318;
  const double param2 = 5.2788E-3;
  const double param3 = 2.36E-5;
  const double param4 = 1.092E-6;
  const double param5 = -1.82E-15;
  const double param6 = 2.279E-10;
  const double param7 = 2.2512E-5;
  const double param8 = 9.72659;

  // input pressure
  ioda::ObsDataVector<float> pressure(in.obsspace(), invars_[0].toOopsObsVariables());
  in.get(invars_[0], pressure);
  // Get latitude
  std::vector<float> lats(nlocs);
  in.get(Variable("MetaData/latitude"), lats);

  // output depth
  std::vector<float> &depth = out[0];
  std::fill(depth.begin(), depth.end(), missingFloat);

  // compute depth as function of pressure and latitude
  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (pressure[0][loc] != missingFloat) {
      const float p_db = pressure[0][loc] / 10000.0;  // convert Pa to decibars
      const double sinsq_latitude = pow(std::sin(lats[loc] * Constants::deg2rad), 2);
      // latitude-dependent gravitational constant:
      const double g = param1 * (1.0 + (param2 + param3 * sinsq_latitude) * sinsq_latitude) +
                  param4 * p_db;
      depth[loc] = ((((param5 * p_db + param6) * p_db - param7) * p_db + param8) * p_db) / g;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & OceanPressureToDepth::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
