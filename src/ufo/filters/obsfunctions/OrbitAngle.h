/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_ORBITANGLE_H_
#define UFO_FILTERS_OBSFUNCTIONS_ORBITANGLE_H_

#include <vector>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Compute the orbital angle of observations (in degrees)
/// using the emphemeris lats and longs
/// which record the position and movement of the satellite.
///
/// References:
/// * `Ops_OrbitAngle` (subroutine in the Met Office OPS system): original source code
/// * `ORBIT_ANGLE` (subroutine in AAPP)
///
class OrbitAngle : public ObsFunctionBase<float> {
 public:
  explicit OrbitAngle(const eckit::LocalConfiguration &
                                     = eckit::LocalConfiguration());

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  float vectormod(const std::vector<float> &a) const;
  float dotproduct(const std::vector<float> &a,
                  const std::vector<float> &b) const;
  std::vector<float> ll_to_xyz(const float &lat, const float &lon) const;
  std::vector<float> crossprod_xyz(const std::vector<float> &a,
                                   const std::vector<float> &b) const;

 private:
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_ORBITANGLE_H_
