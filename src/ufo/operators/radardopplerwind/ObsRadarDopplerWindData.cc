/*
 * (C) Copyright 2024 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/missingValues.h"

#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindData.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsRadarDopplerWindData::ObsRadarDopplerWindData
(const ioda::ObsSpace & odb, const Parameters_ & parameters)
  : odb_(odb),
    verticalCoordinate_uv_(oops::Variable{parameters.verticalCoordinate_uv.value()}),
    verticalCoordinate_w_(oops::Variable{parameters.verticalCoordinate_w.value()})
{
  const std::vector <std::string> GeoVaLnames{
      "eastward_wind",
      "northward_wind",
      "upward_air_velocity"};
  requiredVars_ += oops::Variables(GeoVaLnames);
  requiredVars_.push_back(verticalCoordinate_uv_);
  requiredVars_.push_back(verticalCoordinate_w_);
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindData::cacheVertCoordGeoVaLs(const GeoVaLs & gv) const {
  // The GeoVaL caching must only be performed once.
  if (!vec_z_uv_.empty() || !vec_z_w_.empty())
    return;
  const std::size_t nlocs = gv.nlocs();
  const std::size_t nlevs_z_uv = gv.nlevs(verticalCoordinate_uv_);
  const std::size_t nlevs_z_w = gv.nlevs(verticalCoordinate_w_);
  std::vector<double> vec_z_uv(nlevs_z_uv);
  std::vector<double> vec_z_w(nlevs_z_w);
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    gv.getAtLocation(vec_z_uv, verticalCoordinate_uv_, jloc);
    gv.getAtLocation(vec_z_w, verticalCoordinate_w_, jloc);
    vec_z_uv_.push_back(vec_z_uv);
    vec_z_w_.push_back(vec_z_w);
  }
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindData::fillHofX(const ioda::ObsSpace & odb,
                                       const GeoVaLs & gv, ioda::ObsVector & ovec) const {
  // The vertical coordinate GeoVaLs must have been cached prior to running this function.
  ASSERT(!vec_z_uv_.empty() && !vec_z_w_.empty());

  const float missingFloat = util::missingValue<float>();
  const double missingDouble = util::missingValue<double>();

  const std::size_t nlocs = odb.nlocs();
  const std::size_t nvars = ovec.nvars();
  ASSERT(nvars == 1);  // Only one simulated variable (radialVelocity) is permitted.
  const std::size_t nlevs_uv = gv.nlevs(oops::Variable{"eastward_wind"});
  const std::size_t nlevs_w = gv.nlevs(oops::Variable{"upward_air_velocity"});

  // Obtain beam geometry variables.
  // These must have been precomputed before running the observation operator.
  std::vector<float> sinTilt(nlocs);
  std::vector<float> cosAzimuthCosTilt(nlocs);
  std::vector<float> sinAzimuthCosTilt(nlocs);
  std::vector<float> height(nlocs);
  odb.get_db("MetaData", "sinTilt", sinTilt);
  odb.get_db("MetaData", "cosAzimuthCosTilt", cosAzimuthCosTilt);
  odb.get_db("MetaData", "sinAzimuthCosTilt", sinAzimuthCosTilt);
  odb.get_db("MetaData", "height", height);

  // Determine H(x) at each location.
  std::vector<double> vec_u(nlevs_uv);
  std::vector<double> vec_v(nlevs_uv);
  std::vector<double> vec_w(nlevs_w);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    gv.getAtLocation(vec_u, oops::Variable{"eastward_wind"}, jloc);
    gv.getAtLocation(vec_v, oops::Variable{"northward_wind"}, jloc);
    gv.getAtLocation(vec_w, oops::Variable{"upward_air_velocity"}, jloc);
    const std::vector<double> & vec_z_uv = vec_z_uv_[jloc];
    const std::vector<double> & vec_z_w = vec_z_w_[jloc];

    const float z_ob = height[jloc];

    // Calculate H(x) for u, v, and w.
    const double u = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_uv, vec_u, z_ob);
    const double v = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_uv, vec_v, z_ob);
    const double w = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_w, vec_w, z_ob);

    // Combine u, v and w with beam geometry variable to produced final H(x).
    if (ufo::anyEqualTo(missingDouble, u, v, w) ||
        ufo::anyEqualTo(missingFloat, sinAzimuthCosTilt[jloc],
                        cosAzimuthCosTilt[jloc],
                        sinTilt[jloc])) {
      ovec[jloc] = missingDouble;
    } else {
      ovec[jloc] = sinAzimuthCosTilt[jloc] * u +
        cosAzimuthCosTilt[jloc] * v +
        sinTilt[jloc] * w;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
