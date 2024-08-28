/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWind.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindParameters.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadarDopplerWind> obsRadarDopplerWindMaker_("RadarDopplerWind");
// -----------------------------------------------------------------------------

ObsRadarDopplerWind::ObsRadarDopplerWind(const ioda::ObsSpace & odb,
                                         const Parameters_ & params)
  : ObsOperatorBase(odb), odb_(odb), params_(params)
{
  const std::vector<std::string> expectedVariables{"radialVelocity"};
  if (odb_.assimvariables().variables() != expectedVariables) {
    throw eckit::UserError("This operator can only be used to simulate "
                           "radialVelocity");
  }

  const std::vector <std::string> GeoVaLnames{
      "eastward_wind",
      "northward_wind",
      "upward_air_velocity",
      params.verticalCoordinate_uv.value(),
      params.verticalCoordinate_w.value()};
  requiredVars_ += oops::Variables(GeoVaLnames);
  oops::Log::trace() << "ObsRadarDopplerWind constructed" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarDopplerWind::~ObsRadarDopplerWind() {
  oops::Log::trace() << "ObsRadarDopplerWind destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWind::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                      ObsDiagnostics &,
                                      const QCFlags_t &) const {
  // This routine does the following:
  // - Retrieves GeoVaLs of u, v, w and z.
  // - Interpolates u, v and w to observation z.
  // - Computes H(x).

  oops::Log::trace() << "ObsRadarDopplerWind: simulateObs entered" << std::endl;

  const float missingFloat = util::missingValue<float>();
  const double missingDouble = util::missingValue<double>();

  const std::size_t nlocs = odb_.nlocs();
  const std::size_t nlevs_uv = gv.nlevs(oops::Variable{"eastward_wind"});
  const std::size_t nlevs_w = gv.nlevs(oops::Variable{"upward_air_velocity"});

  // Obtain beam geometry variables.
  // These must have been precomputed before running the observation operator.
  std::vector<float> sinTilt(nlocs);
  std::vector<float> cosAzimuthCosTilt(nlocs);
  std::vector<float> sinAzimuthCosTilt(nlocs);
  std::vector<float> height(nlocs);
  odb_.get_db("MetaData", "sinTilt", sinTilt);
  odb_.get_db("MetaData", "cosAzimuthCosTilt", cosAzimuthCosTilt);
  odb_.get_db("MetaData", "sinAzimuthCosTilt", sinAzimuthCosTilt);
  odb_.get_db("MetaData", "height", height);

  // Compute radial velocity H(x) at each location:
  std::vector<double> vec_z_uv(nlevs_uv);
  std::vector<double> vec_z_w(nlevs_w);
  std::vector<double> vec_u(nlevs_uv);
  std::vector<double> vec_v(nlevs_uv);
  std::vector<double> vec_w(nlevs_w);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    const float z_ob = height[jloc];
    if (ufo::anyEqualTo(missingFloat,
                        z_ob,
                        sinAzimuthCosTilt[jloc],
                        cosAzimuthCosTilt[jloc],
                        sinTilt[jloc])) {
      ovec[jloc] = missingDouble;
      continue;
    }

    // Retrieve GeoVaLs at this location.
    gv.getAtLocation(vec_z_uv, oops::Variable{params_.verticalCoordinate_uv.value()}, jloc);
    gv.getAtLocation(vec_z_w, oops::Variable{params_.verticalCoordinate_w.value()}, jloc);
    gv.getAtLocation(vec_u, oops::Variable{"eastward_wind"}, jloc);
    gv.getAtLocation(vec_v, oops::Variable{"northward_wind"}, jloc);
    gv.getAtLocation(vec_w, oops::Variable{"upward_air_velocity"}, jloc);

    // Interpolate u, v, and w to observation height.
    const double u = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_uv, vec_u, z_ob);
    const double v = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_uv, vec_v, z_ob);
    const double w = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_w, vec_w, z_ob);

    // Combine u, v and w with beam geometry variable to produce the final H(x).
    if (ufo::anyEqualTo(missingDouble, u, v, w)) {
      ovec[jloc] = missingDouble;
    } else {
      ovec[jloc] = sinAzimuthCosTilt[jloc] * u +
        cosAzimuthCosTilt[jloc] * v +
        sinTilt[jloc] * w;
    }
  }

  oops::Log::trace() << "ObsRadarDopplerWind: simulateObs finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWind::print(std::ostream & os) const {
  os << "ObsRadarDopplerWind" << std::endl;
  os << params_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
