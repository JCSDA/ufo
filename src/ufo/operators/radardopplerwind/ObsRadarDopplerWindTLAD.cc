/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindTLAD.h"

#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindParameters.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadarDopplerWindTLAD> makerRadarDopplerWindTL_("RadarDopplerWind");
// -----------------------------------------------------------------------------

ObsRadarDopplerWindTLAD::ObsRadarDopplerWindTLAD(const ioda::ObsSpace & odb,
                                                 const Parameters_ & params)
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())),
    odb_(odb), params_(params)
{
  const std::vector<std::string> expectedVariables{"radialVelocity"};
  if (odb_.assimvariables().variables() != expectedVariables) {
    throw eckit::UserError("This operator can only be used to simulate "
                           "radialVelocity");
  }

  const std::vector <std::string> GeoVaLnames{
      "eastward_wind",
      "northward_wind",
      "upward_air_velocity"};
  requiredVars_ += oops::Variables(GeoVaLnames);

  oops::Log::trace() << "ObsRadarDopplerWindTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarDopplerWindTLAD::~ObsRadarDopplerWindTLAD() {
  oops::Log::trace() << "ObsRadarDopplerWindTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::setTrajectory(const GeoVaLs & gv,
                                            ObsDiagnostics &,
                                            const QCFlags_t &) {
  // This routine does the following:
  // - Retrieves and stores GeoVaLs of z for use in the TL/AD routines.

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: setTrajectory entered" << std::endl;

  const std::size_t nlocs = odb_.nlocs();

  vec_z_uv_.resize(nlocs);
  vec_z_w_.resize(nlocs);

  const std::size_t nlevs_uv = gv.nlevs(oops::Variable{"eastward_wind"});
  const std::size_t nlevs_w = gv.nlevs(oops::Variable{"upward_air_velocity"});

  std::vector<double> vec_z_uv(nlevs_uv);
  std::vector<double> vec_z_w(nlevs_w);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    gv.getAtLocation(vec_z_uv, oops::Variable{params_.verticalCoordinate_uv.value()}, jloc);
    gv.getAtLocation(vec_z_w, oops::Variable{params_.verticalCoordinate_w.value()}, jloc);
    vec_z_uv_[jloc] = vec_z_uv;
    vec_z_w_[jloc] = vec_z_w;
  }

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: setTrajectory finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::simulateObsTL(const GeoVaLs & dx,
                                            ioda::ObsVector & dy,
                                            const QCFlags_t &) const {
  // This routine does the following:
  // - Retrieves GeoVaLs of u, v and w increments.
  // - Uses stored GeoVaLs of z.
  // - Interpolates u, v and w increments to observation z.
  // - Computes tangent linear H(x) (linearized about z).

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsTL entered" << std::endl;

  const float missingFloat = util::missingValue<float>();
  const double missingDouble = util::missingValue<double>();

  const std::size_t nlocs = odb_.nlocs();
  const std::size_t nlevs_uv = dx.nlevs(oops::Variable{"eastward_wind"});
  const std::size_t nlevs_w = dx.nlevs(oops::Variable{"upward_air_velocity"});

  // Obtain beam geometry variables.
  // These must have been precomputed before running the observation operator.
  std::vector<float> height(nlocs);
  std::vector<float> sinTilt(nlocs);
  std::vector<float> cosAzimuthCosTilt(nlocs);
  std::vector<float> sinAzimuthCosTilt(nlocs);
  odb_.get_db("MetaData", "height", height);
  odb_.get_db("MetaData", "sinTilt", sinTilt);
  odb_.get_db("MetaData", "cosAzimuthCosTilt", cosAzimuthCosTilt);
  odb_.get_db("MetaData", "sinAzimuthCosTilt", sinAzimuthCosTilt);

  // Compute tangent linear operator at each location:
  std::vector<double> inc_u(nlevs_uv);
  std::vector<double> inc_v(nlevs_uv);
  std::vector<double> inc_w(nlevs_w);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    const float z_ob = height[jloc];
    if (ufo::anyEqualTo(missingFloat,
                        z_ob,
                        sinAzimuthCosTilt[jloc],
                        cosAzimuthCosTilt[jloc],
                        sinTilt[jloc])) {
      dy[jloc] = missingDouble;
      continue;
    }

    // Get stored GeoVaL heights.
    const std::vector<double> & vec_z_uv = vec_z_uv_[jloc];
    const std::vector<double> & vec_z_w = vec_z_w_[jloc];

    // Get GeoVaL increments.
    dx.getAtLocation(inc_u, oops::Variable{"eastward_wind"}, jloc);
    dx.getAtLocation(inc_v, oops::Variable{"northward_wind"}, jloc);
    dx.getAtLocation(inc_w, oops::Variable{"upward_air_velocity"}, jloc);

    // Interpolate increments to observation heights.
    const double du = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_uv, inc_u, z_ob);
    const double dv = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_uv, inc_v, z_ob);
    const double dw = ufo::PiecewiseLinearInterpolation::interpolate(vec_z_w, inc_w, z_ob);

    // Combine u, v and w with beam geometry variable to produced final dH.
    if (ufo::anyEqualTo(missingDouble, du, dv, dw)) {
      dy[jloc] = missingDouble;
    } else {
      dy[jloc] = sinAzimuthCosTilt[jloc] * du +
        cosAzimuthCosTilt[jloc] * dv +
        sinTilt[jloc] * dw;
    }
  }

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsTL finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy,
                                           const QCFlags_t &) const {
  // This routine does the following:
  // - Retrieves GeoVaLs of u, v and w increments.
  // - Uses stored GeoVaLs of z.
  // - Determines interpolation indices and weights.
  // - Computes adjoint H(x) (linearized about z) and updates the increment GeoVaLs.

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsAD entered" << std::endl;

  const float missingFloat = util::missingValue<float>();
  const double missingDouble = util::missingValue<double>();

  const std::size_t nlocs = odb_.nlocs();
  const std::size_t nlevs_uv = dx.nlevs(oops::Variable{"eastward_wind"});
  const std::size_t nlevs_w = dx.nlevs(oops::Variable{"upward_air_velocity"});

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

  // Compute linear operator adjoint at each location:
  std::vector<double> inc_u(nlevs_uv);
  std::vector<double> inc_v(nlevs_uv);
  std::vector<double> inc_w(nlevs_w);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    const float z_ob = height[jloc];

    // Do not update increments if any of the required variables are missing.
    if (dy[jloc] == missingDouble ||
        ufo::anyEqualTo(missingFloat,
                        z_ob,
                        sinAzimuthCosTilt[jloc],
                        cosAzimuthCosTilt[jloc],
                        sinTilt[jloc])) {
      continue;
    }

    // Get stored GeoVaL heights.
    const std::vector<double> & vec_z_uv = vec_z_uv_[jloc];
    const std::vector<double> & vec_z_w = vec_z_w_[jloc];

    // Get GeoVaL increments.
    dx.getAtLocation(inc_u, oops::Variable{"eastward_wind"}, jloc);
    dx.getAtLocation(inc_v, oops::Variable{"northward_wind"}, jloc);
    dx.getAtLocation(inc_w, oops::Variable{"upward_air_velocity"}, jloc);

    // Update increment components.
    const auto[idx_uv, weight_uv] =
      ufo::PiecewiseLinearInterpolation::interpolationIndexAndWeight(vec_z_uv, z_ob);
    inc_u[idx_uv]     += dy[jloc] * sinAzimuthCosTilt[jloc] * weight_uv;
    inc_u[idx_uv + 1] += dy[jloc] * sinAzimuthCosTilt[jloc] * (1.0 - weight_uv);
    inc_v[idx_uv]     += dy[jloc] * cosAzimuthCosTilt[jloc] * weight_uv;
    inc_v[idx_uv + 1] += dy[jloc] * cosAzimuthCosTilt[jloc] * (1.0 - weight_uv);

    const auto[idx_w, weight_w] =
      ufo::PiecewiseLinearInterpolation::interpolationIndexAndWeight(vec_z_w, z_ob);
    inc_w[idx_w]     += dy[jloc] * sinTilt[jloc] * weight_w;
    inc_w[idx_w + 1] += dy[jloc] * sinTilt[jloc] * (1.0 - weight_w);

    // Store new increments.
    dx.putAtLocation(inc_u, oops::Variable{"eastward_wind"}, jloc);
    dx.putAtLocation(inc_v, oops::Variable{"northward_wind"}, jloc);
    dx.putAtLocation(inc_w, oops::Variable{"upward_air_velocity"}, jloc);
  }

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsAD finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::print(std::ostream & os) const {
  os << "ObsRadarDopplerWindTLAD" << std::endl;
  os << params_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
