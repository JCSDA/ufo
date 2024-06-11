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
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())), odb_(odb),
    data_(odb, params)
{
  oops::Log::trace() << "ObsRadarDopplerWindTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadarDopplerWindTLAD::~ObsRadarDopplerWindTLAD() {
  oops::Log::trace() << "ObsRadarDopplerWindTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::setTrajectory(const GeoVaLs & geovals,
                                            ObsDiagnostics & ydiags,
                                            const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsRadarDopplerWindTLAD: setTrajectory entered" << std::endl;

  data_.cacheVertCoordGeoVaLs(geovals);

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: setTrajectory finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::simulateObsTL(const GeoVaLs & dx, ioda::ObsVector & dy,
                                            const QCFlags_t &) const {
  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsTL entered" << std::endl;

  data_.fillHofX(odb_, dx, dy);

  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsTL finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadarDopplerWindTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy,
                                           const QCFlags_t &) const {
  oops::Log::trace() << "ObsRadarDopplerWindTLAD: simulateObsAD entered" << std::endl;

  // The vertical coordinate GeoVaLs must have been cached prior to running this function.
  ASSERT(!data_.heightGeoVaLs_uv().empty() &&
         !data_.heightGeoVaLs_w().empty());

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

  // Determine adjoint at each location.
  std::vector<double> inc_u(nlevs_uv);
  std::vector<double> inc_v(nlevs_uv);
  std::vector<double> inc_w(nlevs_w);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Do not update increments if any of the required variables are missing.
    if (dy[jloc] == missingDouble ||
        ufo::anyEqualTo(missingFloat,
                        sinAzimuthCosTilt[jloc],
                        cosAzimuthCosTilt[jloc],
                        sinTilt[jloc])) {
      continue;
    }

    // Get existing increments.
    dx.getAtLocation(inc_u, oops::Variable{"eastward_wind"}, jloc);
    dx.getAtLocation(inc_v, oops::Variable{"northward_wind"}, jloc);
    dx.getAtLocation(inc_w, oops::Variable{"upward_air_velocity"}, jloc);

    // Update increment components.
    const std::vector<double> & vec_z_uv = data_.heightGeoVaLs_uv()[jloc];
    const std::vector<double> & vec_z_w = data_.heightGeoVaLs_w()[jloc];

    const float z_ob = height[jloc];

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
  data_.print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
