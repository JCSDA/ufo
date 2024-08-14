/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <string>
#include <unordered_map>

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/MetOfficeReflectivity.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

static ReflectivityAlgorithmMaker<MetOfficeReflectivity>
makerMetOfficeReflectivity_("Met Office reflectivity");

MetOfficeReflectivity::MetOfficeReflectivity(const Parameters_ & params,
                                             const ioda::ObsSpace & odb,
                                             const int idxRefl,
                                             oops::Variables & reqvars,
                                             oops::Variables & reqvarsTL)
  : ReflectivityAlgorithmBase(params, odb, idxRefl, reqvars, reqvarsTL),
    params_(params)
{
  // Required variables for the nonlinear operator.
  reqvars += oops::Variables
    ({"air_pressure",
      "air_temperature",
      "height",
      "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
      "qrain"});

  // Required variables for the TL/AD operator.
  std::vector<std::string> rvTL{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
      "qrain"};
  reqvarsTL += oops::Variables(rvTL);
}

void MetOfficeReflectivity::simulateObsImpl(const GeoVaLs & gv,
                                            ioda::ObsVector & ovec,
                                            ObsDiagnostics &,
                                            const QCFlags_t &) const {
  // This routine does the following:
  // - Retrieves GeoVals of p, T, z, qrain and qice.
  // - Interpolates p, T, qrain and qice to observation z.
  // - Computes H(x).

  oops::Log::trace() << "MetOfficeReflectivity::simulateObsImpl starting" << std::endl;

  const std::size_t nlocs = obsdb_.nlocs();
  const std::size_t nlevs = gv.nlevs(oops::Variable{"height"});
  const float missingFloat = util::missingValue<float>();
  const float missingDouble = util::missingValue<double>();

  // Retrieve gate heights.
  // These must have been precomputed before running the observation operator.
  std::vector<float> height(nlocs);
  obsdb_.get_db("MetaData", "height", height);

  // Vectors that store GeoVaLs at each location.
  std::vector<double> vec_p(nlevs);
  std::vector<double> vec_T(nlevs);
  std::vector<double> vec_z(nlevs);
  std::vector<double> vec_qrain(nlevs);
  std::vector<double> vec_qice(nlevs);

  // Constants
  const double rd = ufo::Constants::rd;
  const double t0c = ufo::Constants::t0c;
  const double rainMultiplier = params_.rainMultiplier.value();
  const double rainExponent = params_.rainExponent.value();
  const double iceMultiplier = params_.iceMultiplier.value();
  const double iceAdditive = params_.iceAdditive.value();
  const double iceExponent = params_.iceExponent.value();
  const double lowerBound = params_.lowerBound.value();

  // Compute reflectivity H(x) at each location.
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const std::size_t idx = jloc * ovec.nvars() + idxReflectivity_;
    const float z_ob = height[jloc];

    if (z_ob == missingFloat) {
      ovec[idx] = missingDouble;
      continue;
    }

    gv.getAtLocation(vec_p, oops::Variable{"air_pressure"}, jloc);
    gv.getAtLocation(vec_T, oops::Variable{"air_temperature"}, jloc);
    gv.getAtLocation(vec_z, oops::Variable{"height"}, jloc);
    gv.getAtLocation(vec_qrain, oops::Variable{"qrain"}, jloc);
    gv.getAtLocation(vec_qice,
                     oops::Variable{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"},
                     jloc);

    // Interpolate each GeoVaL to the observation location.
    const double p = ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_p, z_ob);
    const double T = ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_T, z_ob);
    const double qrain =
      ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_qrain, z_ob);
    const double qice =
      ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_qice, z_ob);

    // Compute rain contribution to squared reflectivity.
    const double air_density = p / (rd * T);
    const double rainWaterContent = 1000.0 * air_density * qrain;
    const double rainComponent =
      rainWaterContent > 0.0 ? rainMultiplier * std::pow(rainWaterContent, rainExponent) : 0.0;

    // Compute ice contribution to squared reflectivity.
    const double gamma = iceMultiplier * (T - t0c) + iceAdditive;
    const double iceWaterContent = 1000.0 * air_density * qice;
    const double iceComponent =
      iceWaterContent > 0.0 ? std::pow(10.0, gamma) * std::pow(iceWaterContent, iceExponent) : 0.0;

    // Compute total squared reflectivity, including an optional lower bound.
    const double squaredReflectivity = rainComponent + iceComponent + lowerBound;

    // Fill reflectivity H(x).
    ovec[idx] = std::sqrt(squaredReflectivity);
  }

  oops::Log::trace() << "MetOfficeReflectivity::simulateObsImpl done" << std::endl;
}

void MetOfficeReflectivity::setTrajectoryImpl(const GeoVaLs & gv,
                                              ObsDiagnostics &,
                                              const QCFlags_t &) {
  // This routine does the following:
  // - Retrieves GeoVals of p, T, z, qrain and qice.
  // - Stores GeoVaLs of z.
  // - Stores values of p, T, qrain and qice interpolated to observation z.

  oops::Log::trace() << "MetOfficeReflectivity::setTrajectoryImpl starting" << std::endl;

  const std::size_t nlocs = obsdb_.nlocs();
  const std::size_t nlevs = gv.nlevs(oops::Variable{"height"});
  const float missingFloat = util::missingValue<float>();
  const float missingDouble = util::missingValue<double>();

  std::vector<float> height(nlocs);
  obsdb_.get_db("MetaData", "height", height);

  // Stored interpolated values which are used in the TL/AD routines.
  vec_gv_z_.resize(nlocs);
  vec_interp_T_.resize(nlocs);
  vec_interp_p_.resize(nlocs);
  vec_interp_qrain_.resize(nlocs);
  vec_interp_qice_.resize(nlocs);

  std::vector<double> vec_p(nlevs);
  std::vector<double> vec_T(nlevs);
  std::vector<double> vec_z(nlevs);
  std::vector<double> vec_qrain(nlevs);
  std::vector<double> vec_qice(nlevs);
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const float z_ob = height[jloc];
    gv.getAtLocation(vec_z, oops::Variable{"height"}, jloc);
    vec_gv_z_[jloc] = vec_z;
    if (z_ob == missingFloat) {
      vec_interp_T_[jloc] = missingDouble;
      vec_interp_p_[jloc] = missingDouble;
      vec_interp_qrain_[jloc] = missingDouble;
      vec_interp_qice_[jloc] = missingDouble;
    } else {
      gv.getAtLocation(vec_p, oops::Variable{"air_pressure"}, jloc);
      gv.getAtLocation(vec_T, oops::Variable{"air_temperature"}, jloc);
      gv.getAtLocation(vec_qrain, oops::Variable{"qrain"}, jloc);
      gv.getAtLocation(vec_qice,
         oops::Variable{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"}, jloc);
      vec_interp_T_[jloc] = ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_T, z_ob);
      vec_interp_p_[jloc] = ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_p, z_ob);
      vec_interp_qrain_[jloc] =
        ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_qrain, z_ob);
      vec_interp_qice_[jloc] =
        ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_qice, z_ob);
    }
  }
  oops::Log::trace() << "MetOfficeReflectivity::setTrajectoryImpl done" << std::endl;
}

void MetOfficeReflectivity::simulateObsTLImpl(const GeoVaLs & dx,
                                              ioda::ObsVector & dy,
                                              const QCFlags_t &) const {
  // This routine does the following:
  // - Retrieves GeoVaLs of qrain and qice increments.
  // - Uses stored GeoVaLs of z.
  // - Interpolates qrain and qice increments to observation z.
  // - Uses stored values of p, T, qrain and qice interpolated to observation z.
  // - Computes TL H(x).

  oops::Log::trace() << "MetOfficeReflectivity::simulateObsTLImpl starting" << std::endl;

  const std::size_t nlocs = obsdb_.nlocs();
  const std::size_t nlevs = dx.nlevs(oops::Variable{"qrain"});
  const float missingFloat = util::missingValue<float>();
  const float missingDouble = util::missingValue<double>();

  // Retrieve gate heights.
  // These must have been precomputed before running the observation operator.
  std::vector<float> height(nlocs);
  obsdb_.get_db("MetaData", "height", height);

  // Vectors that store increment GeoVaLs at each location.
  std::vector<double> vec_dqrain(nlevs);
  std::vector<double> vec_dqice(nlevs);

  // Constants
  const double rd = ufo::Constants::rd;
  const double t0c = ufo::Constants::t0c;
  const double rainMultiplier = params_.rainMultiplier.value();
  const double rainExponent = params_.rainExponent.value();
  const double iceMultiplier = params_.iceMultiplier.value();
  const double iceAdditive = params_.iceAdditive.value();
  const double iceExponent = params_.iceExponent.value();
  const double lowerBound = params_.lowerBound.value();

  // Compute reflectivity dH(x) at each location.
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const std::size_t idx = jloc * dy.nvars() + idxReflectivity_;
    const float z_ob = height[jloc];

    if (z_ob == missingFloat) {
      dy[idx] = missingDouble;
      continue;
    }

    // Retrieve stored GeoVaLs.
    const std::vector<double> & vec_z = vec_gv_z_[jloc];

    // Retrieve stored interpolated values.
    const double p = vec_interp_p_[jloc];
    const double T = vec_interp_T_[jloc];
    const double qrain = vec_interp_qrain_[jloc];
    const double qice = vec_interp_qice_[jloc];

    // Compute rain contribution to squared reflectivity.
    const double air_density = p / (rd * T);
    const double rainWaterContent = 1000.0 * air_density * qrain;
    const double rainComponent =
      rainWaterContent > 0.0 ? rainMultiplier * std::pow(rainWaterContent, rainExponent) : 0.0;

    // Compute ice contribution to squared reflectivity.
    const double gamma = iceMultiplier * (T - t0c) + iceAdditive;
    const double iceWaterContent = 1000.0 * air_density * qice;
    const double iceComponent =
      iceWaterContent > 0.0 ? std::pow(10.0, gamma) * std::pow(iceWaterContent, iceExponent) : 0.0;

    // Compute total squared reflectivity, including an optional lower bound.
    const double squaredReflectivity = rainComponent + iceComponent + lowerBound;

    // Retrieve vectors of model rain and ice increments.
    dx.getAtLocation(vec_dqrain, oops::Variable{"qrain"}, jloc);
    dx.getAtLocation(vec_dqice,
                     oops::Variable{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"},
                     jloc);
    // Interpolate increments to observation height.
    const double dqrain =
      ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_dqrain, z_ob);
    const double dqice =
      ufo::PiecewiseLinearInterpolation::interpolate(vec_z, vec_dqice, z_ob);

    // Fill dy.
    double dH = 0.0;
    if (rainWaterContent > 0.0) {
      dH += 0.5 * std::pow(squaredReflectivity, -0.5) * rainMultiplier * rainExponent *
        std::pow(rainWaterContent, rainExponent - 1.0) * 1000.0 * air_density * dqrain;
    }
    if (iceWaterContent > 0.0) {
      dH += 0.5 * std::pow(squaredReflectivity, -0.5) * std::pow(10.0, gamma) * iceExponent *
        std::pow(iceWaterContent, iceExponent - 1.0) * 1000.0 * air_density * dqice;
    }
    dy[idx] = dH;
  }

  oops::Log::trace() << "MetOfficeReflectivity::simulateObsTLImpl done" << std::endl;
}

void MetOfficeReflectivity::simulateObsADImpl(GeoVaLs & dx,
                                              const ioda::ObsVector & dy,
                                              const QCFlags_t &) const {
  // This routine does the following:
  // - Retrieves GeoVaLs of qrain and qice increments.
  // - Uses stored GeoVaLs of z.
  // - Determines interpolation indices and weights.
  // - Uses stored values of p, T, qrain and qice interpolated to observation z.
  // - Computes AD H(x) and updates the increment GeoVaLs.

  oops::Log::trace() << "MetOfficeReflectivity::simulateObsADImpl starting" << std::endl;

  const std::size_t nlocs = obsdb_.nlocs();
  const std::size_t nlevs = dx.nlevs(oops::Variable{"qrain"});
  const float missingFloat = util::missingValue<float>();
  const float missingDouble = util::missingValue<double>();

  // Retrieve gate heights.
  // These must have been precomputed before running the observation operator.
  std::vector<float> height(nlocs);
  obsdb_.get_db("MetaData", "height", height);

  // Vectors that store increment GeoVaLs at each location.
  std::vector<double> vec_dqrain(nlevs);
  std::vector<double> vec_dqice(nlevs);

  // Constants
  const double rd = ufo::Constants::rd;
  const double t0c = ufo::Constants::t0c;
  const double rainMultiplier = params_.rainMultiplier.value();
  const double rainExponent = params_.rainExponent.value();
  const double iceMultiplier = params_.iceMultiplier.value();
  const double iceAdditive = params_.iceAdditive.value();
  const double iceExponent = params_.iceExponent.value();
  const double lowerBound = params_.lowerBound.value();

  // Increment vectors which are used to update dx.
  std::vector<double> inc_qrain(nlevs);
  std::vector<double> inc_qice(nlevs);

  // Compute adjoint at each location.
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const std::size_t idx = jloc * dy.nvars() + idxReflectivity_;
    const float z_ob = height[jloc];

    // Do not update increments if dy is missing.
    if (dy[idx] == missingDouble || z_ob == missingFloat) {
      continue;
    }

    // Retrieve stored GeoVaLs.
    const std::vector<double> & vec_z = vec_gv_z_[jloc];

    // Retrieve stored interpolated values.
    const double p = vec_interp_p_[jloc];
    const double T = vec_interp_T_[jloc];
    const double qrain = vec_interp_qrain_[jloc];
    const double qice = vec_interp_qice_[jloc];

    // Compute rain contribution to squared reflectivity.
    const double air_density = p / (rd * T);
    const double rainWaterContent = 1000.0 * air_density * qrain;
    const double rainComponent =
      rainWaterContent > 0.0 ? rainMultiplier * std::pow(rainWaterContent, rainExponent) : 0.0;

    // Compute ice contribution to squared reflectivity.
    const double gamma = iceMultiplier * (T - t0c) + iceAdditive;
    const double iceWaterContent = 1000.0 * air_density * qice;
    const double iceComponent =
      iceWaterContent > 0.0 ? std::pow(10.0, gamma) * std::pow(iceWaterContent, iceExponent) : 0.0;

    // Compute total squared reflectivity, including an optional lower bound.
    const double squaredReflectivity = rainComponent + iceComponent + lowerBound;

    // Get existing increments.
    dx.getAtLocation(inc_qrain, oops::Variable{"qrain"}, jloc);
    dx.getAtLocation(inc_qice,
       oops::Variable{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"}, jloc);

    const auto[idx_interp, weight_interp] =
      ufo::PiecewiseLinearInterpolation::interpolationIndexAndWeight(vec_z, z_ob);

    // Multiplicative factors used to produce increments.
    const double mult_qrain = 0.5 * std::pow(squaredReflectivity, -0.5) * rainMultiplier *
      rainExponent * std::pow(rainWaterContent, rainExponent - 1.0) * 1000.0 * air_density;
    const double mult_qice = 0.5 * std::pow(squaredReflectivity, -0.5) * std::pow(10.0, gamma) *
      iceExponent * std::pow(iceWaterContent, iceExponent - 1.0) * 1000.0 * air_density;

    // Fill increments at relevant locations.
    inc_qice[idx_interp]     += dy[idx] * weight_interp * mult_qice;
    inc_qice[idx_interp + 1] += dy[idx] * (1.0 - weight_interp) * mult_qice;
    inc_qrain[idx_interp]     += dy[idx] * weight_interp * mult_qrain;
    inc_qrain[idx_interp + 1] += dy[idx] * (1.0 - weight_interp) * mult_qrain;

    // Update dx.
    dx.putAtLocation(inc_qice,
       oops::Variable{"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"}, jloc);
    dx.putAtLocation(inc_qrain, oops::Variable{"qrain"}, jloc);
  }

  oops::Log::trace() << "MetOfficeReflectivity::simulateObsADImpl done" << std::endl;
}

void MetOfficeReflectivity::printImpl(std::ostream & os) const {
  os << "Met Office reflectivity operator." << std::endl;
  os << "Parameters: " << params_ << std::endl;
}


}  // namespace ufo
