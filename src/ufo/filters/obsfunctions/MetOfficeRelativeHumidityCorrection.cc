/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/obsfunctions/MetOfficeRelativeHumidityCorrection.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/utils/ufo_utils.interface.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<MetOfficeRelativeHumidityCorrection>
 makerMetOfficeRelativeHumidityCorrection_("MetOfficeRelativeHumidityCorrection");

// -----------------------------------------------------------------------------

MetOfficeRelativeHumidityCorrection::MetOfficeRelativeHumidityCorrection
(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "MetOfficeRelativeHumidityCorrection constructor" << std::endl;
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  // Create the scaling factor depending on whether observation RH is in percentage or fraction.
  // If in percentage, the scaling factor is 100.0, otherwise it is 1.0.
  observationRelativeHumidityAtSaturation_ = (
    options_.output_relative_humidity_units.value() ==
    RelativeHumidityUnits::PERCENTAGE) ? 100.0f : 1.0f;

  // Observed air pressure
  invars_ += Variable(options_.observed_pressure.value());

  // GeoVaLs
  invars_ += Variable(std::string("GeoVaLs/air_pressure"));
  invars_ += Variable(std::string("GeoVaLs/water_vapor_mixing_ratio_wrt_moist_air"));
  invars_ += Variable(std::string("GeoVaLs/relative_humidity"));
  invars_ += Variable(std::string("GeoVaLs/air_temperature"));
}

// -----------------------------------------------------------------------------

MetOfficeRelativeHumidityCorrection::~MetOfficeRelativeHumidityCorrection() {
  oops::Log::trace() << "MetOfficeRelativeHumidityCorrection destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void MetOfficeRelativeHumidityCorrection::compute(const ObsFilterData & in,
                                        ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "MetOfficeRelativeHumidityCorrection compute start" << std::endl;

  // Missing float value.
  const float missing = util::missingValue<float>();

  // ObsSpace.
  const ioda::ObsSpace & obsdb = in.obsspace();

  // Number of locations.
  const size_t nlocs = obsdb.nlocs();

  // Observed air pressure
  std::vector<float> obs_p(nlocs);
  in.get(Variable(options_.observed_pressure.value()), obs_p);

  // GeoVaLs.
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Number of model levels.
  const int nlevs = gv->nlevs(oops::Variable{"water_vapor_mixing_ratio_wrt_moist_air"});

  // Vectors of GeoVaLs.
  std::vector<double> gv_rh(nlevs);
  std::vector<float> gv_p(nlevs);
  std::vector<float> gv_q(nlevs);
  std::vector<float> gv_t(nlevs);

  // Log of model pressure.
  std::vector<double> gv_log_p(nlevs);

  // Saturation mixing ratio.
  std::vector<float> gv_qsat(nlevs);

  // Relative humidity calculated from q, T and p.
  std::vector<float> rh_from_q(nlevs);

  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Do not proceed if observed pressure is missing.
    if (obs_p[jloc] == missing) {
      out[0][jloc] = missing;
      continue;
    }

    // Log(observed pressure).
    const double obs_log_p = std::log(obs_p[jloc]);

    // Retrieve GeoVaLs at this location.
    gv->getAtLocation(gv_p, oops::Variable{"air_pressure"}, jloc);
    gv->getAtLocation(gv_rh, oops::Variable{"relative_humidity"}, jloc);
    gv->getAtLocation(gv_q, oops::Variable{"water_vapor_mixing_ratio_wrt_moist_air"}, jloc);
    gv->getAtLocation(gv_t, oops::Variable{"air_temperature"}, jloc);

    // Log(model pressure).
    std::transform(gv_p.cbegin(), gv_p.cend(), gv_log_p.begin(),
                   [](double p) -> double {return std::log(p);});

    // Interpolation utility for rh GeoVaL.
    ufo::PiecewiseLinearInterpolation interpolate_rh(gv_log_p, gv_rh);

    // Interpolated value of rh GeoVaL.
    // (This is the value that is obtained with the VertInterp operator.)
    float interp_rh = interpolate_rh(obs_log_p);

    // Determine saturation mixing ratio at each model level.
    ufo_ops_qsat_f90(gv_qsat.data(), gv_t.data(), gv_p.data(), nlevs);

    // Compute relative humidity from specific humidity and saturation mixing ratio.
    std::transform(gv_q.cbegin(), gv_q.cend(), gv_qsat.cbegin(), rh_from_q.begin(),
                   [](double q, double qsat) -> double {
                    return q / qsat;});

    // Interpolation utility for calculated rh.
    ufo::PiecewiseLinearInterpolation interpolate_rh_from_q
      (gv_log_p, std::vector<double>(rh_from_q.begin(), rh_from_q.end()));

    // Interpolate calculated rh.
    float interp_rh_from_q = interpolate_rh_from_q(obs_log_p);

    // Cap supersaturated relative humidity if required.
    if (options_.capsupersat) {
      if (interp_rh > 1.0f)
        interp_rh = 1.0f;
      if (interp_rh_from_q > 1.0f)
        interp_rh_from_q = 1.0f;
    }

    // Output correction factor.
    out[0][jloc] = (interp_rh - interp_rh_from_q) * observationRelativeHumidityAtSaturation_;
  }

  oops::Log::trace() << "MetOfficeRelativeHumidityCorrection compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & MetOfficeRelativeHumidityCorrection::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
