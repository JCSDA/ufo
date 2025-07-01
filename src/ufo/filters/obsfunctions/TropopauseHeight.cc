/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/obsfunctions/TropopauseHeight.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/ufo_utils.interface.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<TropopauseHeight>
 makerTropopauseHeight_("TropopauseHeight");

// -----------------------------------------------------------------------------

TropopauseHeight::TropopauseHeight
(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "TropopauseHeight constructor" << std::endl;
  // Validate and deserialize options.
  options_.validateAndDeserialize(conf);

  // GeoVaLs.
  invars_ += Variable(std::string("GeoVaLs/") + options_.model_pressure.value());
  invars_ += Variable(std::string("GeoVaLs/") + options_.model_specific_humidity.value());
  invars_ += Variable(std::string("GeoVaLs/") + options_.model_temperature.value());
}

// -----------------------------------------------------------------------------

TropopauseHeight::~TropopauseHeight() {
  oops::Log::trace() << "TropopauseHeight destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void TropopauseHeight::compute(const ObsFilterData & in,
                                        ioda::ObsDataVector<int> & out) const {
  oops::Log::trace() << "TropopauseHeight compute start" << std::endl;

  // Missing int value.
  const int missing = util::missingValue<int>();

  // Number of locations.
  const size_t nlocs = in.nlocs();

  // GeoVaLs.
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Number of model levels.
  const int nlevs = gv->nlevs(oops::Variable{options_.model_specific_humidity.value()});

  // Vectors of GeoVaLs.
  std::vector<float> gv_p(nlevs);
  std::vector<double> gv_q(nlevs);
  std::vector<float> gv_t(nlevs);

  // Log of model pressure.
  std::vector<double> gv_log_p(nlevs);

  // Log of model temperature.
  std::vector<double> gv_log_t(nlevs);

  // Convert the tropopause lapse rate from dT/dz to d(lnT)/d(lnp).
  const double trop_lapse_rate_log_coords =
      (Constants::trop_lapse_rate * Constants::rd) / Constants::grav;

  // Saturation specific humidity.
  std::vector<float> gv_qsat(nlevs);

  // Relative humidity calculated from q, T and p.
  std::vector<double> rh_from_q(nlevs);


  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Retrieve GeoVaLs at this location.
    gv->getAtLocation(gv_p, oops::Variable{options_.model_pressure.value()}, jloc);
    gv->getAtLocation(gv_q, oops::Variable{options_.model_specific_humidity.value()}, jloc);
    gv->getAtLocation(gv_t, oops::Variable{options_.model_temperature.value()}, jloc);

    // Initially set the tropopause level to a missing value.
    out[0][jloc] = missing;

    // If the user has not requested the relative humidity method (default is false) then the lapse
    // rate method will be used.
    if (!options_.rhmethod) {
      // Log(model pressure).
      std::transform(gv_p.cbegin(), gv_p.cend(), gv_log_p.begin(),
                    [](double p) -> double {return std::log(p);});

      // Log(model temperature).
      std::transform(gv_t.cbegin(), gv_t.cend(), gv_log_t.begin(),
                    [](double t) -> double {return std::log(t);});

      // Calculate bottom of atmosphere layer lapse rates.
      double lapse_upper_bound = (gv_log_t[nlevs-3] - gv_log_t[nlevs-2]) /
                                 (gv_log_p[nlevs-3] - gv_log_p[nlevs-2]);
      double lapse_lower_bound = (gv_log_t[nlevs-2] - gv_log_t[nlevs-1]) /
                                 (gv_log_p[nlevs-2] - gv_log_p[nlevs-1]);

      // Proceed through the atmosphere from the bottom upwards.
      for (int ilev = nlevs-2; ilev >= 2; --ilev) {
        // Calculate lapse rates
        const double lapse_previous = lapse_lower_bound;
        lapse_lower_bound = lapse_upper_bound;
        lapse_upper_bound = (gv_log_t[ilev -2] - gv_log_t[ilev-1]) /
                            (gv_log_p[ilev -2] - gv_log_p[ilev-1]);

        // Check to see if the tropopause has been found.
        if (lapse_upper_bound  < trop_lapse_rate_log_coords &&
          lapse_lower_bound  < trop_lapse_rate_log_coords &&
          lapse_previous > 0.0 &&
          gv_p[ilev] < options_.trop_max_pressure.value() &&
          gv_p[ilev] > options_.trop_min_pressure.value()) {
            // If so, set the tropopause level and exit.
            out[0][jloc] = ilev;
            break;
        }
      }
    }

    // If the tropopause height hasn't been found by the lapse rate method, the relative humidity
    // method will be called. It will also be called if the user requests it by setting rhmethod
    // to True.
    if (out[0][jloc] == missing || options_.rhmethod) {
      // Determine saturation specific humidity at each model level.
      ufo_ops_qsat_f90(gv_qsat.data(), gv_t.data(), gv_p.data(), nlevs);

      // Compute relative humidity from specific humidity.
      std::transform(gv_q.cbegin(), gv_q.cend(), gv_qsat.cbegin(), rh_from_q.begin(),
                    [](double q, double qsat) -> double {return q / qsat;});

      // Loop down through the atmosphere until the humidity rises to a critical value. There are
      // three exit methods; either the RH rises to a critical value, or the loop descends to too
      // high a pressure, or the loop reaches the surface.
      int ilev = 0;
      while (ilev < nlevs - 1 &&
            rh_from_q[ilev] < options_.trop_relative_humidity.value() &&
            gv_p[ilev] < options_.trop_max_pressure.value()) {
              ilev++;
      }
      // Set the tropopause level.
      out[0][jloc] = ilev;
    }
  }

  oops::Log::trace() << "TropopauseHeight compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & TropopauseHeight::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
