/*
 * (C) Copyright 2026 NOAA/OAR/GSL
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ModelHeightAdjustedSpecificHumidity.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

static ObsFunctionMaker<ModelHeightAdjustedSpecificHumidity>
        makerModelHeightAdjustedSpecificHumidity_("ModelHeightAdjustedSpecificHumidity");

// -----------------------------------------------------------------------------

ModelHeightAdjustedSpecificHumidity::ModelHeightAdjustedSpecificHumidity(
        const eckit::LocalConfiguration & conf): invars_() {
  oops::Log::trace() << "ModelHeightAdjustedSpecificHumidity constructor" << std::endl;

  parameters_.validateAndDeserialize(conf);

  // Observed specific humidity
  invars_ += parameters_.observedSpecificHumidity.value();

  // Observed temperature at observation station height
  invars_ += parameters_.observedTemperature.value();

  // Adjusted temperature, i.e. observed temperature adjusted to model height
  invars_ += parameters_.adjustedTemperature.value();

  // Observed station pressure at observation station height
  invars_ += parameters_.stationPressure.value();

  // Model pressure at observed station height
  invars_ += parameters_.modelSurfacePressure.value();
}

// -----------------------------------------------------------------------------

void ModelHeightAdjustedSpecificHumidity::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ModelHeightAdjustedSpecificHumidity compute start" << std::endl;
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

  // Verify that the adjusted temperature variable exists (requires a prior Variable Assignment
  // using ModelHeightAdjustedAirTemperature)
  const Variable & adjTempVar = parameters_.adjustedTemperature.value();
  if (!in.has(adjTempVar)) {
    throw eckit::UserError("ModelHeightAdjustedSpecificHumidity: adjusted temperature variable '"
                           + adjTempVar.fullName() + "' not found. "
                           "Ensure a Variable Assignment filter using "
                           "ModelHeightAdjustedAirTemperature runs before this ObsFunction.",
                           Here());
  }

  // Read input data
  std::vector<float> qObs(nlocs);
  std::vector<float> tObs(nlocs);
  std::vector<float> tAdj(nlocs);
  std::vector<float> pStation(nlocs);
  std::vector<float> pModelSurface(nlocs);

  in.get(parameters_.observedSpecificHumidity.value(), qObs);
  in.get(parameters_.observedTemperature.value(), tObs);
  in.get(adjTempVar, tAdj);
  in.get(parameters_.stationPressure.value(), pStation);
  in.get(parameters_.modelSurfacePressure.value(), pModelSurface);

  // Compute adjusted specific humidity preserving relative humidity
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (qObs[jj] == missing || tObs[jj] == missing || tAdj[jj] == missing ||
        pModelSurface[jj] == missing) {
      out[0][jj] = missing;
      continue;
    }

    // Determine pressure at station height:
    // Use observed station pressure if available, otherwise fall back to model surface pressure
    float pAtStation = pStation[jj];
    if (pAtStation == missing) {
      pAtStation = pModelSurface[jj];
    }

    // Compute saturation specific humidity at station (T_obs, P_station)
    float esatObs = formulas::SatVaporPres_fromTemp(tObs[jj],
                        formulas::Formulation::Sonntag);
    if (esatObs == missing || esatObs <= 0.0f) {
      out[0][jj] = missing;
      continue;
    }
    float qsatObs = formulas::Qsat_From_Psat(esatObs, pAtStation,
                        formulas::Formulation::GillUKMO);

    // Compute relative humidity at station
    float rh = 0.0f;
    if (qsatObs > 0.0f) {
      rh = qObs[jj] / qsatObs;
    }
    // Clamp RH to [0, 1] to prevent unphysical values
    rh = std::max(0.0f, std::min(1.0f, rh));

    // Compute saturation specific humidity at model surface (T_adj, P_model_surface)
    float esatAdj = formulas::SatVaporPres_fromTemp(tAdj[jj],
                        formulas::Formulation::Sonntag);
    if (esatAdj == missing || esatAdj <= 0.0f) {
      out[0][jj] = missing;
      continue;
    }
    float qsatAdj = formulas::Qsat_From_Psat(esatAdj, pModelSurface[jj],
                        formulas::Formulation::GillUKMO);

    // Adjusted q preserving RH
    out[0][jj] = rh * qsatAdj;
  }

  oops::Log::trace() << "ModelHeightAdjustedSpecificHumidity compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ModelHeightAdjustedSpecificHumidity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
