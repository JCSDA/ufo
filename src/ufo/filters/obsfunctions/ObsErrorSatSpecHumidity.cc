/*
 * (C) Copyright 2023 NASA GMAO
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorSatSpecHumidity.h"

#include <float.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/PropertiesOfNVectors.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/Constants.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorSatSpecHumidity> maker_("ObsErrorSatSpecHumidity");

// -----------------------------------------------------------------------------

ObsErrorSatSpecHumidity::ObsErrorSatSpecHumidity(const eckit::Configuration &config)
  : invars_() {
  // Initialize options
  options_.reset(new ObsErrorSatSpecHumidityParameters());
  options_->deserialize(config);

  // Include list of required data for the estimation of saturation specific humidity
  invars_ += Variable("MetaData/pressure");
  invars_ += Variable("GeoVaLs/air_pressure");
  invars_ += Variable("GeoVaLs/saturation_specific_humidity");
}

// -----------------------------------------------------------------------------

ObsErrorSatSpecHumidity::~ObsErrorSatSpecHumidity() {
    oops::Log::debug() << "ObsErrorSatSpecHumidity: destructing "  << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorSatSpecHumidity::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue<float>();
  float logp_ob, satSpecificHumidity;
  double d_err;

  // Get dimensions
  const size_t nlocs = data.nlocs();
  const size_t nlevs = data.nlevs(Variable("GeoVaLs/air_pressure"));

  const std::string inflatevars = options_->inflatevars.value();
  const std::string errgrp = options_->testObserr.value();
  const std::string inputerr_name = options_->inputerr_name.value();
  std::vector<float> error(nlocs);
  data.get(Variable(errgrp+"/"+inflatevars), error);

  std::vector<int> itype(nlocs);
  data.get(Variable("ObsType/"+inflatevars), itype);

  // Get obs pressure
  std::vector<float> ob_pressure(nlocs);
  data.get(Variable("MetaData/pressure"), ob_pressure);

  std::vector<float> inputErr;
  data.get(Variable(inputerr_name+"/"+inflatevars), inputErr);

  // Get GeoVaLs
  const ufo::GeoVaLs * gvals = data.getGeoVaLs();
  // Vectors storing GeoVaL column for one obs location.
  std::vector<double> airTemperature_gval(nlevs, 0.0);
  std::vector<double> pressure_gval(nlevs, 0.0), logp(nlevs);
  std::vector<double> q_profile(nlevs, 0.0);

  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    gvals->getAtLocation(pressure_gval, oops::Variable{"air_pressure"}, jobs);
    gvals->getAtLocation(q_profile, oops::Variable{"saturation_specific_humidity"}, jobs);

    // Convert pressure to log(pressure)
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      logp[ilev] = log(pressure_gval[ilev]);
    }
    logp_ob = log(ob_pressure[jobs]);

    ufo::PiecewiseLinearInterpolation vert_interp_model(logp, q_profile);
    if ((itype[jobs] >= 180) && (itype[jobs] <= 184)) {
      satSpecificHumidity = q_profile[0];
    } else {
      satSpecificHumidity = vert_interp_model(logp_ob);
    }

    // Output error is the product of interpolated RH error and forecast
    // saturation specific humidity
    d_err = static_cast<double>(inputErr[jobs]);
    obserr[0][jobs] = d_err * satSpecificHumidity;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorSatSpecHumidity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
