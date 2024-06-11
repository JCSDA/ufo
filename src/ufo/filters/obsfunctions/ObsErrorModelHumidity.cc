/*
 * (C) Copyright 2022 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorModelHumidity.h"

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

static ObsFunctionMaker<ObsErrorModelHumidity> maker_("ObsErrorModelHumidity");

// -----------------------------------------------------------------------------

ObsErrorModelHumidity::ObsErrorModelHumidity(const eckit::Configuration &config)
  : invars_() {
  oops::Log::debug() << "ObsErrorModelHumidity: config = " << config << std::endl;

  // Initialize options
  options_.deserialize(config);

  // Get piece-wise parameters from options.
  const std::vector<float> &xvals = options_.xvals.value();
  const std::vector<float> &errors = options_.errors.value();

  // Ensure same size vectors (xvals and errors)
  if (!oops::allVectorsSameNonZeroSize(xvals, errors)) {
      std::ostringstream errString;
      errString << "At least one vector is the wrong size "
                << std::endl << "Vector sizes of xvals, errors: "
                << oops::listOfVectorSizes(xvals, errors) << std::endl;
      oops::Log::error() << errString.str();
      throw eckit::BadValue(errString.str());
  }

  // Initialize x-variable
  const Variable &xvar = options_.xvar.value();
  ASSERT(xvar.size() == 1);
  invars_ += xvar;

  // Check that all errors >= 0
  for (size_t i = 0; i < errors.size(); ++i) {
    ASSERT(errors[i] >= 0.0);
  }

  if (errors.size() > 1) {
    // In order to check for beyond range of xvals, determine if ascending or descending.
    // Also ensure that the entire vector is consistent throughout and no consecutive elements
    // of xval are equal.
    if (xvals.back() < xvals.front()) {
      isAscending_ = false;
    }
    for (size_t kv = 1; kv < xvals.size(); ++kv) {
      if ((xvals[kv] >= xvals[kv-1] && !isAscending_) ||
          (xvals[kv] <= xvals[kv-1] && isAscending_)) {
        std::ostringstream errString;
        errString << "The xvals vector of elements is NOT internally consistent."
                  << " It must be entirely either ascending or descending order,"
                  << " or two consecutive values are the same." << std::endl;
        oops::Log::error() << errString.str();
        throw eckit::BadValue(errString.str());
      }
    }
  }
  oops::Log::debug() << "ObsErrorModelHumidity: config (constructor) = "
                     << config << std::endl;

  // Include list of required data for the estimation of saturation specific humidity
  invars_ += Variable("MetaData/pressure");
  invars_ += Variable("GeoVaLs/air_pressure");
  invars_ += Variable("GeoVaLs/air_temperature");
}

// -----------------------------------------------------------------------------

ObsErrorModelHumidity::~ObsErrorModelHumidity() {
    oops::Log::debug() << "ObsErrorModelHumidity: destructing "  << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorModelHumidity::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue<float>();
  float logp_ob, temp, satVaporPres, satSpecificHumidity;

  // Get dimensions
  const size_t nlocs = data.nlocs();
  const size_t nlevs = data.nlevs(Variable("GeoVaLs/air_pressure"));

  // Get obs pressure
  std::vector<float> ob_pressure(nlocs);
  data.get(Variable("MetaData/pressure"), ob_pressure);

  // Get GeoVaLs
  const ufo::GeoVaLs * gvals = data.getGeoVaLs();
  // Vectors storing GeoVaL column for one obs location.
  std::vector<double> airTemperature_gval(nlevs, 0.0);
  std::vector<double> pressure_gval(nlevs, 0.0), logp(nlevs);

  // Get formulation for saturated vapor pressure
  formulas::MethodFormulation formulation;
  formulation = formulas::resolveFormulations(options_.Formulation, options_.Method);

  // Linearly interpolate from y0 to y1 at xstar between x0 and x1 to arrive at error
  float x0, x1, y0, y1;
  float xstar, error;

  // Get the x-variable name and piece-wise parameters from options
  const Variable &xvar = options_.xvar.value();
  const std::vector<float> &xvals = options_.xvals.value();
  const std::vector<float> &errors = options_.errors.value();
  oops::Log::debug() << "  ObsErrorModelHumidity, x-variable name: " << xvar.variable()
                     << "  and group: " << xvar.group() << std::endl;

  // Populate the testdata array.  xstar is just the 0..nloc-1 value of testvar[iv]
  // At each nloc, find matching (x0,x1) and (y0,y1) pair for linear interp.
  ioda::ObsDataVector<float> testdata(data.obsspace(), xvar.toOopsObsVariables());
  data.get(xvar, testdata);

  // The 1st index of data should have size 1 and 2nd index should be size nlocs.
  int iv = 0;
  if (testdata[iv].size() != obserr[iv].size()) {
    std::ostringstream errString;
    errString << "Something is wrong, testdata size not equal obserr size."
              << " Sizes: " << testdata[iv].size() << " and " << obserr[iv].size() << std::endl;
    oops::Log::error() << errString.str();
    throw eckit::BadValue(errString.str());
  }

  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    error = 0.0;
    obserr[iv][jobs] = missing;
    if (ob_pressure[jobs] == missing) continue;

    // Step 1: get input RH error
    if (errors.size() == 1) {
      error = errors[0];
    } else {
      if (testdata[iv][jobs] == missing) continue;
      // interpolate input RH errors to obs pressure (from ObsErrorModelStepwiseLinear)
      xstar = testdata[iv][jobs];
      if ((xstar <= xvals[0] && isAscending_) || (xstar >= xvals[0] && !isAscending_)) {
        error = errors[0];
      } else if ((xstar >= xvals.back() && isAscending_)
                || (xstar <= xvals.back() && !isAscending_)) {
        error = errors[errors.size()-1];
      } else {
        for (size_t kv = 1; kv < xvals.size(); ++kv) {
          if ((isAscending_ && (xstar > xvals[kv-1]) && (xstar <= xvals[kv])) ||
                (!isAscending_ && (xstar < xvals[kv-1]) && (xstar >= xvals[kv]))) {
            x0 = xvals[kv-1];
            x1 = xvals[kv];
            y0 = errors[kv-1];
            y1 = errors[kv];
            error = y0 + (xstar-x0)*((y1-y0)/(x1-x0));
            break;
          }
        }
      }
    }

    // Step 2: derive forecast saturation specific humidity
    // Get GeoVaLs at current obs location
    gvals->getAtLocation(airTemperature_gval, oops::Variable{"air_temperature"}, jobs);
    gvals->getAtLocation(pressure_gval, oops::Variable{"air_pressure"}, jobs);

    // Convert pressure to log(pressure)
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      logp[ilev] = log(pressure_gval[ilev]);
    }
    logp_ob = log(ob_pressure[jobs]);

    // Interpolate geoval temperature to obs pressure
    ufo::PiecewiseLinearInterpolation vert_interp_model(logp, airTemperature_gval);
    temp = vert_interp_model(logp_ob);

    // Calculate saturated vapor pressure using selected method/formulation
    satVaporPres = formulas::SatVaporPres_fromTemp(temp, formulation);

    // Convert saturated vapor pressure to saturation specific humidity
    satSpecificHumidity = Constants::epsilon * satVaporPres / ob_pressure[jobs];

    // Output error is the product of interpolated RH error and forecast
    // saturation specific humidity
    obserr[iv][jobs] = error * std::max(1.0e-12f, satSpecificHumidity);
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorModelHumidity::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
