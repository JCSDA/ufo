/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include "ufo/filters/obsfunctions/SatwindIndivErrors.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"

namespace ufo {

static ObsFunctionMaker<SatwindIndivErrors> makerSatwindIndivErrors_("SatwindIndivErrors");

// -----------------------------------------------------------------------------

SatwindIndivErrors::SatwindIndivErrors(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get parameters from options.
  std::string const profile = options_.profile.value();
  std::string const obs_vcoord = options_.obs_vcoord.value();
  std::string const vcoord = options_.vcoord.value();

  // Include list of required data from ObsSpace
  invars_ += Variable("MetaData/" + obs_vcoord);
  invars_ += Variable("HofX/" + profile);
  invars_ += options_.pressure_error;
  invars_ += options_.quality_index;

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/eastward_wind");
  invars_ += Variable("GeoVaLs/northward_wind");
  invars_ += Variable("GeoVaLs/" + vcoord);
}

// -----------------------------------------------------------------------------

SatwindIndivErrors::~SatwindIndivErrors() {}

// -----------------------------------------------------------------------------
/*! \brief Function to calculate situation dependent observation errors for satwinds
 *
 * \details The errors are calculated by combining an estimate of the error in
 * the vector, with an estimate of the error in vector due to an error in the
 * pressure (i.e. height assignment). The latter will depend on the AMV height
 * error and the background vertical wind shear. In the future we aim to use
 * estimates of the vector error and height error from the data producers.
 *
 * Currently the vector error estimate \f$E_{vector}\f$
 * is based on the quality index (QI) (ideally model-independent) and is calculated as:
 * \f[
 * E_{vector} = \text{EuMult}\left(QI \times 0.01\right) + \text{EuAdd}
 * \f]
 * The defaults are EuMult=-5.0 and EuAdd=7.5, which gives errors in the range
 * from 2.5 m/s (at QI=100) to 4.5 m/s (at QI=60).
 *
 * The height error estimate \f$E_{p}\f$ is currently set by a look up table
 * (dependent on e.g. satellite, channel, pressure level).
 * The values are based on the RMS of model best-fit pressure minus AMV observed
 * pressure distributions. These are calculated from several months of data.
 *
 * The error in vector due to the height error, \f$E_{vpress}\f$, is calculated as:
 * \f[
 * E_{vpress} = \sqrt{\frac{\sum{W_{i}\left(v_{i}-v_{n}\right)^{2}}}
 *                     {\sum{W_{i}}}
 *               }
 * \f]
 * where
 * \f[
 * W_{i} = \exp\left(- \left( \left( p_{i} - p_{n} \right)^2 / 2E_{p}^2 \right)\right) \times dP_{i}
 * \f]
 * i = model level \n
 * N = number of model levels (sum over) \n
 * \f$v_i\f$ = wind component on model level \n
 * \f$v_n\f$ = wind component at observation location \n
 * \f$p_i\f$ = pressure on model level \n
 * \f$p_n\f$ = pressure at observation location \n
 * \f$E_p\f$ = error in height assignment \n
 * \f$dP_i\f$ = layer thickness \n
 *
 * This is calculated separately for the u and v components giving separate
 * u and v component errors.
 * \f[
 * E_{total}^2 = E_{vector}^2 + E_{vpress}^2
 * \f]
 *
 * \author J.Cotton (Met Office)
 *
 * \date 05/01/2021: Created
 */

void SatwindIndivErrors::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & obserr) const {
  // Get obs space
  auto & obsdb = in.obsspace();

  // Get parameters from options
  float const eu_add = options_.eu_add.value();
  float const eu_mult = options_.eu_mult.value();
  float const min_press = options_.min_press.value();  // Pa
  std::string const profile = options_.profile.value();
  std::string const obs_vcoord = options_.obs_vcoord.value();
  oops::Variable const vcoord{options_.vcoord.value()};
  oops::Log::debug() << "Wind profile for calculating observation errors is "
                     << profile << std::endl;
  oops::Log::debug() << "Observation vertical coordinate is " << obs_vcoord << std::endl;
  oops::Log::debug() << "Model vertical coordinate is " << vcoord << std::endl;

  std::ostringstream errString;

  // check profile name matches one of windEastward or windNorthward
  if ( profile != "windEastward" && profile != "windNorthward" ) {
    errString << "Wind component must be one of windEastward or windNorthward" << std::endl;
    throw eckit::BadValue(errString.str(), Here());
  }

  // check vcoord name matches air_pressure, air_pressure_levels, or air_pressure_levels_minus_one
  if ( vcoord.name() != "air_pressure" && vcoord.name() != "air_pressure_levels" &&
       vcoord.name() != "air_pressure_levels_minus_one") {
    errString << "Vertical coordinate not recognised" << std::endl;
    throw eckit::BadValue(errString.str(), Here());
  }

  // Get dimensions
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/" + vcoord.name()));

  // local variables
  const float missing = util::missingValue<float>();

  // Get variables from ObsSpace if present. If not, throw an exception
  std::vector<float> ob_p;
  std::vector<float> bg_windcomponent;
  std::vector<float> pressure_error;
  std::vector<float> ob_qi;
  in.get(Variable("MetaData/" + obs_vcoord), ob_p);
  in.get(Variable("HofX/" + profile), bg_windcomponent);
  in.get(options_.pressure_error, pressure_error);
  in.get(options_.quality_index, ob_qi);

  // Get GeoVaLs
  const ufo::GeoVaLs * gvals = in.getGeoVaLs();
  // Vectors storing GeoVaL column for current location.
  std::vector <double> cx_p(nlevs, 0.0);
  std::vector <double> cx_windcomponent(nlevs, 0.0);

  // diagnostic variables to be summed over all processors at the end of the routine
  std::unique_ptr<ioda::Accumulator<size_t>> countQiAccumulator =
    obsdb.distribution()->createAccumulator<size_t>();

  // Loop through locations
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    // Get GeoVaLs at the specified location.
    gvals->getAtLocation(cx_p, vcoord, iloc);
    // Temporary mapping for variable names
    if ( profile == "windEastward" ) {
      gvals->getAtLocation(cx_windcomponent, oops::Variable{"eastward_wind"}, iloc);
    } else if ( profile == "windNorthward" ) {
      gvals->getAtLocation(cx_windcomponent, oops::Variable{"northward_wind"}, iloc);
    }
    // Check GeoVaLs are in correct vertical order
    if (cx_p.front() > cx_p.back()) {
      throw eckit::BadValue("GeoVaLs are not ordered from model top to bottom", Here());
    }
    // Initialize at each location
    float error_press = 0.0;  // default wind error contribution from error in pressure, ms-1
    float error_vector = 3.5;  // default wind error contribution from error in vector, ms-1
    double weight = 0.0;
    double sum_top = 0.0;
    double sum_weight = 0.0;
    obserr[0][iloc] = missing;

    // Check for valid pressure error estimate.
    // The choice of minimum pressure error is somewhat arbitrary.
    // Having a minimum helps catch cases where we may have forgotten to specify in Pa
    // and instead used hPa, e.g. 60 Pa instead of 6000 Pa.
    float const min_pressure_error =   500;  //   5 hPa
    float const max_pressure_error = 50000;  // 500 hPa
    if ( pressure_error[iloc] == missing ||
         pressure_error[iloc] < min_pressure_error ||
         pressure_error[iloc] > max_pressure_error ) {
        errString << "Pressure error invalid: " << pressure_error[iloc] << " Pa" << std::endl;
        throw eckit::BadValue(errString.str(), Here());
    }

    // Calculate vector error using QI
    if ( (ob_qi[iloc] != missing) &&
         (ob_qi[iloc] > 0.0) &&
         (ob_qi[iloc] <= 100.0)) {
      error_vector = eu_mult * (ob_qi[iloc] * 0.01) + eu_add;
    } else {
      countQiAccumulator->addTerm(iloc, 1);
    }

    // Calculate the error in vector due to error in pressure
    // Loop over levels starting from highest pressure (bottom to top)
    // Start at level number 2 (ilev=nlevs-2) as need to difference with ilev-1
    for (int ilev = nlevs-2; ilev >= 0; ilev--) {
      // ignore contribution above min_press in atmosphere
      if (cx_p[ilev] < min_press) {
        continue;
      }
      // Calculate weight for each background level, avoiding zero divide.
      if (pressure_error[iloc] > 0) {
        weight = exp(-0.5 * pow(cx_p[ilev] - ob_p[iloc], 2) /
                            pow(pressure_error[iloc], 2) )
                 * std::abs(cx_p[ilev] - cx_p[ilev + 1]);
      } else {
          weight = 0.0;
      }

      // ignore level if weight is very small
      double const min_weight = 0.00001;
      if (weight < min_weight) {
        continue;
      }

      sum_top += weight * pow(cx_windcomponent[ilev] - bg_windcomponent[iloc], 2);
      sum_weight += weight;
    }

    if (sum_weight > 0) {
      error_press = sqrt(sum_top / sum_weight);
    }

    obserr[0][iloc] = sqrt(pow(error_vector, 2) +
                           pow(error_press, 2) );
  }
  // sum number of bad QI values
  const std::size_t countQi = countQiAccumulator->computeResult();
  if (countQi > 0) {
    oops::Log::warning() << "Satwind Indiv Errors: " << countQi
      << " observations with bad/missing QI. Using default vector err in these cases" << std::endl;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SatwindIndivErrors::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
