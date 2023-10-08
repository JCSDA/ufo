/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/AircraftBiasCorrectionTerm.h"

#include <algorithm>
#include <cmath>
#include <valarray>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<AircraftBiasCorrectionTerm> makerObsFuncAircraftBiasCorrectionTerm\
_("AircraftBiasCorrectionTerm");

// -----------------------------------------------------------------------------

AircraftBiasCorrectionTerm::AircraftBiasCorrectionTerm(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "AircraftBiasCorrectionTerm: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);
  // We need to retrieve the aircraft bias correction coefficient
  // Draw variable from group/name supplied as an option
  const std::string coeff_grpvarname = options_.coeff_grpvarname.value();
  invars_ += Variable(coeff_grpvarname);
  // We need to retrieve the aircraft bias correction predictor
  // Draw variable from group/name supplied as an option
  const std::string predi_grpvarname = options_.predi_grpvarname.value();
  invars_ += Variable(predi_grpvarname);
  // We need to retrieve the order to which the predictor is to be raised
  // (e.g., 1. for a 1st order predictor, 2. for a squared predictor, etc)
  // Supplied as an option, drawn directly from options_ below
}

// -----------------------------------------------------------------------------

AircraftBiasCorrectionTerm::~AircraftBiasCorrectionTerm() {}

// -----------------------------------------------------------------------------

void AircraftBiasCorrectionTerm::compute(const ObsFilterData & in,
                                       ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();
  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);
  // Retrieve observation aircraft bias correction coefficient
  std::vector<float> BCcoeff;
  const std::string &coeff_grpvarname = options_.coeff_grpvarname.value();
  in.get(Variable(coeff_grpvarname), BCcoeff);
  // Retrieve observation aircraft bias correction predictor
  std::vector<float> BCpredi;
  const std::string &predi_grpvarname = options_.predi_grpvarname.value();
  in.get(Variable(predi_grpvarname), BCpredi);
  // Retrieve predictor's order from options_
  const float ORDpredi = options_.predi_order.value();

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (BCcoeff[jj] != missing && BCpredi[jj] != missing) {
      // BC term = BCcoeff[jj] * std::pow(BCpredi[jj], ORDpredi)
      out[0][jj] = BCcoeff[jj] * std::pow(BCpredi[jj], ORDpredi);
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & AircraftBiasCorrectionTerm::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
