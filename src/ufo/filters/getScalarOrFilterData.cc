/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/getScalarOrFilterData.h"

#include <sstream>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

// -----------------------------------------------------------------------------

std::vector<float> getScalarOrFilterData(const std::string & strfactor,
                                         const ObsFilterData & data) {
  std::istringstream iss(strfactor);
  std::vector<float> factors(data.nlocs());
  float factor;
  iss >> factor;
// Check float was read:
  if (iss.eof() && !iss.fail()) {
//  single float in the config:
    oops::Log::debug() << "processing a float: " << factor << std::endl;
    std::fill(factors.begin(), factors.end(), factor);
  } else {
//  it's a string; get from ObsFilterData
    Variable var(strfactor);
    oops::Log::debug() << "processing data: " << var << std::endl;
    if (!data.has(var)) {
      oops::Log::error() << "getScalarOrFilterData: either a value or a valid variable from "
                         << "data available to filter should be specified instead of "
                         << strfactor << std::endl;
      ABORT("getScalarOrFilterData: either a value or a valid variable should be specified");
    }
    data.get(var, factors);
  }
  return factors;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
