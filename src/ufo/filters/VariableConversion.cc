/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/calculate/CalculateBase.h"
#include "ufo/filters/VariableConversion.h"
#include "ufo/filters/VariableConversionParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------

VariableConversion::VariableConversion(
    ioda::ObsSpace& obsdb, const eckit::Configuration& config,
    std::shared_ptr<ioda::ObsDataVector<int>> flags,
    std::shared_ptr<ioda::ObsDataVector<float>> obserr)
    : FilterBase(obsdb, config, flags, obserr)
{
  options_.reset(new VariableConversionParameters());
  options_->deserialize(config);
  allvars_ += Variables(filtervars_);
}

// -----------------------------------------------------------------------------

VariableConversion::~VariableConversion() {}

// -----------------------------------------------------------------------------

void VariableConversion::applyFilter(
    const std::vector<bool>& apply, const Variables& filtervars,
    std::vector<std::vector<bool>>& flagged) const {
  print(oops::Log::trace());
  std::cout << " --> In VariableConversion::applyFilter" << std::endl;
  std::cout << "     --> set Calculate object" << std::endl;

  // Run all calculations requested
  for (const auto& cal_ : options_->Calculate.value()) {
    std::cout << "         estimate: " << cal_ << std::endl;
    std::unique_ptr<CalculateBase> calculate =
        CalculateFactory::create(cal_, *options_, obsdb_, flags_);

    calculate->runCalculate();
  }
}

// -----------------------------------------------------------------------------

void VariableConversion::print(std::ostream& os) const {
  os << "VariableConversion: config = " << config_ << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
