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

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/VariableTransforms.h"
#include "ufo/filters/VariableTransformsParameters.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

VariableTransforms::VariableTransforms(
    ioda::ObsSpace& obsdb, const eckit::Configuration& config,
    std::shared_ptr<ioda::ObsDataVector<int>> flags,
    std::shared_ptr<ioda::ObsDataVector<float>> obserr)
    : FilterBase(obsdb, config, flags, obserr)
{
  options_.reset(new VariableTransformsParameters());
  options_->deserialize(config);
  allvars_ += Variables(filtervars_);

  oops::Log::debug() << "VariableTransforms: config = " << config << std::endl;
}

// -----------------------------------------------------------------------------

VariableTransforms::~VariableTransforms() {}

// -----------------------------------------------------------------------------

void VariableTransforms::applyFilter(
    const std::vector<bool>& apply, const Variables& filtervars,
    std::vector<std::vector<bool>>& flagged) const {
  print(oops::Log::trace());
  std::cout << " --> In variabletransforms::applyFilter" << std::endl;
  std::cout << "     --> set Transform object" << std::endl;

  // Run all calculations requested
  for (const auto& cal_ : options_->Transform.value()) {
    std::cout << "         estimate: " << cal_ << std::endl;
    std::unique_ptr<TransformBase> Transform =
        TransformFactory::create(cal_, *options_, obsdb_, flags_, apply);

    Transform->runTransform();
  }
}

// -----------------------------------------------------------------------------

void VariableTransforms::print(std::ostream& os) const {
  os << "VariableTransforms: config = " << config_ << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
