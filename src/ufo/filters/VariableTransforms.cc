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
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/filters/VariableTransforms.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

VariableTransforms::VariableTransforms(
    ioda::ObsSpace& obsdb, const eckit::Configuration& config,
    std::shared_ptr<ioda::ObsDataVector<int>> flags,
    std::shared_ptr<ioda::ObsDataVector<float>> obserr)
    : FilterBase(obsdb, config, flags, obserr), parameters_()
{
  // Create parameters for this transformation
  parameters_ = TransformFactory::createParameters(config.getString("Transform"));
  parameters_->validateAndDeserialize(config);

  // Create the transform
  std::unique_ptr<TransformBase> transform =
    TransformFactory::create(parameters_->Transform.value(), *parameters_, data_, flags_, obserr_);

  // Add any required transform variables to allvars_
  allvars_ += Variables(filtervars_);
  allvars_ += transform->requiredVariables();

  oops::Log::debug() << this << std::endl;
}

// -----------------------------------------------------------------------------

VariableTransforms::~VariableTransforms() {}

// -----------------------------------------------------------------------------

void VariableTransforms::applyFilter(
    const std::vector<bool>& apply, const Variables&,
    std::vector<std::vector<bool>>&) const {
  print(oops::Log::trace());
  oops::Log::debug() << " --> In variabletransforms::applyFilter" << std::endl;

  // Do not perform transformation if there are no observations present.
  if (parameters_->SkipWhenNoObs.value() && data_.nlocs() == 0) {
    oops::Log::debug() << " --> No observations present. "
                       << "Transformation will not be performed." << std::endl;
    return;
  }

  // Create the transform again.  This is because data_ is updated by the filter
  // but not by the copy held by the variable transform.
  oops::Log::debug() << "     --> set Transform object" << std::endl;
  std::unique_ptr<TransformBase> transform =
    TransformFactory::create(parameters_->Transform.value(), *parameters_, data_, flags_, obserr_);

  // Run all calculations requested
  oops::Log::debug() << "         estimate: " << parameters_->Transform.value() << std::endl;
  transform->runTransform(apply);
}

// -----------------------------------------------------------------------------

void VariableTransforms::print(std::ostream& os) const {
  os << "VariableTransforms: config = " << *parameters_ << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
