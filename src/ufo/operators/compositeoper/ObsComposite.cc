/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/compositeoper/ObsComposite.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/compositeoper/ObsCompositeParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsComposite> obsCompositeMaker_("Composite");
// -----------------------------------------------------------------------------

ObsComposite::ObsComposite(const ioda::ObsSpace & odb, const Parameters_ & parameters)
  : ObsOperatorBase(odb), odb_(odb)
{
  oops::Log::trace() << "ObsComposite constructor starting" << std::endl;

  for (const eckit::LocalConfiguration &operatorConfig : parameters.components.value()) {
    ObsOperatorParametersWrapper operatorParams;
    operatorParams.validateAndDeserialize(operatorConfig);
    std::unique_ptr<ObsOperatorBase> op(
          ObsOperatorFactory::create(odb, operatorParams.operatorParameters));
    requiredVars_ += op->requiredVars();
    components_.push_back(std::move(op));
  }

  oops::Log::trace() << "ObsComposite constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsComposite::~ObsComposite() {
  oops::Log::trace() << "ObsComposite destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsComposite::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                               ObsDiagnostics & ydiags, const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsComposite: simulateObs entered" << std::endl;

  for (const std::unique_ptr<ObsOperatorBase> &component : components_)
  {
    component->simulateObs(gv, ovec, ydiags, qc_flags);
  }
  oops::Log::trace() << "ObsComposite: simulateObs exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsVariables ObsComposite::simulatedVars() const {
  // Merge the lists of variables simulated by all components, ensuring that there are
  // no overlaps.
  oops::ObsVariables vars;
  for (const std::unique_ptr<ObsOperatorBase> &component : components_) {
    oops::ObsVariables componentVars;
    // We use += rather than = to make sure componentVars contains no duplicate entries.
    componentVars += component->simulatedVars();
    const size_t oldSize = vars.size();
    vars += componentVars;
    if (vars.size() != oldSize + componentVars.size())
      // We don't want multiple components to write to the same row in the H(x) array.
      throw eckit::UserError("Multiple components simulate the same variables", Here());
  }
  return vars;
}

// -----------------------------------------------------------------------------

void ObsComposite::print(std::ostream & os) const {
  os << "ObsComposite with the following components:\n";
  for (size_t i = 0; i < components_.size(); ++i) {
    os << "  - " << *components_[i];
    if (i != components_.size() - 1)
      os << '\n';
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
