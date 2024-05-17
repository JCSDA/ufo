/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/compositeoper/ObsCompositeTLAD.h"

#include <ostream>
#include <utility>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/compositeoper/ObsCompositeParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsCompositeTLAD> makerCompositeTL_("Composite");
// -----------------------------------------------------------------------------

ObsCompositeTLAD::ObsCompositeTLAD(const ioda::ObsSpace & odb, const Parameters_ & parameters)
  : LinearObsOperatorBase(odb)
{
  oops::Log::trace() << "ObsCompositeTLAD constructor starting" << std::endl;

  for (const eckit::LocalConfiguration &operatorConfig : parameters.components.value()) {
    LinearObsOperatorParametersWrapper operatorParams;
    operatorParams.validateAndDeserialize(operatorConfig);
    std::unique_ptr<LinearObsOperatorBase> op(
          LinearObsOperatorFactory::create(odb, operatorParams.operatorParameters));
    requiredVars_ += op->requiredVars();
    components_.push_back(std::move(op));
  }

  oops::Log::trace() << "ObsCompositeTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsCompositeTLAD::~ObsCompositeTLAD() {
  oops::Log::trace() << "ObsCompositeTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCompositeTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & ydiags,
                                     const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsCompositeTLAD: setTrajectory entered" << std::endl;

  for (const std::unique_ptr<LinearObsOperatorBase> &component : components_) {
    component->setTrajectory(geovals, ydiags, qc_flags);
  }
  oops::Log::trace() << "ObsCompositeTLAD: setTrajectory exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCompositeTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                     const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsCompositeTLAD: simulateObsTL entered" << std::endl;

  for (const std::unique_ptr<LinearObsOperatorBase> &component : components_) {
    component->simulateObsTL(geovals, ovec, qc_flags);
  }
  oops::Log::trace() << "ObsCompositeTLAD: simulateObsTL exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCompositeTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                     const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsCompositeTLAD: simulateObsAD entered" << std::endl;

  for (const std::unique_ptr<LinearObsOperatorBase> &component : components_) {
    component->simulateObsAD(geovals, ovec, qc_flags);
  }
  oops::Log::trace() << "ObsCompositeTLAD: simulateObsAD exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsVariables ObsCompositeTLAD::simulatedVars() const {
  // Merge the lists of variables simulated by all components, ensuring that there are
  // no overlaps.
  oops::ObsVariables vars;
  for (const std::unique_ptr<LinearObsOperatorBase> &component : components_) {
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

void ObsCompositeTLAD::print(std::ostream & os) const {
  os << "ObsComposite with the following components:\n";
  for (size_t i = 0; i < components_.size(); ++i) {
    os << "  - " << *components_[i];
    if (i != components_.size() - 1)
      os << '\n';
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
