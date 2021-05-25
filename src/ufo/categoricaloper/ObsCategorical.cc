/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/categoricaloper/ObsCategorical.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/categoricaloper/ObsCategoricalParameters.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsCategorical> obsCategoricalMaker_("Categorical");
// -----------------------------------------------------------------------------

ObsCategorical::ObsCategorical(const ioda::ObsSpace & odb,
                               const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), odb_(odb)
{
  oops::Log::trace() << "ObsCategorical constructor starting" << std::endl;

  ObsCategoricalParameters parameters;
  parameters.validateAndDeserialize(config);

  // Get categorical variable from ObsSpace (and throw an exception if it is not present).
  // In the ObsSpace, the categorical variable can be either a vector of strings or
  // a vector of integers; if the latter, it is converted to a vector of strings here.
  const std::string &categoricalVariableName = parameters.categoricalVariable.value();
  oops::Log::debug() << "categorical variable: " << categoricalVariableName << std::endl;
  categoricalVariable_.assign(odb_.nlocs(), "");
  if (odb_.has("MetaData", categoricalVariableName)) {
    const ioda::ObsDtype dtype = odb_.dtype("MetaData", categoricalVariableName);
    if (dtype == ioda::ObsDtype::String) {
      odb_.get_db("MetaData", categoricalVariableName, categoricalVariable_);
    } else if (dtype == ioda::ObsDtype::Integer) {
      std::vector <int> categoricalVariableInt(odb_.nlocs());
      odb_.get_db("MetaData", categoricalVariableName, categoricalVariableInt);
      std::transform(categoricalVariableInt.cbegin(), categoricalVariableInt.cend(),
                     categoricalVariable_.begin(), [](int i) {return std::to_string(i);});
    } else {
      throw eckit::UserError("The categorical variable must be a vector of "
                             "either strings or integers", Here());
    }
  } else {
    throw eckit::UserError("The categorical variable " + categoricalVariableName +
                           " does not exist", Here());
  }

  // Name of fallback operator.
  fallbackOperatorName_ = parameters.fallbackOperatorName.value();
  oops::Log::debug() << "Fallback operator: " << fallbackOperatorName_ << std::endl;

  // Map of categorised operator names.
  categorisedOperatorNames_ = parameters.categorisedOperatorNames.value();
  oops::Log::debug() << "Categorised operators: " << categorisedOperatorNames_ << std::endl;

  // Create list of component operators.
  for (const eckit::LocalConfiguration &operatorConfig :
         parameters.operatorConfigurations.value()) {
    std::unique_ptr<ObsOperatorBase> op(ObsOperatorFactory::create(odb, operatorConfig));
    requiredVars_ += op->requiredVars();
    components_.emplace(std::make_pair(operatorConfig.getString("name"), std::move(op)));
  }

  // Check the fallback operator has been configured.
  if (components_.find(fallbackOperatorName_) == components_.end())
    throw eckit::UserError("The operator " + fallbackOperatorName_ +
                           " has not been configured", Here());

  // Check the categorised operators have been configured.
  for (const auto &operName : categorisedOperatorNames_)
    if (components_.find(operName.second) == components_.end())
      throw eckit::UserError("The operator " + operName.second +
                             " has not been configured", Here());

  oops::Log::trace() << "ObsCategorical constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsCategorical::~ObsCategorical() {
  oops::Log::trace() << "ObsCategorical destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategorical::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                 ObsDiagnostics & ydiags) const {
  oops::Log::trace() << "ObsCategorical: simulateObs entered" << std::endl;

  // Container of ObsVectors produced by each operator.
  std::map <std::string, ioda::ObsVector> ovecs;
  oops::Log::debug() << "Running operators" << std::endl;
  // Run each operator and store output in ovecs.
  for (const auto& component : components_) {
    ioda::ObsVector ovecTemp(ovec);
    component.second->simulateObs(gv, ovecTemp, ydiags);
    ovecs.insert({component.first, ovecTemp});
  }

  // Insert values into ovec according to the categorical variable.
  // Use the fallback operator when necessary.
  oops::Log::debug() << "Producing final ObsVector" << std::endl;
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    auto it_operName = categorisedOperatorNames_.find(categoricalVariable_[jloc]);
    const auto &operName = (it_operName != categorisedOperatorNames_.end() ?
                            it_operName->second :
                            fallbackOperatorName_);
    oops::Log::debug() << "Location " << jloc << ": operator name = " << operName << std::endl;
    const auto &ovecloc = ovecs.at(operName);
    // Loop over each variable at this location.
    for (size_t jvar = 0; jvar < ovec.nvars(); ++jvar) {
      const size_t idx = jloc * ovec.nvars() + jvar;
      ovec[idx] = ovecloc[idx];
    }
  }

  oops::Log::trace() << "ObsCategorical: simulateObs exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategorical::print(std::ostream & os) const {
  os << "ObsCategorical operator:" << std::endl;
  os << "- Fallback operator: " << fallbackOperatorName_ << std::endl;
  os << "- Categorised operators: " << categorisedOperatorNames_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
