/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICALDATA_H_
#define UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICALDATA_H_

#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "oops/util/Logger.h"

#include "ufo/operators/categoricaloper/ObsCategoricalParameters.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

// Define some traits that enable ObsCategoricalData to work with both linear and nonlinear
// operators.

class LinearObsOperatorBase;
class LinearObsOperatorFactory;
class LinearObsOperatorParametersWrapper;

class ObsOperatorBase;
class ObsOperatorFactory;
class ObsOperatorParametersWrapper;

template <typename OPBASE>
struct ObsOperatorTraits {};

template <>
struct ObsOperatorTraits<ObsOperatorBase> {
  typedef ObsOperatorFactory Factory_;
  typedef ObsOperatorParametersWrapper ParametersWrapper_;
};

template <>
struct ObsOperatorTraits<LinearObsOperatorBase> {
  typedef LinearObsOperatorFactory Factory_;
  typedef LinearObsOperatorParametersWrapper ParametersWrapper_;
};

/// \brief Data handler for the Categorical observation operator and TL/AD code.
template <typename OPBASE>
class ObsCategoricalData  {
 public:
  /// Get all information related to the configuration of the Categorical operator and TL/AD code.
  void configure(const ioda::ObsSpace & odb,
                 const ObsCategoricalParameters & parameters)
  {
    // Get categorical variable from ObsSpace (and throw an exception if it is not present).
    // In the ObsSpace, the categorical variable can be either a vector of strings or
    // a vector of integers; if the latter, it is converted to a vector of strings here.
    const std::string &categoricalVariableName = parameters.categoricalVariable.value();
    oops::Log::debug() << "categorical variable: " << categoricalVariableName << std::endl;
    categoricalVariable_.assign(odb.nlocs(), "");
    if (odb.has("MetaData", categoricalVariableName)) {
      const ioda::ObsDtype dtype = odb.dtype("MetaData", categoricalVariableName);
      if (dtype == ioda::ObsDtype::String) {
        odb.get_db("MetaData", categoricalVariableName, categoricalVariable_);
      } else if (dtype == ioda::ObsDtype::Integer) {
        std::vector <int> categoricalVariableInt(odb.nlocs());
        odb.get_db("MetaData", categoricalVariableName, categoricalVariableInt);
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

    // Fill vector of operator names.
    for (size_t jloc = 0; jloc < odb.nlocs(); ++jloc) {
      auto it_operName = categorisedOperatorNames_.find(categoricalVariable_[jloc]);
      const auto &operName = (it_operName != categorisedOperatorNames_.end() ?
                              it_operName->second :
                              fallbackOperatorName_);
      locOperNames_.emplace_back(operName);
    }

    // Create list of component operators.
    const std::vector<eckit::LocalConfiguration> & operatorConfigs =
      parameters.operatorConfigurations.value();

    // If a list of labels has been specified, ensure it is the correct length.
    if (parameters.operatorLabels.value() != boost::none &&
        operatorConfigs.size() != parameters.operatorLabels.value().value().size())
        throw eckit::UserError("Incorrect number of operator labels specified", Here());

    // Configure each component operator.
    for (std::size_t jop = 0; jop < operatorConfigs.size(); ++jop) {
      const eckit::LocalConfiguration & operatorConfig = operatorConfigs[jop];
      typedef typename ObsOperatorTraits<OPBASE>::Factory_ Factory_;
      typedef typename ObsOperatorTraits<OPBASE>::ParametersWrapper_ ParametersWrapper_;

      ParametersWrapper_ operatorParams;
      operatorParams.validateAndDeserialize(operatorConfig);
      std::unique_ptr<OPBASE> op(Factory_::create(odb, operatorParams.operatorParameters));
      requiredVars_ += op->requiredVars();

      const std::string operatorLabel = parameters.operatorLabels.value() != boost::none ?
        parameters.operatorLabels.value().value()[jop] : operatorConfig.getString("name");
      components_.emplace(std::make_pair(operatorLabel, std::move(op)));
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

    // Check that there are no duplicate component operators.
    if (components_.size() != operatorConfigs.size()) {
      throw eckit::UserError("There are at least two duplicate component operators. Consider using "
                             "the 'operator labels' configuration option to differentiate between "
                             "them", Here());
    }
  }

  /// Return required variables for the operator.
  const oops::Variables & requiredVars() const {return requiredVars_;}

  /// Return component operators.
  const std::map<std::string, std::unique_ptr<OPBASE>> & components() const {return components_;}

  /// Return list of operator names to use at each location.
  const std::vector<std::string> & locOperNames() const {return locOperNames_;}

  /// Fill final H(x) vector from a list of components.
  void fillHofX(const std::map <std::string, ioda::ObsVector> & ovecs,
                ioda::ObsVector & ovec) const {
    // Insert values into ovec according to the categorical variable.
    // Use the fallback operator when necessary.
    for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
      const auto &ovecloc = ovecs.at(locOperNames_[jloc]);
      // Loop over each variable at this location.
      for (size_t jvar = 0; jvar < ovec.nvars(); ++jvar) {
        const size_t idx = jloc * ovec.nvars() + jvar;
        ovec[idx] = ovecloc[idx];
      }
    }
  }

  void print(std::ostream & os) const {
    os << "ObsCategorical operator:" << std::endl;
    os << "- Fallback operator: " << fallbackOperatorName_ << std::endl;
    os << "- Categorised operators: " << categorisedOperatorNames_ << std::endl;
  }

 private:
  /// Required variables.
  oops::Variables requiredVars_;

  /// Observation operators which are run by the Categorical operator.
  std::map<std::string, std::unique_ptr<OPBASE>> components_;

  /// Value of the categorical variable in the ObsSpace.
  std::vector <std::string> categoricalVariable_;

  /// Name of the fallback observation operator.
  std::string fallbackOperatorName_;

  /// Names of the categorised observation operators.
  std::map<std::string, std::string> categorisedOperatorNames_;

  /// Operator name at each location.
  std::vector<std::string> locOperNames_;
};

}  // namespace ufo
#endif  // UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICALDATA_H_
