/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PERFORMACTION_H_
#define UFO_FILTERS_PERFORMACTION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters controlling the action performed on observations flagged by a filter.
///
/// The action must be specified explicitly (there is no default action).
class RequiredFilterActionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(RequiredFilterActionParameters, Parameters)

 public:
  /// After deserialization, holds an instance of a subclass of FilterActionParametersBase
  /// controlling the behavior of a filter action. The type of the subclass is determined
  /// by the value of the `name` key in the Configuration object from which this object
  /// is deserialized.
  oops::RequiredPolymorphicParameter<FilterActionParametersBase, FilterActionFactory>
    actionParameters{"name", this};
};

/// \brief Parameters controlling the Perform Action filter.
class PerformActionParameters : public FilterParametersBaseWithAbstractAction {
  OOPS_CONCRETE_PARAMETERS(PerformActionParameters, FilterParametersBaseWithAbstractAction)

 public:
  /// Parameters controlling the action performed on observations flagged by the filter.
  const FilterActionParametersBase &action() const override {
    return action_.value().actionParameters.value();
  }

 private:
  /// Parameters controlling the action performed on observations flagged by the filter.
  oops::RequiredParameter<RequiredFilterActionParameters> action_{"action", this};
};

/// \brief Perform the action specified in the `action` YAML section on each observation selected by
/// the `where` statement.
///
/// In contrast to other filters, this filter requires the action to be specified explicitly in the
/// YAML file; there is no default action.
class PerformAction : public FilterBase,
                      private util::ObjectCounter<PerformAction> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef PerformActionParameters Parameters_;

  static const std::string classname() {return "ufo::PerformAction";}

  PerformAction(ioda::ObsSpace &, const Parameters_ &,
                std::shared_ptr<ioda::ObsDataVector<int> >,
                std::shared_ptr<ioda::ObsDataVector<float> >);

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::black;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PERFORMACTION_H_
