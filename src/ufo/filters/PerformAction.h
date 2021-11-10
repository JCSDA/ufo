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

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters controlling the Perform Action filter.
class PerformActionParameters : public FilterParametersBaseWithAbstractActions {
  OOPS_CONCRETE_PARAMETERS(PerformActionParameters, FilterParametersBaseWithAbstractActions)

 public:
  // Import both overloads of deserialize() from the base class. We will override one of them.
  using FilterParametersBaseWithAbstractActions::deserialize;

  /// \brief Load the values of all previously registered parameters from the (not necessarily
  /// top-level) configuration \p config.
  ///
  /// The base class implementation is overridden to throw an exception unless exactly one of the
  /// `action` and `actions` keys is present in the input configuration `config`.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  /// \brief Return parameters specifying the actions to be performed on observations flagged by the
  /// filter.
  ///
  /// This implementation returns the actions specified in the `action` or `actions` option,
  /// depending on which is present.
  std::vector<std::unique_ptr<FilterActionParametersBase>> actions() const override;
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
