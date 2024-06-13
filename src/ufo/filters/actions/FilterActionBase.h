/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_FILTERACTIONBASE_H_
#define UFO_FILTERS_ACTIONS_FILTERACTIONBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ioda {
template <typename T> class ObsDataVector;
}  // namespace ioda

namespace ufo {

class FilterActionFactory;
class ObsFilterData;
class Variables;

// -----------------------------------------------------------------------------
/// Parameters controlling a filter action.
class FilterActionParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(FilterActionParametersBase, Parameters)

 public:
  /// \brief Name of the action to be performed.
  ///
  /// \note This parameter is marked as optional because it may not be required in certain
  /// circumstances, e.g. in tests instantiating a particular FilterAction directly (rather than
  /// via the FilterActionFactory). FilterActionFactory will throw an exception if this parameter is
  /// not provided.
  oops::OptionalParameter<std::string> name{"name", this};
};

// -----------------------------------------------------------------------------
/// Base class for actions performed on observations flagged by filters.
///
/// Note: each concrete implementation should typedef `Parameters_` to the name of a subclass of
/// FilterActionParametersBase encapsulating its configuration options. It should also provide
/// a constructor with the following signature:
///
///     FilterActionBase(const Parameters_ &);
class FilterActionBase : private boost::noncopyable {
 public:
  FilterActionBase() {}
  virtual ~FilterActionBase() {}

  /// \brief Perform the action.
  ///
  /// \param vars
  ///   The list of filter variables.
  /// \param flagged
  ///   If flagged[i][j] is true, it means that the action should be performed on jth observation
  ///   of ith filter variable.
  /// \param data
  ///   Accessor to obs filter data.
  /// \param filterQCflag
  ///   QC flag identifying observations rejected by the type of filter performing the action.
  ///   (Relevant only for actions rejecting observations.)
  /// \param flags
  ///   QC flags of all "simulated variables".
  /// \param obserr
  ///   Obs error estimates of all "simulated variables".
  virtual void apply(const ufo::Variables &vars, const std::vector<std::vector<bool>> &flagged,
                     ObsFilterData &data, int filterQCflag,
                     ioda::ObsDataVector<int> &flags, ioda::ObsDataVector<float> &obserr) const = 0;

  /// \brief Return the list of variables required by the action.
  ///
  /// This list must in particular contain any required variables that become available to
  /// filters only after FilterBase::priorFilter() or FilterBase::postFilter() is called; this
  /// includes `GeoVaLs`, `HofX` and `ObsDiag`.
  virtual const Variables & requiredVariables() const = 0;

  /// \brief Return true if this action modifies QC flags.
  ///
  /// When a filter executes multiple actions, only the last is allowed to modify QC flags.
  virtual bool modifiesQCFlags() const = 0;
};

// -----------------------------------------------------------------------------

/// Filter action factory.
class FilterActionFactory {
 public:
  virtual ~FilterActionFactory() = default;

  /// \brief Create and return a new filter action.
  ///
  /// The action's type is determined by the \c name attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of FilterActionParametersBase
  /// associated with that action type, otherwise an exception will be thrown.
  static std::unique_ptr<FilterActionBase> create(const FilterActionParametersBase &parameters);

  /// \brief Create and return an instance of the subclass of FilterActionParametersBase
  /// storing parameters of actions of the specified type.
  static std::unique_ptr<FilterActionParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all actions that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

 protected:
  /// \brief Register a maker able to create actions of type \p name.
  explicit FilterActionFactory(const std::string &name);

 private:
  virtual std::unique_ptr<FilterActionBase> make(const FilterActionParametersBase &) = 0;

  virtual std::unique_ptr<FilterActionParametersBase> makeParameters() const = 0;

  static std::map < std::string, FilterActionFactory * > & getMakers() {
    static std::map < std::string, FilterActionFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class FilterActionMaker : public FilterActionFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<FilterActionBase> make(const FilterActionParametersBase & parameters) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return boost::make_unique<T>(stronglyTypedParameters);
  }

  std::unique_ptr<FilterActionParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit FilterActionMaker(const std::string & name) : FilterActionFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_FILTERACTIONBASE_H_
