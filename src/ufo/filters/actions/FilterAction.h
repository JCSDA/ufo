/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_FILTERACTION_H_
#define UFO_FILTERS_ACTIONS_FILTERACTION_H_

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

namespace ioda {
  template<typename DATATYPE> class ObsDataVector;
}

namespace ufo {
  class FilterActionBase;
  class FilterActionParametersBase;
  class ObsFilterData;
  class Variables;

// -----------------------------------------------------------------------------

class FilterAction : private boost::noncopyable {
 public:
  explicit FilterAction(const FilterActionParametersBase &);
  ~FilterAction();

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
  void apply(const ufo::Variables &vars, const std::vector<std::vector<bool>> &flagged,
             ObsFilterData &data, int filterQCflag,
             ioda::ObsDataVector<int> &flags, ioda::ObsDataVector<float> &obserr) const;
  const ufo::Variables & requiredVariables() const;

  /// \brief Return true if this action modifies QC flags.
  ///
  /// When a filter executes multiple actions, only the last is allowed to modify QC flags.
  bool modifiesQCFlags() const;

 private:
  std::unique_ptr<FilterActionBase> action_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_FILTERACTION_H_
