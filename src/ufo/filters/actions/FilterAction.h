/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_FILTERACTION_H_
#define UFO_FILTERS_ACTIONS_FILTERACTION_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "ufo/filters/actions/FilterActionBase.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

class ObsFilterData;

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
             const ObsFilterData &data, int filterQCflag,
             ioda::ObsDataVector<int> &flags, ioda::ObsDataVector<float> &obserr) const;
  virtual const ufo::Variables & requiredVariables() const;

 private:
  std::unique_ptr<FilterActionBase> action_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_FILTERACTION_H_
