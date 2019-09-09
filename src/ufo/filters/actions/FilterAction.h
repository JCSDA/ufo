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

#include "ioda/ObsDataVector.h"
#include "ufo/filters/actions/FilterActionBase.h"
namespace eckit {
  class Configuration;
}

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class FilterAction : private boost::noncopyable {
 public:
  explicit FilterAction(const eckit::Configuration &);
  ~FilterAction();

  void apply(const oops::Variables &, const std::vector<std::vector<bool>> &,
             const ObsFilterData &,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const;
 private:
  std::unique_ptr<FilterActionBase> action_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_FILTERACTION_H_
