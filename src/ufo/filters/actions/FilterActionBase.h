/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_FILTERACTIONBASE_H_
#define UFO_FILTERS_ACTIONS_FILTERACTIONBASE_H_

#include <map>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------
/// Base class for computing obs diagnostics

class FilterActionBase : private boost::noncopyable {
 public:
  FilterActionBase() {}
  virtual ~FilterActionBase() {}

/// compute the diagnostic
  virtual void apply(const oops::Variables &, const std::vector<std::vector<bool>> &,
                     const ObsFilterData &,
                     ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const = 0;
  virtual const Variables & requiredVariables() const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class FilterActionFactory {
 public:
  static FilterActionBase * create(const eckit::Configuration &);
  virtual ~FilterActionFactory() { getMakers().clear(); }
 protected:
  explicit FilterActionFactory(const std::string &);
 private:
  virtual FilterActionBase * make(const eckit::Configuration &) = 0;
  static std::map < std::string, FilterActionFactory * > & getMakers() {
    static std::map < std::string, FilterActionFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class FilterActionMaker : public FilterActionFactory {
  virtual FilterActionBase * make(const eckit::Configuration & conf)
    { return new T(conf); }
 public:
  explicit FilterActionMaker(const std::string & name) : FilterActionFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_FILTERACTIONBASE_H_
