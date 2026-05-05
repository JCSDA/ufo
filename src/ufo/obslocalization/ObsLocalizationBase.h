/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace ufo {

/// Base class for observation-space localization.
template<typename ITERATOR>
class ObsLocalizationBase : public util::Printable,
                            private boost::noncopyable {
 public:
  ObsLocalizationBase() = default;
  virtual ~ObsLocalizationBase() = default;

  /// compute obs-space localization: update \p locfactor with observation-space
  /// localization values between observations and \p point in model-space.
  /// Set \p locfactor to missing value for observations that are not local.
  virtual void computeLocalization(const ITERATOR & point, ioda::ObsVector & locfactor) const = 0;

  /// compute localization between an observation and a point in model-space, or
  /// between two observations. Return values should be between 0.0 and 1.0, inclusive,
  /// with 0 indicating that the points are outside of localization.
  virtual double computeLocalization(const eckit::geometry::Point3 &,
                                     const eckit::geometry::Point3 &) const {
    // NOTE, this really should be a pure virtual method, but since not all
    // ObsLocalizationBase derived classes will be implementing this method yet,
    // we are providing a default implementation.
    ABORT("computeLocalization(Point3,Point3) not implemented for this ObsLocalization class");
    return 0.0;}
};

// =============================================================================

/// ObsLocalization Factory
template <typename ITERATOR>
class ObsLocalizationFactory {
 public:
  static std::unique_ptr<ObsLocalizationBase<ITERATOR>> create(const eckit::Configuration &,
                                                               ioda::ObsSpace &);
 protected:
  explicit ObsLocalizationFactory(const std::string &);
  virtual ~ObsLocalizationFactory() = default;
 private:
  virtual ObsLocalizationBase<ITERATOR> * make(const eckit::Configuration &,
                                               ioda::ObsSpace &) = 0;
  static std::map < std::string, ObsLocalizationFactory<ITERATOR> * > & getMakers() {
    static std::map < std::string, ObsLocalizationFactory<ITERATOR> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class ITERATOR, class T>
class ObsLocalizationMaker : public ObsLocalizationFactory<ITERATOR> {
  virtual ObsLocalizationBase<ITERATOR> * make(const eckit::Configuration & conf,
                                               ioda::ObsSpace & obspace)
    { return new T(conf, obspace); }
 public:
  explicit ObsLocalizationMaker(const std::string & name) :
    ObsLocalizationFactory<ITERATOR>(name) {}
};

// -----------------------------------------------------------------------------

template <typename ITERATOR>
ObsLocalizationFactory<ITERATOR>::ObsLocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs localization factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename ITERATOR>
std::unique_ptr<ObsLocalizationBase<ITERATOR>> ObsLocalizationFactory<ITERATOR>::create(
                              const eckit::Configuration & conf, ioda::ObsSpace & obspace) {
  oops::Log::trace() << "ObsLocalizationBase<ITERATOR>::create starting" << std::endl;
  const std::string id = conf.getString("localization method");
  typename std::map<std::string, ObsLocalizationFactory<ITERATOR>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in obs localization factory.");
  }
  std::unique_ptr<ObsLocalizationBase<ITERATOR>> ptr(jloc->second->make(conf, obspace));
  oops::Log::trace() << "ObsLocalizationBase<ITERATOR>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
