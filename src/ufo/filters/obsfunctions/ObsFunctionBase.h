/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONBASE_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

class Variables;

// -----------------------------------------------------------------------------
/// Base class for computing functions on observation data

class ObsFunctionBase : private boost::noncopyable {
 public:
  explicit ObsFunctionBase(const eckit::LocalConfiguration conf = eckit::LocalConfiguration()) {}
  virtual ~ObsFunctionBase() {}

/// compute the result of the function
  virtual void compute(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const = 0;

/// geovals required to compute the function
  virtual const ufo::Variables & requiredVariables() const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Function Factory
class ObsFunctionFactory {
 public:
  static ObsFunctionBase * create(const Variable &);
  virtual ~ObsFunctionFactory() { getMakers().clear(); }
  static bool functionExists(const std::string &);
 protected:
  explicit ObsFunctionFactory(const std::string &);
 private:
  virtual ObsFunctionBase * make(const eckit::LocalConfiguration conf) = 0;
  static std::map < std::string, ObsFunctionFactory * > & getMakers() {
    static std::map < std::string, ObsFunctionFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsFunctionMaker : public ObsFunctionFactory {
  virtual ObsFunctionBase * make(const eckit::LocalConfiguration conf)
    { return new T(conf); }
 public:
  explicit ObsFunctionMaker(const std::string & name)
    : ObsFunctionFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONBASE_H_
