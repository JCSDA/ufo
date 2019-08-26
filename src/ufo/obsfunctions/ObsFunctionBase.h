/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSFUNCTIONS_OBSFUNCTIONBASE_H_
#define UFO_OBSFUNCTIONS_OBSFUNCTIONBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Base class for computing functions on observation data

class ObsFunctionBase : private boost::noncopyable {
 public:
  ObsFunctionBase() {}
  virtual ~ObsFunctionBase() {}

/// compute the result of the function
  virtual void compute(const ioda::ObsDataVector<float> &,
                       const ioda::ObsDataVector<float> &,
                       ioda::ObsDataVector<float> &) const = 0;

/// observed variables required to compute the function
  virtual const oops::Variables & requiredObsData() const = 0;
/// metadata requited to compute the function
  virtual const oops::Variables & requiredMetaData() const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Function Factory
class ObsFunctionFactory {
 public:
  static ObsFunctionBase * create(const std::string &);
  virtual ~ObsFunctionFactory() { getMakers().clear(); }
 protected:
  explicit ObsFunctionFactory(const std::string &);
 private:
  virtual ObsFunctionBase * make() = 0;
  static std::map < std::string, ObsFunctionFactory * > & getMakers() {
    static std::map < std::string, ObsFunctionFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsFunctionMaker : public ObsFunctionFactory {
  virtual ObsFunctionBase * make()
    { return new T(); }
 public:
  explicit ObsFunctionMaker(const std::string & name) : ObsFunctionFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSFUNCTIONS_OBSFUNCTIONBASE_H_
