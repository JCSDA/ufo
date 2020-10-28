/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_LINEAROBSOPERATORBASE_H_
#define UFO_LINEAROBSOPERATORBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace ioda {
class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsBias;
class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Base class for observation operators

class LinearObsOperatorBase : public util::Printable,
                              private boost::noncopyable {
 public:
  LinearObsOperatorBase() {}
  virtual ~LinearObsOperatorBase() {}

/// Obs Operator
  virtual void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) = 0;
  virtual void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const = 0;
  virtual void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const = 0;

/// Operator input required from Model
  virtual const oops::Variables & requiredVars() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class LinearObsOperatorFactory {
 public:
  static LinearObsOperatorBase * create(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~LinearObsOperatorFactory() = default;
 protected:
  explicit LinearObsOperatorFactory(const std::string &);
 private:
  virtual LinearObsOperatorBase * make(const ioda::ObsSpace &, const eckit::Configuration &) = 0;
  static std::map < std::string, LinearObsOperatorFactory * > & getMakers() {
    static std::map < std::string, LinearObsOperatorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LinearObsOperatorMaker : public LinearObsOperatorFactory {
  virtual LinearObsOperatorBase * make(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & conf) {
    return new T(odb, conf);
  }
 public:
  explicit LinearObsOperatorMaker(const std::string & name) : LinearObsOperatorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_LINEAROBSOPERATORBASE_H_
