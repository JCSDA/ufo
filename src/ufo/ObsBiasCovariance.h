/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIASCOVARIANCE_H_
#define UFO_OBSBIASCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------

class ObsBiasCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasCovariance> {
 public:
  static const std::string classname() {return "ufo::ObsBiasCovariance";}

/// Constructor, destructor
  explicit ObsBiasCovariance(const eckit::Configuration &) {}
  ~ObsBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ObsBias &) {}
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const {}
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const {}
  void randomize(ObsBiasIncrement &) const {}

 private:
  void print(std::ostream &) const {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASCOVARIANCE_H_
