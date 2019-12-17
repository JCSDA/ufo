/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIASCOVARIANCE_H_
#define UFO_OBSBIASCOVARIANCE_H_

#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
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
  ObsBiasCovariance(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ObsBias &);
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void randomize(ObsBiasIncrement &) const;

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream &) const {}
  const eckit::LocalConfiguration conf_;
  std::vector<double> variance_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASCOVARIANCE_H_
