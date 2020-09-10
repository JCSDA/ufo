/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIASCOVARIANCE_H_
#define UFO_OBSBIASCOVARIANCE_H_

#include <map>
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

// Constructor, destructor
  ObsBiasCovariance(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsBiasCovariance() {}

// Linear algebra operators
  void linearize(const ObsBias &, const eckit::Configuration &);
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void randomize(ObsBiasIncrement &) const;

// Utilities
  const eckit::Configuration & config() const {return conf_;}
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &);

 private:
  void print(std::ostream &) const {}
  const eckit::LocalConfiguration conf_;
  const ioda::ObsSpace & odb_;

// QCed obs numbers <channel>
  std::vector<std::size_t> obs_num_prior_;

// Minimal required QCed obs number to add contribution
  std::size_t minimal_required_obs_number_;

// Analysis error variances from prev cycle
  std::vector<double> variances_prior_;

// Active variances
  std::vector<double> variances_;

// Default smallest variance value
  double smallest_variance_ = 1.0e-6;

// Default smallest variance value
  double largest_variance_ = 10.0;

  std::vector<std::string> prednames_;
  std::vector<int> jobs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASCOVARIANCE_H_
