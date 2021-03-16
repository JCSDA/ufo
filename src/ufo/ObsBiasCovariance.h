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

#include "oops/base/Variables.h"
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
  ObsBiasCovariance(ioda::ObsSpace & odb,
                    const eckit::Configuration & biasConf);
  ~ObsBiasCovariance() {}

// Linear algebra operators
  void linearize(const ObsBias &, const eckit::Configuration &);
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void randomize(ObsBiasIncrement &) const;

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &);
  const std::vector<std::string> predictorNames() const {return prednames_;}

 private:
  void print(std::ostream &) const {}

  ioda::ObsSpace & odb_;

// Hessian contribution from Jo bias correction terms
  std::vector<double> ht_rinv_h_;

// preconditioner
  std::vector<double> preconditioner_;

// QCed obs numbers <channel>
  std::vector<std::size_t> obs_num_;

// Minimal required QCed obs number to add contribution
  std::size_t minimal_required_obs_number_;

// Analysis error variances
  std::vector<double> analysis_variances_;

// Error variances
  std::vector<double> variances_;

// Default smallest variance value
  double smallest_variance_ = 1.0e-6;

// Default largest variance value
  double largest_variance_ = 10.0;

// Default largest analysis error variance
  double largest_analysis_variance_ = 10000.0;

// Default stepsize
  double step_size_ = 1.e-4;

  std::vector<std::string> prednames_;

  /// variables for which bias correction coefficients will be updated
  oops::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIASCOVARIANCE_H_
