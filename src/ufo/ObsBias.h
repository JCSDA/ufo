/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_H_
#define UFO_OBSBIAS_H_

#include <memory>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/obsbias/ObsBiasBase.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVals;
  class ObsBiasIncrement;

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

class ObsBias : public util::Printable,
                private boost::noncopyable,
                private util::ObjectCounter<ObsBias> {
 public:
  static const std::string classname() {return "ufo::ObsBias";}

  explicit ObsBias(const eckit::Configuration &);
  ObsBias(const ObsBias &, const bool);
  ~ObsBias() {}

  ObsBias & operator+=(const ObsBiasIncrement &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::size_t size() const;

  const double & operator[](const unsigned int ii) const {return (*biasbase_)[ii];}

/// Obs bias model
  void computeObsBias(const GeoVaLs &,
                      ioda::ObsVector &,
                      const ioda::ObsSpace &) const;

/// Other
  const oops::Variables & variables() const;
  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsBiasBase> biasbase_;
  const eckit::LocalConfiguration conf_;
  oops::Variables vars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
