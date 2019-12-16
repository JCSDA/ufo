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

#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/obsbias/ObsBiasBase.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVals;
  class ObsBiasIncrement;
  class ObsDiagnostics;

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

class ObsBias : public util::Printable,
                private util::ObjectCounter<ObsBias> {
 public:
  static const std::string classname() {return "ufo::ObsBias";}

  ObsBias(ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsBias() {}

  ObsBias & operator+=(const ObsBiasIncrement &);
  ObsBias & operator=(const ObsBias &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::size_t size() const;

/// Bias parameters interface
  const double & operator[](const unsigned int ii) const {return (*biasbase_)[ii];}

/// Obs bias model
  void computeObsBias(const GeoVaLs &, ioda::ObsVector &, const ObsDiagnostics &) const;

/// Obs Bias Predictors
  void computeObsBiasPredictors(const GeoVaLs &, const ObsDiagnostics &,
                                std::unique_ptr<ioda::ObsDataVector<float>> &) const;

/// Other
  const oops::Variables & requiredGeoVaLs() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}
  const eckit::Configuration & config() const {return conf_;}
  const ioda::ObsSpace & obspace() const {return biasbase_->obspace();}

/// Operator
  operator bool() const {return biasbase_.get();}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsBiasBase> biasbase_;
  const eckit::LocalConfiguration conf_;
  oops::Variables geovars_;
  oops::Variables hdiags_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_H_
