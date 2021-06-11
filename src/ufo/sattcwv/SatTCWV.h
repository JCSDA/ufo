/*
 *
 * (C) Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_SATTCWV_SATTCWV_H_
#define UFO_SATTCWV_SATTCWV_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/sattcwv/SatTCWV.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// SatTCWV observation operator
class SatTCWV : public ObsOperatorBase,
                private util::ObjectCounter<SatTCWV> {
 public:
  static const std::string classname() {return "ufo::SatTCWV";}

  SatTCWV(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~SatTCWV();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperSatTCWV_;}
  const int & toFortran() const {return keyOperSatTCWV_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperSatTCWV_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_SATTCWV_SATTCWV_H_
