/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SATTCWV_SATTCWVTLAD_H_
#define UFO_SATTCWV_SATTCWVTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/sattcwv/SatTCWVTLAD.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Precipitable water observation operator
class SatTCWVTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<SatTCWVTLAD> {
 public:
  static const std::string classname() {return "ufo::SatTCWVTLAD";}

  SatTCWVTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~SatTCWVTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperSatTCWV_;}
  const int & toFortran() const {return keyOperSatTCWV_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperSatTCWV_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_SATTCWV_SATTCWVTLAD_H_
