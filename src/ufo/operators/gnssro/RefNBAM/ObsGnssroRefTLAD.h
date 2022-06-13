/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_REFNBAM_OBSGNSSROREFTLAD_H_
#define UFO_OPERATORS_GNSSRO_REFNBAM_OBSGNSSROREFTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/gnssro/RefNBAM/ObsGnssroRefTLAD.interface.h"

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
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// GnssroRef observation operator
class ObsGnssroRefTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGnssroRefTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroRefTLAD";}

  ObsGnssroRefTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroRefTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroRef_;}
  const int & toFortran() const {return keyOperGnssroRef_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroRef_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_REFNBAM_OBSGNSSROREFTLAD_H_
