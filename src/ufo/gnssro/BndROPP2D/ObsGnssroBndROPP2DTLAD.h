/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2DTLAD_H_
#define UFO_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2DTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/gnssro/BndROPP2D/ObsGnssroBndROPP2DTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

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
/// GnssroBndROPP2D observation operator
class ObsGnssroBndROPP2DTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGnssroBndROPP2DTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroBndROPP2DTLAD";}

  ObsGnssroBndROPP2DTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroBndROPP2DTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBndROPP2D_;}
  const int & toFortran() const {return keyOperGnssroBndROPP2D_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBndROPP2D_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_GNSSRO_BNDROPP2D_OBSGNSSROBNDROPP2DTLAD_H_
