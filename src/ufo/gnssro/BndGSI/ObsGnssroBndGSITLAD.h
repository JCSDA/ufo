/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GNSSRO_BNDGSI_OBSGNSSROBNDGSITLAD_H_
#define UFO_GNSSRO_BNDGSI_OBSGNSSROBNDGSITLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/gnssro/BndGSI/ObsGnssroBndGSITLAD.interface.h"
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
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// GnssroBndGSI observation operator
class ObsGnssroBndGSITLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGnssroBndGSITLAD> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroBndGSITLAD";}

  ObsGnssroBndGSITLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroBndGSITLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperGnssroBndGSI_;}
  const int & toFortran() const {return keyOperGnssroBndGSI_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperGnssroBndGSI_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_GNSSRO_BNDGSI_OBSGNSSROBNDGSITLAD_H_
