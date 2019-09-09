/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ATMSFCINTERP_OBSATMSFCINTERPTLAD_H_
#define UFO_ATMSFCINTERP_OBSATMSFCINTERPTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmsfcinterp/ObsAtmSfcInterpTLAD.interface.h"
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

// -----------------------------------------------------------------------------
/// AtmSfcInterp TL/AD observation operator class
class ObsAtmSfcInterpTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAtmSfcInterpTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsAtmSfcInterpTLAD";}

  ObsAtmSfcInterpTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAtmSfcInterpTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const;

  // Other
  const oops::Variables & variables() const {return varin_;}

  int & toFortran() {return keyOperAtmSfcInterp_;}
  const int & toFortran() const {return keyOperAtmSfcInterp_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAtmSfcInterp_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_ATMSFCINTERP_OBSATMSFCINTERPTLAD_H_
