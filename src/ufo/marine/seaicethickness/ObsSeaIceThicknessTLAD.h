/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESSTLAD_H_
#define UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESSTLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/marine/seaicethickness/ObsSeaIceThicknessTLAD.interface.h"

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
/// Sea-ice thickness TL/AD observation operator class
class ObsSeaIceThicknessTLAD : public LinearObsOperatorBase,
                               private util::ObjectCounter<ObsSeaIceThicknessTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsSeaIceThicknessTLAD";}

  ObsSeaIceThicknessTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaIceThicknessTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESSTLAD_H_
