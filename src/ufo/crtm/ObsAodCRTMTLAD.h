/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CRTM_OBSAODCRTMTLAD_H_
#define UFO_CRTM_OBSAODCRTMTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/crtm/ObsAodCRTMTLAD.interface.h"
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
/// AodCRTM (currently only temperature) observation for UFO.
class ObsAodCRTMTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsAodCRTMTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsAodCRTMTLAD";}

  ObsAodCRTMTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAodCRTMTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAodCRTM_;}
  const int & toFortran() const {return keyOperAodCRTM_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAodCRTM_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_CRTM_OBSAODCRTMTLAD_H_
