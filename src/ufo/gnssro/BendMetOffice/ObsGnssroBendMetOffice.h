/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_H_
#define UFO_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/gnssro/BendMetOffice/ObsGnssroBendMetOffice.interface.h"
#include "ufo/ObsOperatorBase.h"
#include "ObsGnssroBendMetOfficeParameters.h"

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
/// GnssroBendMetOffice observation operator
// -----------------------------------------------------------------------------
class ObsGnssroBendMetOffice : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroBendMetOffice> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroBendMetOffice";}

  ObsGnssroBendMetOffice(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroBendMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBendMetOffice_;}
  const int & toFortran() const {return keyOperGnssroBendMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBendMetOffice_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
  ObsGnssroBendMetOfficeParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_H_
