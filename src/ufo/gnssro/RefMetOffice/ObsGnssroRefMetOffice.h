/*
 * (C) British Crown Copyright 2020 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICE_H_
#define UFO_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICE_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/gnssro/RefMetOffice/ObsGnssroRefMetOffice.interface.h"
#include "ufo/ObsOperatorBase.h"

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

/// GnssroRefMetOffice observation operator
class ObsGnssroRefMetOffice : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroRefMetOffice> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroRefMetOffice";}

  ObsGnssroRefMetOffice(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroRefMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroRefMetOffice_;}
  const int & toFortran() const {return keyOperGnssroRefMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroRefMetOffice_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICE_H_
