/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_BNDNBAM_OBSGNSSROBNDNBAM_H_
#define UFO_OPERATORS_GNSSRO_BNDNBAM_OBSGNSSROBNDNBAM_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/BndNBAM/ObsGnssroBndNBAM.interface.h"

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
// Gnssro BndNBAM observation operator
//   -- to reproduce exactly the operational (2019) NBAM
class ObsGnssroBndNBAM : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroBndNBAM> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroBndNBAM";}

  ObsGnssroBndNBAM(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroBndNBAM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBndNBAM_;}
  const int & toFortran() const {return keyOperGnssroBndNBAM_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBndNBAM_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_BNDNBAM_OBSGNSSROBNDNBAM_H_
