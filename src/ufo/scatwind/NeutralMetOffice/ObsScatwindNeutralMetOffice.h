/*
 * (C) British Crown Copyright 2020 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_H_
#define UFO_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/scatwind/NeutralMetOffice/ObsScatwindNeutralMetOffice.interface.h"

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

/// ScatwindNeutralMetOffice observation operator
class ObsScatwindNeutralMetOffice : public ObsOperatorBase,
                      private util::ObjectCounter<ObsScatwindNeutralMetOffice> {
 public:
  static const std::string classname() {return "ufo::ObsScatwindNeutralMetOffice";}

  ObsScatwindNeutralMetOffice(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsScatwindNeutralMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperScatwindNeutralMetOffice_;}
  const int & toFortran() const {return keyOperScatwindNeutralMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperScatwindNeutralMetOffice_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_H_
