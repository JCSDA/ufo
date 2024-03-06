/*
 * (C) British Crown Copyright 2020 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_H_
#define UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/scatwind/NeutralMetOffice/ObsScatwindNeutralMetOffice.interface.h"
#include "ufo/operators/scatwind/NeutralMetOffice/ObsScatwindNeutralMetOfficeParameters.h"

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
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsScatwindNeutralMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> qcflags;
  static const std::string classname() {return "ufo::ObsScatwindNeutralMetOffice";}

  ObsScatwindNeutralMetOffice(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsScatwindNeutralMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperScatwindNeutralMetOffice_;}
  const int & toFortran() const {return keyOperScatwindNeutralMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperScatwindNeutralMetOffice_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  ObsScatwindNeutralMetOfficeParameters parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_SCATWIND_NEUTRALMETOFFICE_OBSSCATWINDNEUTRALMETOFFICE_H_
