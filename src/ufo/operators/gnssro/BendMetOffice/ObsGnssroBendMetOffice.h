/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_H_
#define UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/BendMetOffice/ObsGnssroBendMetOffice.interface.h"
#include "ufo/operators/gnssro/BendMetOffice/ObsGnssroBendMetOfficeParameters.h"

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
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsGnssroBendMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroBendMetOffice";}

  ObsGnssroBendMetOffice(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroBendMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBendMetOffice_;}
  const int & toFortran() const {return keyOperGnssroBendMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBendMetOffice_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICE_H_
