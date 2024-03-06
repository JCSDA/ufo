/*
 * (C) British Crown Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICE_H_
#define UFO_OPERATORS_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICE_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/gnssro/RefMetOffice/ObsGnssroRefMetOffice.interface.h"

#include "ObsGnssroRefMetOfficeParameters.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// GnssroRefMetOffice observation operator
// -----------------------------------------------------------------------------
class ObsGnssroRefMetOffice : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroRefMetOffice> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsGnssroRefMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroRefMetOffice";}

  ObsGnssroRefMetOffice(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroRefMetOffice();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroRefMetOffice_;}
  const int & toFortran() const {return keyOperGnssroRefMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroRefMetOffice_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICE_H_
