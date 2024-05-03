/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICETLAD_H_
#define UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICETLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/gnssro/BendMetOffice/ObsGnssroBendMetOfficeParameters.h"
#include "ufo/operators/gnssro/BendMetOffice/ObsGnssroBendMetOfficeTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// GnssroBendMetOffice observation operator
class ObsGnssroBendMetOfficeTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGnssroBendMetOfficeTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsGnssroBendMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroBendMetOfficeTLAD";}

  ObsGnssroBendMetOfficeTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroBendMetOfficeTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroBendMetOffice_;}
  const int & toFortran() const {return keyOperGnssroBendMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroBendMetOffice_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_BENDMETOFFICE_OBSGNSSROBENDMETOFFICETLAD_H_
