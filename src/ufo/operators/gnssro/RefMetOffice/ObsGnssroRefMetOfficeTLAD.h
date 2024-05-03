/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICETLAD_H_
#define UFO_OPERATORS_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICETLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/gnssro/RefMetOffice/ObsGnssroRefMetOfficeTLAD.interface.h"
#include "ObsGnssroRefMetOfficeParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// GnssroRefMetOffice observation operator - TL/AD
// -----------------------------------------------------------------------------
class ObsGnssroRefMetOfficeTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGnssroRefMetOfficeTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsGnssroRefMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGnssroRefMetOfficeTLAD";}

  ObsGnssroRefMetOfficeTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGnssroRefMetOfficeTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGnssroRefMetOffice_;}
  const int & toFortran() const {return keyOperGnssroRefMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGnssroRefMetOffice_;
  std::unique_ptr<const oops::Variables> varin_;
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_GNSSRO_REFMETOFFICE_OBSGNSSROREFMETOFFICETLAD_H_
