/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICETLAD_H_
#define UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICETLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/groundgnss/ZenithTotalDelayMetOffice/ObsGroundgnssMetOfficeParameters.h"
#include "ufo/operators/groundgnss/ZenithTotalDelayMetOffice/ObsGroundgnssMetOfficeTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// GroundgnssMetOffice observation operator
class ObsGroundgnssMetOfficeTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGroundgnssMetOfficeTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsGroundgnssMetOfficeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGroundgnssMetOfficeTLAD";}

  ObsGroundgnssMetOfficeTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGroundgnssMetOfficeTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGroundgnssMetOffice_;}
  const int & toFortran() const {return keyOperGroundgnssMetOffice_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGroundgnssMetOffice_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICETLAD_H_
