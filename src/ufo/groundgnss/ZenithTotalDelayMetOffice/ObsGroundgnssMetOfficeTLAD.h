/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICETLAD_H_
#define UFO_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICETLAD_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/groundgnss/ZenithTotalDelayMetOffice/ObsGroundgnssMetOfficeTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
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
/// GroundgnssMetOffice observation operator
class ObsGroundgnssMetOfficeTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsGroundgnssMetOfficeTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsGroundgnssMetOfficeTLAD";}

  ObsGroundgnssMetOfficeTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGroundgnssMetOfficeTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

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
#endif  // UFO_GROUNDGNSS_ZENITHTOTALDELAYMETOFFICE_OBSGROUNDGNSSMETOFFICETLAD_H_
