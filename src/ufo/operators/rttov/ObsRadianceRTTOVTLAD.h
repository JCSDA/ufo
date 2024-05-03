/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RTTOV_OBSRADIANCERTTOVTLAD_H_
#define UFO_OPERATORS_RTTOV_OBSRADIANCERTTOVTLAD_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/rttov/ObsRadianceRTTOVParameters.h"
#include "ufo/operators/rttov/ObsRadianceRTTOVTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadianceRTTOV (currently only temperature) observation for UFO.
class ObsRadianceRTTOVTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsRadianceRTTOVTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsRadianceRTTOVTLAD";}
  typedef ObsRadianceRTTOVParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  ObsRadianceRTTOVTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsRadianceRTTOVTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperRadianceRTTOV_;}
  const int & toFortran() const {return keyOperRadianceRTTOV_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperRadianceRTTOV_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RTTOV_OBSRADIANCERTTOVTLAD_H_
