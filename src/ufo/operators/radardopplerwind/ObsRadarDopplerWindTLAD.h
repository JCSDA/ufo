/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDTLAD_H_
#define UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDTLAD_H_

#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/radardopplerwind/ObsRadarDopplerWindData.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Radar Doppler wind observation operator TL/AD code.
/// Please refer to the equivalent observation operator for further documentation.
class ObsRadarDopplerWindTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsRadarDopplerWindTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsRadarDopplerWindParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() { return "ufo::ObsRadarDopplerWindTLAD"; }

  ObsRadarDopplerWindTLAD(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsRadarDopplerWindTLAD() override;

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return data_.requiredVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the RadarDopplerWind operator and TL/AD code.
  ObsRadarDopplerWindData data_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARDOPPLERWIND_OBSRADARDOPPLERWINDTLAD_H_
