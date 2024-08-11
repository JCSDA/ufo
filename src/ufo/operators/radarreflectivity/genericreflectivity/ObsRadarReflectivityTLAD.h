/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITYTLAD_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITYTLAD_H_

#include <memory>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ObsRadarReflectivityParameters.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ReflectivityAlgorithmBase.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Radar reflectivity observation operator TL/AD code.
/// Please refer to the equivalent observation operator for further documentation.
class ObsRadarReflectivityTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsRadarReflectivityTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsRadarReflectivityParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() { return "ufo::ObsRadarReflectivityTLAD"; }

  ObsRadarReflectivityTLAD(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsRadarReflectivityTLAD() override;

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVarsTL_; }

  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Placeholder for required variables for the nonlinear operator.
  /// Ensures the reflectivity algorithm constructor has the correct signature.
  oops::Variables requiredVarsPlaceholder_;

  /// Required variables for the TL operator.
  oops::Variables requiredVarsTL_;

  /// Operator variables.
  oops::ObsVariables operatorVars_;

  /// Pointer to the reflectivity algorithm class.
  std::unique_ptr<ReflectivityAlgorithmBase> reflectivityAlgorithm_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITYTLAD_H_
