/*
 * (C) Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITY_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITY_H_

#include <memory>
#include <string>

#include "ioda/ObsDataVector.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ObsRadarReflectivityParameters.h"
#include "ufo/operators/radarreflectivity/genericreflectivity/ReflectivityAlgorithmBase.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Radar reflectivity observation operator.
///
/// This operator computes the model equivalent of reflectivity [dBZ].
/// The `algorithm.name` parameter controls which algorithm is used.
/// Each algorithm can implement its own set of parameters if desired.
///
/// An example yaml configuration that uses an algorithm called
/// `abc` is as follows:
///
/// - obs operator:
///     name: RadarReflectivity
///     algorithm:
///       name: abc
///
/// More concrete examples can be found in the header files
/// of each subclass of ReflectivityAlgorithmBase.
/// The ReflectivityAlgorithmBase header file contains further information
/// on how its subclasses should be created.
///
/// In addition, this operator can be used as part of a Composite operator.
/// Further details can be found in the ReflectivityAlgorithmBase code.
class ObsRadarReflectivity : public ObsOperatorBase,
  private util::ObjectCounter<ObsRadarReflectivity> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsRadarReflectivityParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsRadarReflectivity";}

  ObsRadarReflectivity(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsRadarReflectivity() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t&) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Required variables for the nonlinear operator.
  oops::Variables requiredVars_;

  /// Placeholder for required variables for the TL operator.
  /// Ensures the reflectivity algorithm constructor has the correct signature.
  oops::Variables requiredVarsTLPlaceholder_;

  /// Operator variables.
  oops::ObsVariables operatorVars_;

  /// Pointer to the reflectivity algorithm class.
  std::unique_ptr<ReflectivityAlgorithmBase> reflectivityAlgorithm_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARREFLECTIVITY_GENERICREFLECTIVITY_OBSRADARREFLECTIVITY_H_
