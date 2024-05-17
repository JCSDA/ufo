/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_ORBITALANGLE_H_
#define UFO_PREDICTORS_ORBITALANGLE_H_

#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

enum class FourierTermType {
  SIN, COS
};

struct FourierTermTypeParameterTraitsHelper {
  typedef FourierTermType EnumType;
  static constexpr char enumTypeName[] = "FourierTermType";
  static constexpr util::NamedEnumerator<FourierTermType> namedValues[] = {
    { FourierTermType::SIN, "sin" },
    { FourierTermType::COS, "cos" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::FourierTermType> :
    public EnumParameterTraits<ufo::FourierTermTypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

// -----------------------------------------------------------------------------

/// Configuration parameters of the OrbitalAngle predictor.
class OrbitalAngleParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(OrbitalAngleParameters, PredictorParametersBase);

 public:
  /// Order of the Fourier term.
  oops::RequiredParameter<int> order{"order", this};
  /// Type of the Fourier term (either `sin` or `cos`).
  oops::RequiredParameter<FourierTermType> component{"component", this};
};

// -----------------------------------------------------------------------------
/**
 *This orbital angle predictor is used to fit residual errors as a function of satellite
 *orbital angle using a Fourier series. The data must contain this orbital angle time
 *series in the variable "MetaData/satelliteOrbitalAngle". Two member variables are 
 *used to store the order of the term in the series being calculated (order_) and the 
 *Fourier component (cos or sin), respectively. These are read from the yaml configuration
 *file.
 */

class OrbitalAngle : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef OrbitalAngleParameters Parameters_;

  OrbitalAngle(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  int order_;
  FourierTermType component_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_ORBITALANGLE_H_
