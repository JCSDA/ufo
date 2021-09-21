/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_H_
#define UFO_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsDiagnostics;

/// \brief Options controlling the ObsBackgroundErrorVertInterp observation operator.
class ObsBackgroundErrorVertInterpParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsBackgroundErrorVertInterpParameters, ObsOperatorParametersBase)

 public:
  /// Simulated variables whose background errors may be calculated by this operator.
  /// If not specified, defaults to the list of all simulated variables in the ObsSpace.
  oops::OptionalParameter<std::vector<Variable>> variables{"variables", this};

  /// Name of the ufo variable (from the `MetaData` group) storing the vertical coordinate of
  /// observation locations.
  oops::RequiredParameter<std::string> observationVerticalCoordinate{
    "observation vertical coordinate", this};

  /// Name of the GeoVaL storing the interpolation levels of background errors.
  oops::RequiredParameter<std::string> verticalCoordinate{"vertical coordinate", this};
};

/// \brief An observation operator calculating ObsDiagnostics representing vertically interpolated
/// background errors of simulated variables.
///
/// It should be used as a component of the `Composite` observation operator (with another
/// component handling the calculation of model equivalents of observations). It populates all
/// requested ObsDiagnostics called `<var>_background_error`, where `<var>` is the name of a
/// simulated variable, by vertically interpolating the `<var>_background_error` GeoVaL at the
/// observation locations. Element (i, j) of this GeoVaL is interpreted as the background error
/// estimate of variable `<var>` at the ith observation location and the vertical position read from
/// the (i, j)th element of the GeoVaL specified in the `interpolation level` option of the
/// ObsOperator.
///
/// If the `variables` option is present, the operator does not calculate the background errors
/// of all simulated variables in the ObsSpace, but only those listed in the `variables` option.
///
/// See ObsBackgroundErrorVertInterpParameters for the description of YAML configuration options
/// accepted by this operator.
///
/// Example configuration:
///
///     obs operator:
///       name: Composite
///       components:
///       # operator used to evaluate H(x)
///       - name: VertInterp
///         vertical coordinate: geopotential_height  # coordinate used for obs value interpolation
///       # operator used to evaluate background errors
///       - name: BackgroundErrorVertInterp
///         # vertical coordinate of observation locations
///         observation vertical coordinate: geopotential_height
///         # GeoVaL storing interpolation levels of background errors
///         vertical coordinate: background_error_geopotential_height
///
class ObsBackgroundErrorVertInterp : public ObsOperatorBase,
                                     private util::ObjectCounter<ObsBackgroundErrorVertInterp> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsBackgroundErrorVertInterpParameters Parameters_;

  static const std::string classname() {return "ufo::ObsBackgroundErrorVertInterp";}

  ObsBackgroundErrorVertInterp(const ioda::ObsSpace &, const Parameters_ &);

  virtual ~ObsBackgroundErrorVertInterp();

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

  const oops::Variables & requiredVars() const override;

  oops::Variables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

 private:
  const ioda::ObsSpace& odb_;
  Parameters_ parameters_;
  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_H_
