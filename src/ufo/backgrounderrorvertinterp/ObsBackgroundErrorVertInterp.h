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
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsDiagnostics;

/// \brief Options controlling the ObsBackgroundErrorVertInterp observation operator.
class ObsBackgroundErrorVertInterpParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBackgroundErrorVertInterpParameters, Parameters)

 public:
  /// Name of the ObsOperator. Must be BackgroundErrorVertInterp.
  ///
  /// TODO(wsmigaj): create an ObsOperatorParametersBase class, move this parameter there and
  /// derive ObsBackgroundErrorVertInterpParameters from that class.
  oops::Parameter<std::string> name{"name", "", this};

  /// Simulated variables whose background errors may be calculated by this operator.
  /// If not specified, defaults to the list of all simulated variables in the ObsSpace.
  oops::OptionalParameter<std::vector<Variable>> variables{"variables", this};

  /// Name of the ufo variable (from the `MetaData` group) storing the vertical coordinate of
  /// observation locations.
  oops::Parameter<std::string> verticalCoordinate{"vertical coordinate", "air_pressure", this};

  /// Name of the GeoVaL storing the interpolation levels of background errors.
  ///
  /// If not specified, it is assumed to be the same as the value of the `vertical coordinate`
  /// option.
  oops::OptionalParameter<std::string> interpolationLevels{"interpolation levels", this};
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
///         vertical coordinate: geopotential_height  # coordinate used for bg error interpolation
///         interpolation levels: background_error_geopotential_height
///
class ObsBackgroundErrorVertInterp : public ObsOperatorBase,
                                     private util::ObjectCounter<ObsBackgroundErrorVertInterp> {
 public:
  static const std::string classname() {return "ufo::ObsBackgroundErrorVertInterp";}

  ObsBackgroundErrorVertInterp(const ioda::ObsSpace &, const eckit::Configuration &);

  virtual ~ObsBackgroundErrorVertInterp();

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

  const oops::Variables & requiredVars() const override;

  oops::Variables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

  std::string interpolationLevels() const;

 private:
  const ioda::ObsSpace& odb_;
  ObsBackgroundErrorVertInterpParameters parameters_;
  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_BACKGROUNDERRORVERTINTERP_OBSBACKGROUNDERRORVERTINTERP_H_
