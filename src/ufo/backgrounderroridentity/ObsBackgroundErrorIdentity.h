/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_H_
#define UFO_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

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

/// \brief Options controlling the ObsBackgroundErrorIdentity observation operator.
class ObsBackgroundErrorIdentityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBackgroundErrorIdentityParameters, Parameters)

 public:
  /// Name of the ObsOperator. Must be BackgroundErrorIdentity.
  ///
  /// TODO(wsmigaj): create an ObsOperatorParametersBase class, move this parameter there and
  /// derive ObsBackgroundErrorIdentityParameters from that class.
  oops::Parameter<std::string> name{"name", "", this};

  /// Simulated variables whose background errors may be calculated by this operator.
  /// If not specified, defaults to the list of all simulated variables in the ObsSpace.
  oops::OptionalParameter<std::vector<Variable>> variables{"variables", this};
};

/// \brief An observation operator calculating ObsDiagnostics representing single-level
/// background errors of simulated variables.
///
/// It should be used as a component of the `Composite` observation operator (with another
/// component handling the calculation of model equivalents of observation). It populates all
/// requested ObsDiagnostics called `<var>_background_error`, where `<var>` is the name of a
/// simulated variable, by copying the `<var>_background_error` GeoVaL at the observation
/// locations.
///
/// If the `variables` option is present, the operator does not calculate the background errors
/// of all simulated variables in the ObsSpace, but only those listed in the `variables` option.
///
/// See ObsBackgroundErrorIdentityParameters for the description of YAML configuration options
/// accepted by this operator.
///
/// Example configuration:
///
///     obs operator:
///       name: Composite
///       components:
///       # operator used to evaluate H(x)
///       - name: identity
///         vertical coordinate: geopotential_height  # coordinate used for obs value interpolation
///       # operator used to evaluate background errors
///       - name: BackgroundErrorIdentity
///         vertical coordinate: geopotential_height  # coordinate used for bg error interpolation
///         interpolation levels: background_error_geopotential_height
///
class ObsBackgroundErrorIdentity : public ObsOperatorBase,
                                     private util::ObjectCounter<ObsBackgroundErrorIdentity> {
 public:
  static const std::string classname() {return "ufo::ObsBackgroundErrorIdentity";}

  ObsBackgroundErrorIdentity(const ioda::ObsSpace &, const eckit::Configuration &);

  virtual ~ObsBackgroundErrorIdentity();

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

  const oops::Variables & requiredVars() const override;

  oops::Variables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

 private:
  const ioda::ObsSpace& odb_;
  ObsBackgroundErrorIdentityParameters parameters_;
  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_H_
