/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_H_
#define UFO_OPERATORS_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

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
  class ObsDiagnostics;

/// \brief Options controlling the ObsBackgroundErrorIdentity observation operator.
class ObsBackgroundErrorIdentityParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsBackgroundErrorIdentityParameters, ObsOperatorParametersBase)

 public:
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
class ObsBackgroundErrorIdentity : public ObsOperatorBase,
                                     private util::ObjectCounter<ObsBackgroundErrorIdentity> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsBackgroundErrorIdentityParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsBackgroundErrorIdentity";}

  ObsBackgroundErrorIdentity(const ioda::ObsSpace &, const Parameters_ &);

  virtual ~ObsBackgroundErrorIdentity();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override;

  oops::ObsVariables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

 private:
  const ioda::ObsSpace& odb_;
  Parameters_ parameters_;
  oops::Variables requiredVars_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_BACKGROUNDERRORIDENTITY_OBSBACKGROUNDERRORIDENTITY_H_
