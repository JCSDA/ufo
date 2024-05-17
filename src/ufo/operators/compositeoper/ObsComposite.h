/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITE_H_
#define UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/compositeoper/ObsCompositeParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief A collection of observation operators used to simulate different variables.
///
/// Use this operator to split the list of simulated variables in the ObsSpace into groups,
/// each of which should be simulated with a different operator. For example, if the list of
/// simulated variables contains some upper-air variables that need to be vertically interpolated
/// and some surface variables that don't, you can arrange the `obs operator` section in the input
/// YAML file as in the following example:
///
///     obs operator:
///       name: Composite
///       components:
///       - name: VertInterp
///         variables:
///         - name: relativeHumidity
///           name: windNorthward
///           name: windEastward
///       - name: Identity
///         variables:
///         - name: surfacePressure
///
/// \note Only some operators (currently VertInterp and Identity) currently support the `variables`
/// option and thus can be used to simulate only a subset of variables.
class ObsComposite : public ObsOperatorBase,
                     private util::ObjectCounter<ObsComposite> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsCompositeParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsComposite";}

  ObsComposite(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsComposite() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

 private:
  const ioda::ObsSpace& odb_;
  std::vector<std::unique_ptr<ObsOperatorBase>> components_;
  oops::Variables requiredVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITE_H_
