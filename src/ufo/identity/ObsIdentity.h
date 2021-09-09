/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_IDENTITY_OBSIDENTITY_H_
#define UFO_IDENTITY_OBSIDENTITY_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"

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
  class ObsDiagnostics;

/// \brief Identity observation operator.
///
/// This observation operator transfers model values directly to the H(x) vector, after horizontal
/// interpolation has been performed, with no further processing.
/// For GeoVaLs with more than one vertical level, only the first entry in the GeoVaL is processed
/// in this way.
///
/// An example yaml configuration is:
///
///     obs operator:
///      name: Identity
///
/// This operator also accepts an optional `variables` parameter, which controls which ObsSpace
/// variables will be simulated. This option should only be set if this operator is used as a
/// component of the Composite operator. If `variables` is not set, the operator will simulate
/// all ObsSpace variables. Please see the documentation of the Composite operator for further
/// details.
///
/// The boolean parameter `level index 0 is closest to surface` can be set to `false`
/// for models whose level index 0 is furthest from the Earth's surface.
class ObsIdentity : public ObsOperatorBase,
  private util::ObjectCounter<ObsIdentity> {
 public:
  static const std::string classname() {return "ufo::ObsIdentity";}

  ObsIdentity(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsIdentity() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::Variables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

 private:
  /// Required variables.
  oops::Variables requiredVars_;

  /// Operator variables.
  oops::Variables operatorVars_;

  /// Indices of operator variables.
  std::vector<int> operatorVarIndices_;

  /// Level index 0 is closest to surface.
  bool levelIndexZeroAtSurface_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_IDENTITY_OBSIDENTITY_H_
