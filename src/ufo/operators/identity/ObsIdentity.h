/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_IDENTITY_OBSIDENTITY_H_
#define UFO_OPERATORS_IDENTITY_OBSIDENTITY_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/identity/ObsIdentityParameters.h"

/// Forward declarations
namespace ioda {
class ObsSpace;
class ObsVector;
}  // namespace ioda

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

/// \brief Identity observation operator.
///
/// This observation operator transfers model values directly to the H(x)
/// vector, after horizontal interpolation has been performed, with no further
/// processing. For GeoVaLs with more than one vertical level, only the first
/// entry in the GeoVaL is processed in this way.
///
/// An example yaml configuration is:
///
///     obs operator:
///      name: Identity
///
/// The associated parameters class has details of all available options.
class ObsIdentity : public ObsOperatorBase,
                    private util::ObjectCounter<ObsIdentity> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() { return "ufo::ObsIdentity"; }
  typedef ObsIdentityParameters Parameters_;

  ObsIdentity(const ioda::ObsSpace &, const Parameters_ &);

  ~ObsIdentity() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override { return operatorVars_; }

 private:
  void print(std::ostream &) const override;

  /// Required variables.
  oops::Variables requiredVars_;

  /// Operator variables.
  oops::ObsVariables operatorVars_;

  /// Indices of operator variables.
  std::vector<int> operatorVarIndices_;

  /// Level index 0 is closest to surface.
  bool levelIndexZeroAtSurface_;
};

}  // namespace ufo
#endif  // UFO_OPERATORS_IDENTITY_OBSIDENTITY_H_
