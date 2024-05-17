/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITETLAD_H_
#define UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITETLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/compositeoper/ObsCompositeParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Composite TL/AD observation operator class
class ObsCompositeTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsCompositeTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsCompositeParameters Parameters_;
  typedef typename ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() { return "ufo::ObsCompositeTLAD"; }

  ObsCompositeTLAD(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsCompositeTLAD() override;

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

 private:
  std::vector<std::unique_ptr<LinearObsOperatorBase>> components_;
  oops::Variables requiredVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_COMPOSITEOPER_OBSCOMPOSITETLAD_H_
