/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTIONTLAD_H_
#define UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTIONTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/densityreduction/DensityReductionBase.h"
#include "ufo/operators/densityreduction/ObsDensityReductionParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Density reduction TL/AD observation operator class
class ObsDensityReductionTLAD : public LinearObsOperatorBase,
                                private util::ObjectCounter<ObsDensityReductionTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsDensityReductionParameters Parameters_;
  typedef typename ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() { return "ufo::ObsDensityReductionTLAD"; }

  ObsDensityReductionTLAD(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsDensityReductionTLAD() override;

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override;

 private:
  void print(std::ostream &) const override;

 private:
  const ioda::ObsSpace& odb_;
  std::unique_ptr<LinearObsOperatorBase> operator_;
  oops::Variables requiredVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTIONTLAD_H_
