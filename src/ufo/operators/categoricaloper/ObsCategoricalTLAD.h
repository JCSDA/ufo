/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICALTLAD_H_
#define UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICALTLAD_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/categoricaloper/ObsCategoricalData.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Categorical observation operator TL/AD code.
/// Please refer to the Categorical observation operator for further documentation.
class ObsCategoricalTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsCategoricalTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsCategoricalParameters Parameters_;
  typedef typename ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() { return "ufo::ObsCategoricalTLAD"; }

  ObsCategoricalTLAD(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsCategoricalTLAD() override;

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return data_.requiredVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the Categorical operator and TL/AD code.
  ObsCategoricalData<LinearObsOperatorBase> data_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CATEGORICALOPER_OBSCATEGORICALTLAD_H_
