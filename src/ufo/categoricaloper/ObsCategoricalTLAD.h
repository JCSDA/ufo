/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CATEGORICALOPER_OBSCATEGORICALTLAD_H_
#define UFO_CATEGORICALOPER_OBSCATEGORICALTLAD_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/categoricaloper/ObsCategoricalData.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
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

/// \brief Categorical observation operator TL/AD code.
/// Please refer to the Categorical observation operator for further documentation.
class ObsCategoricalTLAD : public LinearObsOperatorBase,
  private util::ObjectCounter<ObsCategoricalTLAD> {
 public:
  static const std::string classname() { return "ufo::ObsCategoricalTLAD"; }

  ObsCategoricalTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsCategoricalTLAD() override;

  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  const oops::Variables & requiredVars() const override { return data_.requiredVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the Categorical operator and TL/AD code.
  ObsCategoricalData<LinearObsOperatorBase, LinearObsOperatorFactory> data_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_CATEGORICALOPER_OBSCATEGORICALTLAD_H_
