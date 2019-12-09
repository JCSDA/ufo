/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_TIMEOPER_OBSTIMEOPERTLAD_H_
#define UFO_TIMEOPER_OBSTIMEOPERTLAD_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/timeoper/ObsTimeOperUtil.h"


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
  class ObsBias;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// TimeOper observation operator
class ObsTimeOperTLAD : public LinearObsOperatorBase,
                        private util::ObjectCounter<ObsTimeOperTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsTimeOperTLAD";}

  ObsTimeOperTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsTimeOperTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & variables() const override {return actualoperator_->variables();}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<LinearObsOperatorBase> actualoperator_;
  const ioda::ObsSpace& odb_;
  std::vector<std::vector<float>> timeWeights_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_TIMEOPER_OBSTIMEOPERTLAD_H_
