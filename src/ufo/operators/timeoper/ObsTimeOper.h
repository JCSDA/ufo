/*
 * (C) Copyright 2019 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_TIMEOPER_OBSTIMEOPER_H_
#define UFO_OPERATORS_TIMEOPER_OBSTIMEOPER_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/timeoper/ObsTimeOperParameters.h"

/// Forward declarations
namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsDiagnostics;


// -----------------------------------------------------------------------------
/// TimeInterp observation operator class
class ObsTimeOper : public ObsOperatorBase,
                    private util::ObjectCounter<ObsTimeOper> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsTimeOperParameters Parameters_;

  static const std::string classname() {return "ufo::ObsTimeOper";}

  ObsTimeOper(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsTimeOper();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  Locations_ locations() const override;

// Other
  const oops::Variables & requiredVars() const override {return actualoperator_->requiredVars();}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ObsOperatorBase> actualoperator_;
  const ioda::ObsSpace& odb_;
  std::vector<std::vector<float>> timeWeights_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_TIMEOPER_OBSTIMEOPER_H_
