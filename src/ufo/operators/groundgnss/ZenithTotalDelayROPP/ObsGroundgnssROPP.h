/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYROPP_OBSGROUNDGNSSROPP_H_
#define UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYROPP_OBSGROUNDGNSSROPP_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/groundgnss/ZenithTotalDelayROPP/ObsGroundgnssROPP.interface.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

class GroundgnssROPPParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(GroundgnssROPPParameters, ObsOperatorParametersBase)
 public:
  // no parameters
};

// -----------------------------------------------------------------------------

/// GroundgnssROPP observation operator
class ObsGroundgnssROPP : public ObsOperatorBase,
                          private util::ObjectCounter<ObsGroundgnssROPP> {
 public:
  typedef GroundgnssROPPParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsGroundgnssROPP";}

  ObsGroundgnssROPP(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsGroundgnssROPP();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return *varin_;}

  int & toFortran() {return keyOperGroundgnssROPP_;}
  const int & toFortran() const {return keyOperGroundgnssROPP_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperGroundgnssROPP_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_GROUNDGNSS_ZENITHTOTALDELAYROPP_OBSGROUNDGNSSROPP_H_
