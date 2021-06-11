/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GROUNDGNSS_ZENITHTOTALDELAYROPP_OBSGROUNDGNSSROPP_H_
#define UFO_GROUNDGNSS_ZENITHTOTALDELAYROPP_OBSGROUNDGNSSROPP_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/groundgnss/ZenithTotalDelayROPP/ObsGroundgnssROPP.interface.h"
#include "ufo/ObsOperatorBase.h"

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

// -----------------------------------------------------------------------------

/// GroundgnssROPP observation operator
class ObsGroundgnssROPP : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGroundgnssROPP> {
 public:
  static const std::string classname() {return "ufo::ObsGroundgnssROPP";}

  ObsGroundgnssROPP(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGroundgnssROPP();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

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

#endif  // UFO_GROUNDGNSS_ZENITHTOTALDELAYROPP_OBSGROUNDGNSSROPP_H_
