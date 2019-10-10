/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_TIMEOPER_OBSTIMEOPER_H_
#define UFO_TIMEOPER_OBSTIMEOPER_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/Locations.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/timeoper/ObsTimeOper.interface.h"

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
  class Locations;
  class ObsDiagnostics;


// -----------------------------------------------------------------------------
/// TimeInterp observation operator class
class ObsTimeOper : public ObsOperatorBase,
                    private util::ObjectCounter<ObsTimeOper> {
 public:
  static const std::string classname() {return "ufo::ObsTimeOper";}

  ObsTimeOper(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsTimeOper();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

  Locations * locations(const util::DateTime &, const util::DateTime &) const override;

// Other
  const oops::Variables & variables() const override {return actualoperator_->variables();}

  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ObsOperatorBase> actualoperator_;
  F90hop keyTimeOper_;
  const ioda::ObsSpace& odb_;
  util::DateTime time_;
  int timeStencil_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_TIMEOPER_OBSTIMEOPER_H_
