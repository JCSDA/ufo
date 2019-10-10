/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/timeoper/ObsTimeOper.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsTimeOper> makerTimeOpInterp_("TimeOpInterp");
// -----------------------------------------------------------------------------

ObsTimeOper::ObsTimeOper(const ioda::ObsSpace & odb,
                         const eckit::Configuration & config)
  : ObsOperatorBase(odb, config),
    actualoperator_(ObsOperatorFactory::create(
    odb, eckit::LocalConfiguration(config, "ObsOperator"))),
    odb_(odb), timeStencil_(2)
{
  oops::Log::trace() << "ObsTimeOper creating" << std::endl;

  const eckit::Configuration * configc = &config;

  ufo_timeoper_setup_f90(keyTimeOper_, &configc, odb_.nlocs(), timeStencil_);

  util::DateTime windowBegin(odb_.windowStart());
  util::DateTime windowEnd(odb_.windowEnd());
  util::Duration windowSub;
  windowSub = util::Duration(config.getString("windowSub"));

  util::DateTime t0(windowBegin);
  util::DateTime t3, stateTime;

  for (util::DateTime stateTime = windowBegin;
       stateTime <= windowEnd; stateTime = stateTime + windowSub) {
    if (stateTime == windowBegin) {
      t0 = windowBegin;
      t3 = t0 + windowSub;
    } else if (stateTime == windowEnd) {
      t0 = windowEnd - windowSub;
      t3 = windowEnd;
    } else {
      t0 = stateTime - windowSub;
      t3 = stateTime + windowSub;
    }
    const util::DateTime * p0 = &t0;
    const util::DateTime * p3 = &t3;
    const util::DateTime * st = &stateTime;
    const util::Duration * ws = &windowSub;

    oops::Log::debug() << " t0 = " << t0.toString() << std::endl;
    oops::Log::debug() << " t3 = " << t3.toString() << std::endl;
    oops::Log::debug() << " stateTime = " << stateTime.toString() << std::endl;
    oops::Log::debug() << " windowSub = " << windowSub.toString() << std::endl;

    ufo_timeoper_set_timeweight_f90(keyTimeOper_, &configc,
                                    odb_, &p0, &p3, &st, &ws);
  }

  oops::Log::trace() << "ObsTimeOper created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsTimeOper::~ObsTimeOper() {
  oops::Log::trace() << "ObsTimeOper destructed" << std::endl;
}


// -----------------------------------------------------------------------------

Locations * ObsTimeOper::locations(const util::DateTime & t1,
                                   const util::DateTime & t2) const {
  oops::Log::trace() << "entered ObsOperatorTime::locations" << std::endl;

  Locations * locs = new Locations();
  int keylocs = locs->toFortran();

  util::DateTime t0, t3, stateTime;
  util::Duration initial_dt, windowSub;

  util::DateTime windowBegin(odb_.windowStart());
  util::DateTime windowEnd(odb_.windowEnd());

  initial_dt = t2 - t1;

  if ((t1 == windowBegin) && (t2 == windowEnd)) {
    windowSub = initial_dt;
    t0 = t1;
    t3 = t2;
    stateTime = t1 + initial_dt/2;
  } else {
    if ((t1 == windowBegin) || (t2 == windowEnd))
      windowSub = initial_dt * 2;
    else
      windowSub = initial_dt;

    t0 = t1 - (windowSub/2) * (2 * (timeStencil_) - 3);
    t3 = t2 + (windowSub/2) * (2 * (timeStencil_) - 3);
    stateTime = t1 + initial_dt/2;

    if (t1 == windowBegin) {
      t0 = windowBegin;
      stateTime = windowBegin;
    }
    if (t2 == windowEnd) {
      t3 = windowEnd;
      stateTime = windowEnd;
    }
  }
  oops::Log::debug() << " timeStencil_ = " << timeStencil_ << std::endl;
  oops::Log::debug() << " windowBegin = " << windowBegin.toString() << std::endl;
  oops::Log::debug() << " windowEnd = " << windowEnd.toString() << std::endl;
  oops::Log::debug() << " t0 = " << t0.toString() << std::endl;
  oops::Log::debug() << " t1 = " << t1.toString() << std::endl;
  oops::Log::debug() << " t2 = " << t2.toString() << std::endl;
  oops::Log::debug() << " t3 = " << t3.toString() << std::endl;
  oops::Log::debug() << " stateTime = " << stateTime.toString() << std::endl;

  const util::DateTime * p0 = &t0;
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  const util::DateTime * p3 = &t3;
  const util::DateTime * st = &stateTime;

  ufo_timeoper_locs_init_f90(keyTimeOper_, keylocs, odb_, &p0, &p1, &p2, &p3, &st);
  return locs;
}
// -----------------------------------------------------------------------------

void ObsTimeOper::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics & ydiags) const {
  oops::Log::trace() << "ObsTimeOper: simulateObs entered" << std::endl;

  oops::Log::debug() << gv;

  ufo_timeoper_simobs_f90(keyTimeOper_, gv.toFortran(), odb_);

  oops::Log::debug() << gv;

  actualoperator_->simulateObs(gv, ovec, ydiags);

  oops::Log::debug() << gv;

  oops::Log::trace() << "ObsTimeOper: simulateObs exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOper::print(std::ostream & os) const {
  os << "ObsTimeOper::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo


