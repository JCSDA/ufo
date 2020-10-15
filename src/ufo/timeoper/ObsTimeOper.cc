/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/timeoper/ObsTimeOper.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/timeoper/ObsTimeOperUtil.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsTimeOper> makerTimeOper_("TimeOperLinInterp");
// -----------------------------------------------------------------------------

ObsTimeOper::ObsTimeOper(const ioda::ObsSpace & odb,
                         const eckit::Configuration & config)
  : ObsOperatorBase(odb, config),
    actualoperator_(ObsOperatorFactory::create(
    odb, eckit::LocalConfiguration(config, "obs operator"))),
    odb_(odb), timeWeights_(timeWeightCreate(odb, config))
{
  oops::Log::trace() << "ObsTimeOper creating" << std::endl;

  util::DateTime windowBegin(odb_.windowStart());
  util::DateTime windowEnd(odb_.windowEnd());

  util::Duration windowSub;
  windowSub = util::Duration(config.getString("windowSub"));
  util::Duration window = windowEnd - windowBegin;

  if (window == windowSub) {
    ABORT("Time Interpolation of obs not implemented when assimilation window = subWindow");
  }
  oops::Log::trace() << "ObsTimeOper created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsTimeOper::~ObsTimeOper() {
  oops::Log::trace() << "ObsTimeOper destructed" << std::endl;
}


// -----------------------------------------------------------------------------

std::unique_ptr<Locations> ObsTimeOper::locations() const {
  oops::Log::trace() << "entered ObsOperatorTime::locations" << std::endl;

  std::unique_ptr<Locations> locs;

  locs = std::unique_ptr<Locations>(new Locations(odb_));
  *locs += *locs;

  // create concatenation of Locations class
  return locs;
}


void ObsTimeOper::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              ObsDiagnostics & ydiags) const {
  oops::Log::trace() << "ObsTimeOper: simulateObs entered" << std::endl;

  oops::Log::trace() << gv <<  std::endl;

  GeoVaLs gv1(odb_.comm()), gv2(odb_.comm());
  gv.split(gv1, gv2);

  oops::Log::trace() << gv1 << std::endl;
  oops::Log::trace() << gv2 << std::endl;

  gv1 *= timeWeights_[0];
  gv2 *= timeWeights_[1];
  gv1 += gv2;
  actualoperator_->simulateObs(gv1, ovec, ydiags);

  oops::Log::trace() << "ObsTimeOper: simulateObs exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsTimeOper::print(std::ostream & os) const {
  os << "ObsTimeOper::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo


