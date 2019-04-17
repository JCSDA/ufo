/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/BackgroundCheck.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/processWhere.h"
#include "ufo/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, BackgroundCheck> >
  makerBgChk_("Background Check");
// -----------------------------------------------------------------------------

BackgroundCheck::BackgroundCheck(ioda::ObsSpace & os, const eckit::Configuration & config)
  : obsdb_(os), config_(config), threshold_(-1.0), geovars_()
{
  oops::Log::trace() << "BackgroundCheck contructor starting" << std::endl;
  oops::Log::debug() << "BackgroundCheck: config = " << config << std::endl;
  threshold_ = config.getDouble("threshold");
  ASSERT(threshold_ > 0.0);
  ASSERT(geovars_.size() == 0);
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::trace() << "BackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void BackgroundCheck::priorFilter(const GeoVaLs & gv) const {}

// -----------------------------------------------------------------------------

void BackgroundCheck::postFilter(const ioda::ObsVector & hofx) const {
  oops::Log::trace() << "BackgroundCheck postFilter" << std::endl;

  const oops::Variables vars(config_.getStringVector("observed"));
  const std::string qcgrp = config_.getString("QCname");
  const std::string obgrp = "ObsValue";
  const std::string ergrp = "ObsError";
  const float missing = util::missingValue(missing);

  ioda::ObsDataVector<double> obs(obsdb_, vars, obgrp);
  ioda::ObsDataVector<double> err(obsdb_, vars, ergrp);
  ioda::ObsDataVector<int> flags(obsdb_, vars, qcgrp);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
//  Select where the bounds check will apply
//    std::vector<bool> apply = processWhere(obsdb_, gv, bounds[jv]);
    std::vector<bool> apply(obsdb_.nlocs(), true);

    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs] && flags[jv][jobs] == 0) {
        ASSERT(obs[jv][jobs] != missing);
        size_t iobs = vars.size() * jobs + jv;
        if (std::abs(obs[jv][jobs] - hofx[iobs]) > threshold_ * err[jv][jobs]) {
          flags[jv][jobs] = QCflags::fguess;
        }
      }
    }
  }
  flags.save(qcgrp);
}

// -----------------------------------------------------------------------------

void BackgroundCheck::print(std::ostream & os) const {
  os << "BackgroundCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
