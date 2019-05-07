/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/BackgroundCheck.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/GeoVaLs.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, BackgroundCheck> >
  makerBgChk_("Background Check");
// -----------------------------------------------------------------------------

BackgroundCheck::BackgroundCheck(ioda::ObsSpace & os, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(os), config_(config), abs_threshold_(-1.0), threshold_(-1.0), gv_(NULL),
    geovars_(preProcessWhere(config_)), flags_(*flags), obserr_(*obserr)
{
  oops::Log::trace() << "BackgroundCheck contructor starting" << std::endl;
  oops::Log::debug() << "BackgroundCheck: config = " << config << std::endl;
  ASSERT(flags);
  ASSERT(obserr);

  const float missing = util::missingValue(missing);
  threshold_ = config.getFloat("threshold", missing);
  abs_threshold_ = config.getFloat("absolute threshold", missing);
  ASSERT(abs_threshold_ != missing || threshold_ != missing);
  ASSERT(abs_threshold_ == missing || abs_threshold_ > 0.0);
  ASSERT(threshold_ == missing || threshold_ > 0.0);
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::trace() << "BackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void BackgroundCheck::priorFilter(const GeoVaLs & gv) const {
  gv_ = &gv;
}

// -----------------------------------------------------------------------------

void BackgroundCheck::postFilter(const ioda::ObsVector & hofx) const {
  oops::Log::trace() << "BackgroundCheck postFilter" << std::endl;

  const oops::Variables vars(config_.getStringVector("variables"));
  const oops::Variables observed(config_.getStringVector("observed"));
  const float missing = util::missingValue(missing);

  oops::Log::debug() << "BackgroundCheck flags: " << flags_;
  oops::Log::debug() << "BackgroundCheck obserr: " << obserr_;

  ioda::ObsDataVector<float> obs(obsdb_, vars, "ObsValue");
  ioda::ObsDataVector<float> bias(obsdb_, vars, "ObsBias", true);

// Select where the background check will apply
  std::vector<bool> apply = processWhere(obsdb_, *gv_, config_);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = observed.find(vars[jv]);

    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs] && flags_[iv][jobs] == 0) {
        size_t iobs = observed.size() * jobs + iv;
        ASSERT(obserr_[iv][jobs] != util::missingValue(obserr_[iv][jobs]));
        ASSERT(obs[jv][jobs] != util::missingValue(obs[jv][jobs]));
        ASSERT(hofx[iobs] != util::missingValue(hofx[iobs]));

//      Threshold for current observation
        float zz = std::numeric_limits<float>::max();
        if (abs_threshold_ != missing) zz = abs_threshold_;
        if (threshold_ != missing) zz = std::min(zz, threshold_ * obserr_[iv][jobs]);
        ASSERT(zz < std::numeric_limits<float>::max() && zz > 0.0);

//      Apply bias correction
        float yy = obs[jv][jobs] + bias[jv][jobs];

//      Check distance from background
        if (std::abs(static_cast<float>(hofx[iobs]) - yy) > zz) flags_[iv][jobs] = QCflags::fguess;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void BackgroundCheck::print(std::ostream & os) const {
  os << "BackgroundCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
