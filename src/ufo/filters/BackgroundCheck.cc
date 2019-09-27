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
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------

BackgroundCheck::BackgroundCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                 boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr),
    abs_threshold_(config_.getString("absolute threshold", "")),
    threshold_(config_.getString("threshold", ""))
{
  oops::Log::trace() << "BackgroundCheck contructor" << std::endl;
  const oops::Variables vars(config_);
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    allvars_ += vars[jv] + "@HofX";
  }
  ASSERT(abs_threshold_ != "" || threshold_ != "");
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::trace() << "BackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void BackgroundCheck::applyFilter(const std::vector<bool> & apply,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "BackgroundCheck postFilter" << std::endl;
  const oops::Variables vars(config_);
  if (vars.size() == 0) {
    oops::Log::error() << "No variables will be filtered out in filter "
                       << config_ << std::endl;
    ABORT("No variables specified to be filtered out in filter");
  }
  const oops::Variables observed = obsdb_.obsvariables();
  const float missing = util::missingValue(missing);

  oops::Log::debug() << "BackgroundCheck flags: " << flags_;
  oops::Log::debug() << "BackgroundCheck obserr: " << obserr_;

  ioda::ObsDataVector<float> obs(obsdb_, vars, "ObsValue");
  ioda::ObsDataVector<float> bias(obsdb_, vars, "ObsBias", false);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = observed.find(vars[jv]);

//  H(x)
    const std::string varhofx = vars[jv] + "@HofX";
    std::vector<float> hofx;
    data_.get(varhofx, hofx);

//  Threshold for current variable
    std::vector<float> abs_thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
    std::vector<float> thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
    if (abs_threshold_ != "") abs_thr = getScalarOrFilterData(abs_threshold_, data_);
    if (threshold_ != "")     thr     = getScalarOrFilterData(threshold_, data_);

    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs] && flags_[iv][jobs] == 0) {
        size_t iobs = observed.size() * jobs + iv;
        ASSERT(obserr_[iv][jobs] != util::missingValue(obserr_[iv][jobs]));
        ASSERT(obs[jv][jobs] != util::missingValue(obs[jv][jobs]));
        ASSERT(hofx[jobs] != util::missingValue(hofx[jobs]));

//      Threshold for current observation
        float zz = std::min(abs_thr[jobs], thr[jobs] * obserr_[iv][jobs]);
        ASSERT(zz < std::numeric_limits<float>::max() && zz > 0.0);

//      Apply bias correction
        float yy = obs[jv][jobs] + bias[jv][jobs];

//      Check distance from background
        if (std::abs(static_cast<float>(hofx[jobs]) - yy) > zz) flagged[iv][jobs] = true;
      }
    }
  }

// Apply action
  eckit::LocalConfiguration aconf;
  config_.get("action", aconf);
  aconf.set("flag", QCflags::fguess);
  FilterAction action(aconf);
  action.apply(vars, flagged, data_, flags_, obserr_);
}

// -----------------------------------------------------------------------------

void BackgroundCheck::print(std::ostream & os) const {
  os << "BackgroundCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
