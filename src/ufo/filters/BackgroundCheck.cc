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

#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"

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
  allvars_ += Variables(filtervars_, "HofX");
  ASSERT(abs_threshold_ != "" || threshold_ != "");
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::trace() << "BackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void BackgroundCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "BackgroundCheck postFilter" << std::endl;
  const oops::Variables observed = obsdb_.obsvariables();
  const float missing = util::missingValue(missing);

  oops::Log::debug() << "BackgroundCheck obserr: " << obserr_;

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsVariables(), "ObsValue");
  ioda::ObsDataVector<float> bias(obsdb_, filtervars.toOopsVariables(), "ObsBias", false);

  Variables varhofx(filtervars_, "HofX");
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    size_t iv = observed.find(filtervars.variable(jv).variable());
//  H(x)
    std::vector<float> hofx;
    data_.get(varhofx.variable(jv), hofx);

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
        float zz = (thr[jobs] == std::numeric_limits<float>::max()) ? abs_thr[jobs] :
          std::min(abs_thr[jobs], thr[jobs] * obserr_[iv][jobs]);
        ASSERT(zz < std::numeric_limits<float>::max() && zz > 0.0);

//      Apply bias correction
        float yy = obs[jv][jobs] + bias[jv][jobs];

//      Check distance from background
        if (std::abs(static_cast<float>(hofx[jobs]) - yy) > zz) flagged[jv][jobs] = true;
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
