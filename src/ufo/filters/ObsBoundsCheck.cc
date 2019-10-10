/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsBoundsCheck.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBoundsCheck::ObsBoundsCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  if (config_.has("test variables")) {
    eckit::LocalConfiguration testvarconf(config_, "test variables");
    allvars_ += ufo::Variables(testvarconf);
  }
  oops::Log::debug() << "ObsBoundsCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsBoundsCheck::~ObsBoundsCheck() {}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::applyFilter(const std::vector<bool> & apply,
                                 std::vector<std::vector<bool>> & flagged) const {
  const float missing = util::missingValue(missing);

// Find which variables to apply filter on
  oops::Variables filtervars(config_);

// Find which variables are tested and the conditions
  ufo::Variables testvars;
// Use variables specified in test variables for testing, otherwise filter variables
  if (config_.has("test variables")) {
    eckit::LocalConfiguration testvarconf(config_, "test variables");
    testvars += ufo::Variables(testvarconf);
  } else {
    testvars += ufo::Variables(config_, "ObsValue");
  }
  const float vmin = config_.getFloat("minvalue", missing);
  const float vmax = config_.getFloat("maxvalue", missing);

// Sanity checks
  if (filtervars.size() == 0) {
    oops::Log::error() << "No variables will be filtered out in filter "
                       << config_ << std::endl;
    ABORT("No variables specified to be filtered out in filter");
  }
  if (filtervars.size() != testvars.size()) {
    oops::Log::error() << "Filter and test variables in Bounds Check have "
                       << "different sizes: " << filtervars.size() << " and "
                       << testvars.size() << std::endl;
    ABORT("Filter and test variables in Bounds Check have different sizes");
  }
  oops::Log::debug() << "ObsBoundsCheck: filtering " << filtervars << " with "
                     << testvars << std::endl;

// Find which variables are in flags/obserror
  oops::Variables observed = obsdb_.obsvariables();

// Loop over all variables to filter
  for (size_t jv = 0; jv < filtervars.size(); ++jv) {
//  find index of the filtered variable in flags/obserror
    size_t iv = observed.find(filtervars[jv]);
//  get test data for this variable
    std::vector<float> testdata;
    data_.get(testvars[jv], testdata);
//  apply the filter
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs] && flags_[iv][jobs] == 0) {
        ASSERT(testdata[jobs] != missing);
        if (vmin != missing && testdata[jobs] < vmin) flags_[iv][jobs] = QCflags::bounds;
        if (vmax != missing && testdata[jobs] > vmax) flags_[iv][jobs] = QCflags::bounds;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::print(std::ostream & os) const {
  os << "ObsBoundsCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
