/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsBoundsCheck.h"

#include <algorithm>
#include <set>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBoundsCheck::ObsBoundsCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  if (config_.has("test variables") || config_.has("test functions")) {
    std::vector<eckit::LocalConfiguration> testvarconf;
    if (config_.has("test variables")) config_.get("test variables", testvarconf);
    if (config_.has("test functions")) config_.get("test functions", testvarconf);
    allvars_ += ufo::Variables(testvarconf);
  }
  oops::Log::debug() << "ObsBoundsCheck: config (constructor) = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsBoundsCheck::~ObsBoundsCheck() {}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::applyFilter(const std::vector<bool> & apply,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  const float missing = util::missingValue(missing);
  const oops::Variables observed = obsdb_.obsvariables();

// Find which variables are tested and the conditions
  ufo::Variables testvars;
// Use variables specified in test variables/functions for testing, otherwise filter variables
  if (config_.has("test variables") || config_.has("test functions")) {
    std::vector<eckit::LocalConfiguration> varconfs;
    if (config_.has("test variables")) config_.get("test variables", varconfs);
    if (config_.has("test functions")) config_.get("test functions", varconfs);
    testvars += ufo::Variables(varconfs);
  } else {
    testvars += ufo::Variables(filtervars, "ObsValue");
  }
  const float vmin = config_.getFloat("minvalue", missing);
  const float vmax = config_.getFloat("maxvalue", missing);

// Sanity checks
  if (filtervars.nvars() == 0) {
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

// Initialize map from filtervars to observed variables
  std::vector<size_t> filt2obs;
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    filt2obs.push_back(observed.find(filtervars.variable(jv).variable()));
  }

  if (config_.has("test functions")) {
    for (size_t iv = 0; iv < testvars.size(); ++iv) {
      ioda::ObsDataVector<float> testdata(obsdb_, testvars[iv].toOopsVariables(),
                                          "ObsFunction", false);
      data_.get(testvars[iv], testdata);

      std::vector<size_t> test_jv(filtervars[iv].size(), 0);
      if (testvars[iv].size() == filtervars[iv].size()) {
        std::iota(test_jv.begin(), test_jv.end(), 0);
      }

      // Loop over all variables to filter
      for (size_t jv = 0; jv < filtervars[iv].size(); ++jv) {
        for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
          if (apply[jobs] && (*flags_)[filt2obs[jv]][jobs] == QCflags::pass) {
            ASSERT(testdata[test_jv[jv]][jobs] != missing);
            if (vmin != missing && testdata[test_jv[jv]][jobs] < vmin) flagged[jv][jobs] = true;
            if (vmax != missing && testdata[test_jv[jv]][jobs] > vmax) flagged[jv][jobs] = true;
          }
        }
      }
    }
  } else {
    if (filtervars.nvars() != testvars.nvars()) {
      oops::Log::error() << "Filter and test variables in Bounds Check have "
                         << "different sizes: " << filtervars.nvars() << " and "
                         << testvars.nvars() << std::endl;
      ABORT("Filter and test variables in Bounds Check have different sizes");
    }
    // Loop over all variables to filter
    for (size_t jv = 0; jv < testvars.nvars(); ++jv) {
      //  get test data for this variable
      std::vector<float> testdata;
      data_.get(testvars.variable(jv), testdata);
      //  apply the filter
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (apply[jobs] && (*flags_)[filt2obs[jv]][jobs] == QCflags::pass) {
          ASSERT(testdata[jobs] != missing);
          if (vmin != missing && testdata[jobs] < vmin) flagged[jv][jobs] = true;
          if (vmax != missing && testdata[jobs] > vmax) flagged[jv][jobs] = true;
        }
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
