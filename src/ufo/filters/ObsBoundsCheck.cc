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
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  if (config_.has("test variables")) {
    std::vector<eckit::LocalConfiguration> testvarconf;
    config_.get("test variables", testvarconf);
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
  if (config_.has("test variables")) {
    std::vector<eckit::LocalConfiguration> varconfs;
    config_.get("test variables", varconfs);
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

  oops::Log::debug() << "ObsBoundsCheck: filtering " << filtervars << " with "
                     << testvars << std::endl;

  oops::Log::debug() << "      sizes of each: " << filtervars.size() << " and "
                     << testvars.size() << std::endl;

// Initialize map from filtervars to observed variables
  std::vector<size_t> filt2obs;
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    filt2obs.push_back(observed.find(filtervars.variable(jv).variable()));
  }

  if (filtervars.size() == testvars.size()) {
    // Loop over all test variables to get data
    for (size_t iv = 0; iv < testvars.size(); ++iv) {
      const std::string grp = testvars[iv].group();
      ioda::ObsDataVector<float> testdata(obsdb_, testvars[iv].toOopsVariables());
      if (grp == "ObsFunction") {
        data_.get(testvars[iv], testdata);
      } else {
        for (size_t ii = 0; ii < testvars[iv].size(); ++ii) {
          size_t kv = ii + iv * testvars[iv].size();
          data_.get(testvars.variable(kv), testdata[ii]);
        }
      }

      std::vector<size_t> test_jv(filtervars[iv].size(), 0);
      if (testvars[iv].size() == filtervars[iv].size()) {
        std::iota(test_jv.begin(), test_jv.end(), 0);
      }

      // Loop over all variables to filter
      for (size_t jv = 0; jv < filtervars[iv].size(); ++jv) {
        for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
          if (apply[jobs] && (*flags_)[filt2obs[jv]][jobs] == QCflags::pass) {
            ASSERT(testdata[test_jv[jv]][jobs] != missing);
            size_t kv = jv + filtervars[iv].size() * iv;
            if (vmin != missing && testdata[test_jv[jv]][jobs] < vmin) flagged[kv][jobs] = true;
            if (vmax != missing && testdata[test_jv[jv]][jobs] > vmax) flagged[kv][jobs] = true;
          }
        }
      }
    }
  } else {
    int iv = 0;
    if (testvars.size() != 1) {
      oops::Log::error() << "When number filtervars not equal number of test vars, "
                         << "the latter can only be one." << config_ << std::endl;
      ABORT("ONLY one testvar when filtervars>1 because its usage is ambiguous otherwise");
    }

    ioda::ObsDataVector<float> testdata(obsdb_, testvars[iv].toOopsVariables());

    const std::string grp = testvars[iv].group();

    if (grp == "ObsFunction") {
      data_.get(testvars[iv], testdata);
    } else {
      data_.get(testvars.variable(iv), testdata);
    }

    for (size_t jv = 0; jv < filtervars.size(); ++jv) {
      oops::Log::debug() << "ObsBoundsCheck: testing filter var with index " << jv << std::endl;
      if (testvars[iv].size() != filtervars[jv].size()) {
        oops::Log::error() << "Dimension of filtervar, " << filtervars[jv].size()
                  << " does not equal testvar dimension, " << testvars[iv].size() << std::endl;
        ABORT("Aborting, sizes must be equivalent.");
      }
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (apply[jobs] && (*flags_)[filt2obs[jv]][jobs] == QCflags::pass) {
          ASSERT(testdata[iv][jobs] != missing);
          if (vmin != missing && testdata[iv][jobs] < vmin) flagged[jv][jobs] = true;
          if (vmax != missing && testdata[iv][jobs] > vmax) flagged[jv][jobs] = true;
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
