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
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBoundsCheck::ObsBoundsCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
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

// Find which variables are tested and the conditions
  ufo::Variables testvars;
  std::string testname, testvar, testgrp;
// Use variables specified in test variables for testing, otherwise filter variables
  if (config_.has("test variables")) {
    std::vector<eckit::LocalConfiguration> varconfs;
    config_.get("test variables", varconfs);
    testvars += ufo::Variables(varconfs);
    testname = varconfs[0].getString("name");
    splitVarGroup(testname, testvar, testgrp);
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

  if (testgrp != "ObsFunction") {
    if (filtervars.nvars() != testvars.nvars()) {
      oops::Log::error() << "Filter and test variables in Bounds Check have "
                         << "different sizes: " << filtervars.nvars() << " and "
                         << testvars.nvars() << std::endl;
      ABORT("Filter and test variables in Bounds Check have different sizes");
    }
    oops::Log::debug() << "ObsBoundsCheck: filtering " << filtervars << " with "
                       << testvars << std::endl;

    // Loop over all variables to filter
    for (size_t jv = 0; jv < testvars.nvars(); ++jv) {
      //  get test data for this variable
      std::vector<float> testdata;
      data_.get(testvars.variable(jv), testdata);
      //  apply the filter
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (apply[jobs]) {
          ASSERT(testdata[jobs] != missing);
          if (vmin != missing && testdata[jobs] < vmin) flagged[jv][jobs] = true;
          if (vmax != missing && testdata[jobs] > vmax) flagged[jv][jobs] = true;
        }
      }
    }
  } else {
    std::vector<eckit::LocalConfiguration> varconfs;
    config_.get("test variables", varconfs);
    Variable testvar(varconfs[0]);
    ioda::ObsDataVector<float> testdata(obsdb_, testvar.toOopsVariables(), testvar.group(), false);
    data_.get(testvar, testdata);

    // if testdata is 1D variable, apply the same testdata to all filterdata
    // testdata_jv = {0, 0, 0, ..., 0} for all nvars
    std::vector<size_t> testdata_jv(filtervars.nvars(), 0);

    // if multiple variables ar in the testdata,  apply different testdatato different variables
    // testdata_jv = {0, 1, 2, ..., nvars-1}
    if (testvar.size() == filtervars.nvars()) {
      std::iota(testdata_jv.begin(), testdata_jv.end(), 0);
    }

    // Loop over all variables to filter
    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (apply[jobs]) {
          ASSERT(testdata[jv][jobs] != missing);
          if (vmin != missing && testdata[jv][jobs] < vmin) flagged[jv][jobs] = true;
          if (vmax != missing && testdata[jv][jobs] > vmax) flagged[jv][jobs] = true;
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
