/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsErrInflationCheck.h"

#include <algorithm>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ObsErrInflationCheck>>
  mkDomErLst_("ObsErrInflation Check");
// -----------------------------------------------------------------------------

ObsErrInflationCheck::ObsErrInflationCheck(
                               ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(obsdb), data_(obsdb), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    diagvars_(), flags_(*flags), obserr_(*obserr)
{
  oops::Log::debug() << "ObsErrInflationCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "ObsErrInflationCheck: geovars = " << geovars_ << std::endl;

  const float missing = util::missingValue(missing);
}

// -----------------------------------------------------------------------------

ObsErrInflationCheck::~ObsErrInflationCheck() {}

// -----------------------------------------------------------------------------

void ObsErrInflationCheck::priorFilter(const GeoVaLs & gv) const {
  data_.associate(gv);
  const oops::Variables vars(config_);
  if (vars.size() == 0) {
    oops::Log::error() << "No variables will be filtered out in filter "
                       << config_ << std::endl;
    ABORT("No variables specified to be filtered out in filter");
  }
  const oops::Variables observed = obsdb_.obsvariables();

  ioda::ObsDataVector<float> obs(obsdb_, vars, "ObsValue");
  size_t nlocs = obsdb_.nlocs();

// Get masks from yaml
  std::vector<eckit::LocalConfiguration> masks;
  config_.get("where", masks);
  std::vector<float> values(nlocs);
  for (size_t jm = 0; jm < masks.size(); ++jm) {
//    Get variable and group first
      const std::string vargrp(masks[jm].getString("variable"));
      std::string var, grp;
      std::string obgrp = "MetaData";
      splitVarGroup(vargrp, var, grp);
      oops::Log::debug() << "vargrp = " << vargrp << std::endl;
      oops::Log::debug() << "var = " << var << std::endl;
      oops::Log::debug() << "grp = " << grp << std::endl;

//    Compute function
      if (var == "ErrInflationFactor" && grp == "ObsFunction") {
         ioda::ObsDataVector<float> vals(obsdb_, var);
         ObsFunction obsfct(var);
         obsfct.compute(data_, vals);
         for (size_t jj = 0; jj < nlocs; ++jj) {
             values[jj] = vals[var][jj];
             oops::Log::debug() << "ObsErrInflationCheck: computed values = " << values[jj]
                                << std::endl;
         }
      }
}

  std::vector<bool> inside = processWhere(config_, data_);

// Apply inflation factor to observation error
  for (size_t jv = 0; jv < vars.size(); ++jv) {
     size_t iv = observed.find(vars[jv]);
     for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (!inside[jobs] && flags_[iv][jobs] == 0) {
           flags_[iv][jobs] = QCflags::domain;
        }
        float obserr0 = obserr_[iv][jobs];
        obserr_[iv][jobs] = values[jobs] * obserr_[iv][jobs];
        oops::Log::debug() << "ObsErrInflationCheck: unadjusted obserr = " << obserr0 <<
                              "adjusted obserr = " <<  obserr_[iv][jobs] << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

void ObsErrInflationCheck::print(std::ostream & os) const {
  os << "ObsErrInflationCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
