/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsDomainErrCheck.h"

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
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ObsDomainErrCheck>>
  mkDomErLst_("DomainErr Check");
// -----------------------------------------------------------------------------

ObsDomainErrCheck::ObsDomainErrCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(obsdb), data_(obsdb), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    diagvars_(), flags_(*flags), obserr_(*obserr), parameter_(0.0)
{
  oops::Log::debug() << "ObsDomainErrCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "ObsDomainErrCheck: geovars = " << geovars_ << std::endl;
  ASSERT(obserr);

  const float missing = util::missingValue(missing);
  float parameter_ = config.getFloat("infltparameter", missing);
  ASSERT(parameter_ != missing);
}

// -----------------------------------------------------------------------------

ObsDomainErrCheck::~ObsDomainErrCheck() {}

// -----------------------------------------------------------------------------

void ObsDomainErrCheck::priorFilter(const GeoVaLs & gv) const {
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

// compute function
  std::vector<eckit::LocalConfiguration> masks;
  config_.get("where", masks);
  std::vector<float> values(nlocs);
  for (size_t jm = 0; jm < masks.size(); ++jm) {
//  Get variable and group
    const std::string vargrp(masks[jm].getString("variable"));
    std::string fvar, grp;
    std::string obgrp = "MetaData";
    splitVarGroup(vargrp, fvar, grp);
    if (fvar == "Scattering" && grp == "ObsFunction") {
      ioda::ObsDataVector<float> vals(obsdb_, fvar);
      ObsFunction obsdiag(fvar);
      obsdiag.compute(data_, vals);
      for (size_t jj = 0; jj < nlocs; ++jj) {
        values[jj] = vals[fvar][jj];
      }
    }
  }

  std::vector<bool> inside = processWhere(config_, data_);

  size_t count = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = observed.find(vars[jv]);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!inside[jobs] && flags_[iv][jobs] == 0) {
        flags_[iv][jobs] = QCflags::domain;
      } else {
        ASSERT(obserr_[iv][jobs] != util::missingValue(obserr_[iv][jobs]));
        ASSERT(obs[jv][jobs] != util::missingValue(obs[jv][jobs]));
        float bound = 2.5*obserr_[iv][jobs];
        float obserrinc = parameter_ * std::max((values[jobs]-9.0), 0.0) * obserr_[iv][jobs];
        obserrinc = std::max(obserr_[iv][jobs], bound);
        obserr_[iv][jobs] = sqrt(pow(obserr_[iv][jobs], 2) + pow(obserrinc, 2));
        ++count;
        }
      }
  }
  oops::Log::info() << "count=" << count << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDomainErrCheck::print(std::ostream & os) const {
  os << "ObsDomainErrCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
