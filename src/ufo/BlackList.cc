/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/BlackList.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataRows.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/processWhere.h"
#include "ufo/QCflags.h"
#include "ufo/UfoTrait.h"
#include "ufo/utils/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, BlackList>> mkBlkLst_("BlackList");
// -----------------------------------------------------------------------------

BlackList::BlackList(ioda::ObsSpace & obsdb, const eckit::Configuration & config)
  : obsdb_(obsdb), config_(config), geovars_(preProcessWhere(config_))
{
  oops::Log::debug() << "BlackList: config = " << config_ << std::endl;
  oops::Log::debug() << "BlackList: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

BlackList::~BlackList() {}

// -----------------------------------------------------------------------------

void BlackList::priorFilter(const GeoVaLs & gv) const {
  const size_t nobs = obsdb_.nlocs();
  const std::string qcgrp = config_.getString("QCname");
  const std::vector<std::string> vars = config_.getStringVector("observed");

  std::vector<bool> blacklisted = processWhere(obsdb_, gv, config_);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    ioda::ObsDataRows<int> flags(obsdb_, vars[jv], qcgrp);
    int ii = 0;
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      if (blacklisted[jobs] && flags[jobs] == 0) flags[jobs] = QCflags::black;
      if (blacklisted[jobs]) ++ii;
    }
    flags.save();

    oops::Log::debug() << "BlackList: " << obsdb_.obsname() << " " << vars[jv]
                       << " rejected " << ii << " obs." << std::endl;
  }
}

// -----------------------------------------------------------------------------

void BlackList::print(std::ostream & os) const {
  os << "BlackList: config = " << config_ << " , geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
