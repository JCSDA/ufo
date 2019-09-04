/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/BlackList.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------

BlackList::BlackList(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                     boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                     boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), data_(obsdb_), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    diagvars_(), flags_(*flags)
{
  oops::Log::debug() << "BlackList: config = " << config_ << std::endl;
  oops::Log::debug() << "BlackList: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

BlackList::~BlackList() {}

// -----------------------------------------------------------------------------

void BlackList::priorFilter(const GeoVaLs & gv) const {
  const size_t nobs = obsdb_.nlocs();
  const oops::Variables vars(config_);
  if (vars.size() == 0) {
    oops::Log::error() << "No variables will be filtered out in filter "
                       << config_ << std::endl;
    ABORT("No variables specified to be filtered out in filter");
  }
  const oops::Variables observed = obsdb_.obsvariables();

  data_.associate(gv);
  std::vector<bool> blacklisted = processWhere(config_, data_);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = observed.find(vars[jv]);
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      if (blacklisted[jobs] && flags_[iv][jobs] == 0) flags_[iv][jobs] = QCflags::black;
    }
  }
}

// -----------------------------------------------------------------------------

void BlackList::print(std::ostream & os) const {
  os << "BlackList: config = " << config_ << " , geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
