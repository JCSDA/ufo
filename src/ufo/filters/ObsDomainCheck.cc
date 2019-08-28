/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsDomainCheck.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/obsfunctions/ObsFunction.h"
#include "ufo/UfoTrait.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, ObsDomainCheck>>
  mkDomLst_("Domain Check");
// -----------------------------------------------------------------------------

ObsDomainCheck::ObsDomainCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    flags_(*flags)
{
  oops::Variables obsfcts(preProcessWhere(config_, "ObsFunction"));
  for (std::size_t ivar = 0; ivar < obsfcts.size(); ++ivar) {
    std::string var, grp;
    splitVarGroup(obsfcts[ivar], var, grp);
    ObsFunction function(var);
    geovars_ += function.requiredGeoVaLs();
  }
  oops::Log::debug() << "ObsDomainCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "ObsDomainCheck: geovars = " << geovars_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsDomainCheck::~ObsDomainCheck() {}

// -----------------------------------------------------------------------------

void ObsDomainCheck::priorFilter(const GeoVaLs & gv) const {
  const oops::Variables vars(config_);
  if (vars.size() == 0) {
    oops::Log::error() << "No variables will be filtered out in filter "
                       << config_ << std::endl;
    ABORT("No variables specified to be filtered out in filter");
  }
  const oops::Variables observed = obsdb_.obsvariables();

  std::vector<bool> inside = processWhere(obsdb_, gv, config_);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = observed.find(vars[jv]);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!inside[jobs] && flags_[iv][jobs] == 0) flags_[iv][jobs] = QCflags::domain;
    }
  }
}

// -----------------------------------------------------------------------------

void ObsDomainCheck::print(std::ostream & os) const {
  os << "ObsDomainCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
