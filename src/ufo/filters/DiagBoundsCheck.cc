/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/DiagBoundsCheck.h"

#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------

DiagBoundsCheck::DiagBoundsCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), data_(obsdb_), config_(config), geovars_(preProcessWhere(config_, "GeoVaLs")),
    diagvars_(), flags_(*flags)
{
  if (config_.has("ObsDiags")) {
    std::vector<std::string> dvarstr = config_.getStringVector("ObsDiags");

    // Overwrite "variables" (accounts for presence of "channels")
    eckit::LocalConfiguration diagconfig(config_);
    diagconfig.set("variables", dvarstr);
    oops::Variables dvars(diagconfig);

    // Diagnostics setup
    for (size_t jv = 0; jv < dvars.size(); ++jv) {
      std::string var, grp;
      splitVarGroup(dvars[jv], var, grp);
      diagvars_.push_back(var);
    }
  } else {
    oops::Log::error() << "DiagBoundsCheck requires ObsDiags entry: config = "
                       << config_ << std::endl;
  }

  oops::Log::debug() << "DiagBoundsCheck: config = " << config_ << std::endl;
  oops::Log::debug() << "DiagBoundsCheck: geovars = " << geovars_ << std::endl;
  oops::Log::debug() << "DiagBoundsCheck: diagvars = " << diagvars_ << std::endl;
}

// -----------------------------------------------------------------------------

DiagBoundsCheck::~DiagBoundsCheck() {}

// -----------------------------------------------------------------------------

void DiagBoundsCheck::priorFilter(const GeoVaLs & gv) const {
  data_.associate(gv);
}
// -----------------------------------------------------------------------------
void DiagBoundsCheck::postFilter(const ioda::ObsVector & hofx,
                                  const ObsDiagnostics & diags) const {
  const float missing = util::missingValue(missing);

  oops::Variables vars(config_);
  if (vars.size() == 0) {
    oops::Log::error() << "No variables will be filtered out in filter "
                       << config_ << std::endl;
    ABORT("No variables specified to be filtered out in filter");
  }

  data_.associate(diags);

  oops::Variables simulated = hofx.varnames();

  const float vmin = config_.getFloat("minvalue", missing);
  const float vmax = config_.getFloat("maxvalue", missing);

// Select where the bounds check will apply
  std::vector<bool> apply = processWhere(config_, data_);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    size_t iv = simulated.find(vars[jv]);
    std::string dvar = diagvars_[jv] + "@ObsDiag";
    std::vector<float> diagvals = data_.get(dvar);
    for (size_t jobs = 0; jobs < hofx.nlocs(); ++jobs) {
      if (apply[jobs] && flags_[iv][jobs] == 0) {
        ASSERT(diagvals[jobs] != missing);
        if (vmin != missing && diagvals[jobs] < vmin) flags_[iv][jobs] = QCflags::bounds;
        if (vmax != missing && diagvals[jobs] > vmax) flags_[iv][jobs] = QCflags::bounds;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void DiagBoundsCheck::print(std::ostream & os) const {
  os << "DiagBoundsCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
