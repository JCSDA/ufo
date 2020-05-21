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
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

BlackList::BlackList(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                     boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                     boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "BlackList: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

BlackList::~BlackList() {}

// -----------------------------------------------------------------------------

void BlackList::applyFilter(const std::vector<bool> & apply,
                            const Variables & filtervars,
                            std::vector<std::vector<bool>> & flagged) const {
// Initialize map from filtervars to observed variables
  const oops::Variables observed = obsdb_.obsvariables();
  std::vector<size_t> filt2obs;
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    filt2obs.push_back(observed.find(filtervars.variable(jv).variable()));
  }

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      flagged[jv][jobs] = apply[jobs] && (*flags_)[filt2obs[jv]][jobs] == QCflags::pass;
    }
  }
}

// -----------------------------------------------------------------------------

void BlackList::print(std::ostream & os) const {
  os << "BlackList: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
