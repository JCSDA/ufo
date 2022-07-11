/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsDomainCheck.h"

#include <vector>

#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDomainCheck::ObsDomainCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "ObsDomainCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsDomainCheck::~ObsDomainCheck() {}

// -----------------------------------------------------------------------------

void ObsDomainCheck::applyFilter(const std::vector<bool> & inside,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      flagged[jv][jobs] = !inside[jobs];
    }
  }
}

// -----------------------------------------------------------------------------

void ObsDomainCheck::print(std::ostream & os) const {
  os << "ObsDomainCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
