/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/BlackList.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

BlackList::BlackList(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                     std::shared_ptr<ioda::ObsDataVector<int> > flags,
                     std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "BlackList: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

BlackList::~BlackList() {}

// -----------------------------------------------------------------------------

void BlackList::applyFilter(const std::vector<bool> & apply,
                            const Variables & filtervars,
                            std::vector<std::vector<bool>> & flagged) const {
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      flagged[jv][jobs] = apply[jobs];
    }
  }
}

// -----------------------------------------------------------------------------

void BlackList::print(std::ostream & os) const {
  os << "BlackList: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
