/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/Thinning.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

namespace ufo {

// -----------------------------------------------------------------------------

Thinning::Thinning(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                   std::shared_ptr<ioda::ObsDataVector<int> > flags,
                   std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "Thinning: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

Thinning::~Thinning() {}

// -----------------------------------------------------------------------------

void Thinning::applyFilter(const std::vector<bool> & apply,
                           const Variables & filtervars,
                           std::vector<std::vector<bool>> & flagged) const {
  // get local and global number of locations
  const size_t nlocs = obsdb_.nlocs();
  const size_t max_gindex = obsdb_.sourceNumLocs();

  // get global indices of the local locations
  const std::vector<std::size_t> & gindex = obsdb_.index();

  const float amount = parameters_.amount;

  // create random numbers for each observation based on some seed
  unsigned int random_seed = parameters_.randomSeed.value().value_or(std::time(0));
  random_seed += parameters_.member;

  util::UniformDistribution<float> rand(max_gindex, 0.0, 1.0, random_seed);

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      if ( apply[jobs] && rand[gindex[jobs]] < amount ) flagged[jv][jobs] = true;
    }
  }
}

// -----------------------------------------------------------------------------

void Thinning::print(std::ostream & os) const {
  os << "Thinning: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
