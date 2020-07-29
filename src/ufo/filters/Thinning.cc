/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/Thinning.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

namespace ufo {

// -----------------------------------------------------------------------------

Thinning::Thinning(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                   boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                   boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

Thinning::~Thinning() {}

// -----------------------------------------------------------------------------

void Thinning::applyFilter(const std::vector<bool> & apply,
                           const Variables & filtervars,
                           std::vector<std::vector<bool>> & flagged) const {
  // get local and global number of locations
  const size_t nlocs = obsdb_.nlocs();
  const size_t gnlocs = obsdb_.gnlocs();

  // get global indices of the local locations
  const std::vector<std::size_t> & gindex = obsdb_.index();

  const float thinning = config_.getFloat("amount");

  // create random numbers for each observation based on some seed
  unsigned int random_seed = config_.getInt("random seed", std::time(0));
  util::UniformDistribution<float> rand(gnlocs, 0.0, 1.0, random_seed);

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      if ( apply[jobs] && rand[gindex[jobs]] < thinning ) flagged[jv][jobs] = true;
    }
  }
}

// -----------------------------------------------------------------------------

void Thinning::print(std::ostream & os) const {
  os << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
