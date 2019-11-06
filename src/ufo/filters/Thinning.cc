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
  const size_t nobs = obsdb_.nlocs();
  const float thinning = config_.getFloat("amount");

  // create random numbers for each observation based on some seed
  unsigned int random_seed = config_.getInt("random_seed", std::time(0));
  util::UniformDistribution<float> rand(nobs, 0.0, 1.0, random_seed);

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      if ( apply[jobs] && rand[jobs] < thinning ) flagged[jv][jobs] = true;
    }
  }
}

// -----------------------------------------------------------------------------

void Thinning::print(std::ostream & os) const {
  os << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
