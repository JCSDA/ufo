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
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"
#include "ufo/filters/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------

Thinning::Thinning(const ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                   boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                   boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), config_(config), geovars_(), diagvars_(), flags_(*flags)
{
  oops::Log::debug() << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

Thinning::~Thinning() {}

// -----------------------------------------------------------------------------

void Thinning::preProcess() const {
  const size_t nobs = obsdb_.nlocs();
  const oops::Variables vars = obsdb_.obsvariables();
  const float thinning = config_.getFloat("amount");

  // create random numbers for each observation based on some seed
  unsigned int random_seed = config_.getInt("random_seed", std::time(0));
  util::UniformDistribution<float> rand(nobs, 0.0, 1.0, random_seed);

  for (size_t jv = 0; jv < vars.size(); ++jv) {
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      if ( rand[jobs] < thinning && flags_[jv][jobs] == 0) flags_[jv][jobs] = QCflags::thinned;
    }
  }
}

// -----------------------------------------------------------------------------

void Thinning::print(std::ostream & os) const {
  os << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
