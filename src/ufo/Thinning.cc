/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/Thinning.h"

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
#include "ufo/QCflags.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------
static oops::FilterMaker<UfoTrait, oops::ObsFilter<UfoTrait, Thinning>> mkThinning_("Thinning");
// -----------------------------------------------------------------------------

Thinning::Thinning(ioda::ObsSpace & obsdb, const eckit::Configuration & config)
  : obsdb_(obsdb), config_(config), geovars_()
{
  oops::Log::debug() << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

Thinning::~Thinning() {}

// -----------------------------------------------------------------------------

void Thinning::priorFilter(const GeoVaLs & gv) const {
  const size_t nobs = obsdb_.nlocs();
  const std::string qcgrp = config_.getString("QCname");
  const oops::Variables vars(config_.getStringVector("observed"));
  const float thinning = config_.getFloat("amount");

  // create random numbers for each observation based on some seed
  unsigned int random_seed = config_.getInt("random_seed", std::time(0));
  util::UniformDistribution<float> rand(nobs, 0.0, 1.0, random_seed);

  ioda::ObsDataVector<int> flags(obsdb_, vars, qcgrp);
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      if ( rand[jobs] < thinning && flags[jv][jobs] == 0) flags[jv][jobs] = QCflags::thinned;
    }
  }
  flags.save(qcgrp);
}

// -----------------------------------------------------------------------------

void Thinning::print(std::ostream & os) const {
  os << "Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
