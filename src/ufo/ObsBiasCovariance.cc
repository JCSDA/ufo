/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cmath>
#include <memory>
#include <random>
#include <set>

#include "ufo/ObsBiasCovariance.h"

#include "ioda/ObsSpace.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasCovariance::ObsBiasCovariance(const ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : conf_(conf), variance_(0) {
  // Get the number of predictors
  std::vector<eckit::LocalConfiguration> confs;
  if (conf_.has("ObsBias.predictors")) {
    conf_.get("ObsBias.predictors", confs);
  }

  /// Get the jobs(channels)
  std::set<int> jobs;
  if (conf_.has("ObsBias.jobs")) {
    jobs = oops::parseIntSet(conf_.getString("ObsBias.jobs"));
  }

  for (std::size_t ii = 0; ii < confs.size()*jobs.size(); ++ii)
    variance_.push_back(1.0);
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::linearize(const ObsBias &) {
  oops::Log::warning() << "ObsBiasCovariance::linearize is not implmented" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::multiply(const ObsBiasIncrement & bx1, ObsBiasIncrement & bx2) const {
  bx2 = bx1;
  for (std::size_t ii = 0; ii < variance_.size(); ++ii) {
    bx2[ii] *= variance_[ii];
  }
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & bx1,
                                        ObsBiasIncrement & bx2) const {
  bx2 = bx1;
  for (std::size_t ii = 0; ii < variance_.size(); ++ii)
    bx2[ii] /= variance_[ii];
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  static util::NormalDistribution<double> dist(variance_.size(), 0.0, 1.0, 4);
  for (unsigned int jj = 0; jj < variance_.size(); ++jj) {
    dx[jj] = dist[jj] * std::sqrt(variance_[jj]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
