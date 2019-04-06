/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <memory>
#include <vector>

#include "ufo/Locations.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/Locations.interface.h"

namespace ufo {

Locations::Locations(const ioda::ObsSpace & odb,
                     const util::DateTime & t1, const util::DateTime & t2) {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  ufo_locs_init_f90(keyLoc_, odb, &p1, &p2);
}

// -----------------------------------------------------------------------------
/*! UFO Locations Constructor with Configuration
 * 
 * \details This constructor can be used to generate user-specified
 * and/or random locations for use with interpolation or other tests
 *
 * To generate random locations, the relevant parameters specified in 
 * **StateTest.Locations** section of the config file are:
 * 
 * * **lats** user-specified latitudes (degrees)
 * * **lons** user-specified longitudes (degrees)
 * * **Nrandom** number of random locations desired
 * * **random_seed** (optional) random seed for reproducibility of results
 * 
 * \date May, 2018 Created (M. Miesch, JCSDA)
 *
 * \sa ufo::ufo_locs_create() ufo::ufo_loc_test() test::testStateInterpolation() 
 *
 */

Locations::Locations(const eckit::Configuration & conf) {
  std::vector<double> lats = conf.getDoubleVector("lats");
  std::vector<double> lons = conf.getDoubleVector("lons");

  ASSERT(lats.size() == lons.size());
  int nloc = lats.size();

  int rdist = 0;

  if (conf.has("Nrandom")) {
    int Nrandom = conf.getInt("Nrandom");

    unsigned int rseed;
    if (conf.has("random_seed")) {
      rseed = conf.getInt("random_seed");
    } else {
      rseed = std::time(0);
    }

    // random longitudes
    std::vector<double> lonrange;
    if (conf.has("lonrange")) {
      std::vector<double> config_lonrange = conf.getDoubleVector("lonrange");
      ASSERT(config_lonrange.size() == 2);
      lonrange.assign(begin(config_lonrange), end(config_lonrange));
    } else {
      lonrange.push_back(0.0);
      lonrange.push_back(360.0);
    }
    util::UniformDistribution<double> xx(Nrandom, lonrange[0], lonrange[1], rseed);
    for (size_t jj=0; jj < Nrandom; ++jj) lons.push_back(xx[jj]);

    // random latitudes
    std::vector<double> latrange;
    if (conf.has("latrange")) {
      std::vector<double> config_latrange = conf.getDoubleVector("latrange");
      ASSERT(config_latrange.size() == 2);
      latrange.assign(begin(config_latrange), end(config_latrange));
    } else {
      latrange.push_back(-90.0);
      latrange.push_back(90.0);
    }
    util::UniformDistribution<double> yy(Nrandom, latrange[0], latrange[1], rseed);
    for (size_t jj=0; jj < Nrandom; ++jj) lats.push_back(yy[jj]);

    nloc += Nrandom;

    if (conf.has("Rdist")) {
      rdist = conf.getInt("Rdist");
    }
  }

  ufo_locs_create_f90(keyLoc_, nloc, &lats[0], &lons[0], rdist);
}

// -----------------------------------------------------------------------------

Locations::~Locations() {
  ufo_locs_delete_f90(keyLoc_);
}

// -----------------------------------------------------------------------------

int Locations::nobs() const {
  int nobs;
  ufo_locs_nobs_f90(keyLoc_, nobs);
  return nobs;
}

// -----------------------------------------------------------------------------

void Locations::print(std::ostream & os) const {
  int nobs;
  ufo_locs_nobs_f90(keyLoc_, nobs);
  os << "Locations: " << nobs << " locations: ";

  // Write lat and lon to debug stream
  double lat, lon;

  for (int i=0; i < nobs; ++i) {
    ufo_locs_coords_f90(keyLoc_, i, lat, lon);
    oops::Log::debug() << "obs " << i << ": " << std::setprecision(2) << std::fixed
                       << "lat = " << lat << ", lon = " << lon << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

