/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/Locations.h"

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/Locations.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------

Locations::Locations() {
  int nobs = 0;
  ufo_locs_setup_f90(keyLoc_, nobs);
}

// -----------------------------------------------------------------------------

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
  const eckit::LocalConfiguration obsconf(conf, "ObsSpace");
  const util::DateTime bgn = util::DateTime(conf.getString("window_begin"));
  const util::DateTime end = util::DateTime(conf.getString("window_end"));

  ioda::ObsSpace obspace(obsconf, bgn, end);
  const int nlocs = obspace.nlocs();

  std::vector<double> lats(nlocs);
  std::vector<double> lons(nlocs);
  obspace.get_db("MetaData", "latitude", lats.size(), lats.data());
  obspace.get_db("MetaData", "longitude", lons.size(), lons.data());

  ufo_locs_create_f90(keyLoc_, nlocs, &lats[0], &lons[0]);
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

