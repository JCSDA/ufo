/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "Locations.h"
#include "Fortran.h"

namespace ufo {

// -----------------------------------------------------------------------------

Locations::Locations(const eckit::Configuration & conf) {
  std::vector<double> lats = conf.getDoubleVector("lats");
  std::vector<double> lons = conf.getDoubleVector("lons");
  ASSERT(lats.size() == lons.size());
  const int nloc = lats.size();
  ufo_locs_create_f90(keyLoc_, nloc, &lats[0], &lons[0]);
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
  os << "Locations: " << nobs << " locations";

  // Write lat and lon to debug stream
  double lat, lon;

  for (int i=0; i < nobs; ++i) {
    ufo_locs_coords_f90(keyLoc_,i,lat,lon);
    oops::Log::debug() << std::setprecision(2) << "lat = " << lat
		       << ", lon = " << lon << std::endl;
  }
  
}

// -----------------------------------------------------------------------------

}  // namespace UFO

