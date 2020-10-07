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

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/Locations.interface.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

Locations::Locations(const eckit::mpi::Comm & comm) : comm_(comm) {
  int nobs = 0;
  ufo_locs_setup_f90(keyLoc_, nobs);
}

// -------------------------------------------------------------------------------------------------

Locations::Locations(const ioda::ObsSpace & odb) : comm_(odb.comm()) {
  ufo_locs_init_f90(keyLoc_, odb);
}

// -------------------------------------------------------------------------------------------------
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

Locations::Locations(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm) : comm_(comm) {
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  const util::DateTime bgn = util::DateTime(conf.getString("window begin"));
  const util::DateTime end = util::DateTime(conf.getString("window end"));

  ioda::ObsSpace obspace(obsconf, comm, bgn, end, oops::mpi::myself());
  const int nlocs = obspace.nlocs();

  std::vector<double> lats(nlocs);
  std::vector<double> lons(nlocs);
  obspace.get_db("MetaData", "latitude", lats);
  obspace.get_db("MetaData", "longitude", lons);

  ufo_locs_create_f90(keyLoc_, nlocs, obspace, &lats[0], &lons[0]);
}

// -------------------------------------------------------------------------------------------------

Locations::Locations(const Locations & other) : comm_(other.comm_) {
  ufo_locs_copy_f90(keyLoc_, other.toFortran());
}

// -------------------------------------------------------------------------------------------------
Locations & Locations::operator+=(const Locations & other) {
  F90locs otherKeyLoc_ = other.toFortran();
  ufo_locs_concatenate_f90(keyLoc_, otherKeyLoc_);
  return *this;
}

// -------------------------------------------------------------------------------------------------

Locations::~Locations() {
  ufo_locs_delete_f90(keyLoc_);
}

// -------------------------------------------------------------------------------------------------

int Locations::nobs() const {
  int nobs;
  ufo_locs_nobs_f90(keyLoc_, nobs);
  return nobs;
}

// -------------------------------------------------------------------------------------------------

void Locations::print(std::ostream & os) const {
  int nobs, indx, max_indx, i(0);
  ufo_locs_nobs_f90(keyLoc_, nobs);
  ufo_locs_indx_f90(keyLoc_, i, indx, max_indx);
  os << "Locations: " << nobs << " locations: "
                      << max_indx << " maximum indx:";

  // Write lat and lon to debug stream
  double lat, lon;

  for (int i=0; i < nobs; ++i) {
    ufo_locs_indx_f90(keyLoc_, i, indx, max_indx);
    ufo_locs_coords_f90(keyLoc_, i, lat, lon);

    oops::Log::debug() << "obs " << i << ": " << "gv index = " << indx
                       << std::setprecision(2) << std::fixed
                       << " lat = " << lat << ", lon = " << lon << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------

}  // namespace ufo
