/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/Locations.h"

#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

Locations::Locations(const std::vector<float> & lons, const std::vector<float> & lats,
                     const std::vector<util::DateTime> & times,
                     std::shared_ptr<const ioda::Distribution> dist)
  : dist_(std::move(dist)), lons_(lons), lats_(lats), times_(times) {
  ASSERT(lons_.size() == lats_.size());
  ASSERT(lats_.size() == times_.size());
  oops::Log::trace() << "ufo::Locations::Locations done" << std::endl;
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
                     const eckit::mpi::Comm & comm) {
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  const util::DateTime bgn = util::DateTime(conf.getString("window begin"));
  const util::DateTime end = util::DateTime(conf.getString("window end"));
  ioda::ObsTopLevelParameters obsparams;
  obsparams.validateAndDeserialize(obsconf);

  ioda::ObsSpace obspace(obsparams, comm, bgn, end, oops::mpi::myself());
  const size_t nlocs = obspace.nlocs();
  dist_ = obspace.distribution();

  lats_.resize(nlocs);
  lons_.resize(nlocs);
  times_.resize(nlocs);

  obspace.get_db("MetaData", "latitude", lats_);
  obspace.get_db("MetaData", "longitude", lons_);
  obspace.get_db("MetaData", "datetime", times_);
}

// -------------------------------------------------------------------------------------------------
Locations & Locations::operator+=(const Locations & other) {
  lons_.insert(lons_.end(), other.lons_.begin(), other.lons_.end());
  lats_.insert(lats_.end(), other.lats_.begin(), other.lats_.end());
  times_.insert(times_.end(), other.times_.begin(), other.times_.end());
  return *this;
}

// -------------------------------------------------------------------------------------------------
std::vector<bool> Locations::isInTimeWindow(const util::DateTime & t1,
                                            const util::DateTime & t2) const {
  std::vector<bool> isIn(times_.size(), false);
  for (size_t ii = 0; ii < times_.size(); ++ii) {
    if (t1 < times_[ii] && times_[ii] <= t2) isIn[ii] = true;
  }
  return isIn;
}

// -------------------------------------------------------------------------------------------------

size_t Locations::size() const {
  return lats_.size();
}

// -------------------------------------------------------------------------------------------------

void Locations::print(std::ostream & os) const {
  os << "Lat/lon/time locations: " << size() << " locations on this task " << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ufo
