/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/SampledLocations.h"

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

SampledLocations::SampledLocations(
    const std::vector<float> & lons, const std::vector<float> & lats,
    const std::vector<util::DateTime> & times, std::shared_ptr<const ioda::Distribution> dist,
    std::vector<util::Range<size_t>> pathsGroupedByLocation)
  : dist_(std::move(dist)), times_(std::move(times)), lons_(), lats_(),
    pathsGroupedByLocation_(std::move(pathsGroupedByLocation)) {
  oops::Log::trace() << "ufo::SampledLocations::SampledLocations start"
                     << std::endl;
  const size_t npaths = times_.size();
  ASSERT(npaths == lons.size());
  ASSERT(npaths == lats.size());
  // assign() rather than resize()+fill: every element is written from the input, so no
  // element is ever left unset.
  lons_.assign(lons.begin(), lons.end());
  lats_.assign(lats.begin(), lats.end());

  oops::Log::trace() << "ufo::SampledLocations::SampledLocations done"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------
/*! UFO SampledLocations constructor with Configuration
 *
 * \details This constructor can be used to generate user-specified
 * and/or random paths for use with interpolation or other tests
 *
 * To generate random paths, the relevant parameters specified in
 * **StateTest.SampledLocations** section of the config file are:
 *
 * * **lats** user-specified latitudes (degrees)
 * * **lons** user-specified longitudes (degrees)
 * * **Nrandom** number of random paths desired
 * * **random_seed** (optional) random seed for reproducibility of results
 *
 * \date May, 2018 Created (M. Miesch, JCSDA)
 *
 * \sa test::testStateInterpolation()
 *
 */

SampledLocations::SampledLocations(const eckit::Configuration & conf, const eckit::mpi::Comm & comm)
  : dist_(), times_(), lons_(), lats_() {
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const ioda::ObsSpace obspace(obsconf, comm, timeWindow, oops::mpi::myself());
  const size_t nlocs = obspace.nlocs();
  dist_ = obspace.distribution();

  std::vector<float> buffer(nlocs);

  obspace.get_db("MetaData", "longitude", buffer);
  lons_.assign(buffer.begin(), buffer.end());

  obspace.get_db("MetaData", "latitude", buffer);
  lats_.assign(buffer.begin(), buffer.end());

  times_.resize(nlocs);
  obspace.get_db("MetaData", "dateTime", times_);
}

// -------------------------------------------------------------------------------------------------

SampledLocations & SampledLocations::operator+=(
    const SampledLocations & other) {
  times_.insert(times_.end(), other.times_.begin(), other.times_.end());
  lats_.insert(lats_.end(), other.lats_.begin(), other.lats_.end());
  lons_.insert(lons_.end(), other.lons_.begin(), other.lons_.end());

  return *this;
}

// -------------------------------------------------------------------------------------------------

std::vector<bool> SampledLocations::isInTimeWindow(const util::DateTime & t1,
                                                              const util::DateTime & t2) const {
  std::vector<bool> isIn(times_.size(), false);
  for (size_t ii = 0; ii < times_.size(); ++ii) {
    if (t1 < times_[ii] && times_[ii] <= t2) isIn[ii] = true;
  }
  return isIn;
}

// -------------------------------------------------------------------------------------------------

size_t SampledLocations::size() const {
  return times_.size();
}

// -------------------------------------------------------------------------------------------------

std::vector<float> SampledLocations::lons() const {
  return std::vector<float>(lons_.begin(), lons_.end());
}

// -------------------------------------------------------------------------------------------------

std::vector<float> SampledLocations::lats() const {
  return std::vector<float>(lats_.begin(), lats_.end());
}

// -------------------------------------------------------------------------------------------------

size_t SampledLocations::nlocs() const {
  if (pathsGroupedByLocation_.empty())
    return size();
  else
    return pathsGroupedByLocation_.size();
}

// -------------------------------------------------------------------------------------------------

bool SampledLocations::areLocationsSampledOnceAndInOrder() const {
  if (pathsGroupedByLocation_.empty())
    return true;

  for (size_t loc = 0; loc < pathsGroupedByLocation_.size(); ++loc)
    if (pathsGroupedByLocation_[loc].begin != loc || pathsGroupedByLocation_[loc].end != loc + 1)
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

void SampledLocations::print(std::ostream & os) const {
  os << "Lat/lon/time paths: " << size() << " paths on this task " << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ufo
