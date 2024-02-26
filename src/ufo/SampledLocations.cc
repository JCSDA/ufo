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

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Group.h"
#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

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
  lons_.resize(npaths);
  lats_.resize(npaths);
  for (size_t jj = 0; jj < npaths; ++jj) {
    lons_[jj] = lons[jj];
    lats_[jj] = lats[jj];
  }

  initializeObsGroup(npaths);

  // Set float_params that add safety in case of missing values and optimize the
  // potential future use of hdf5's "memory file" as the ObsGroup backend.
  const float floatMissing = util::missingValue<float>();
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;
  float_params.compressWithGZIP();
  float_params.setFillValue<float>(floatMissing);

  const ioda::Variable npathsVar = og_.vars["npaths"];
  og_.vars.createWithScales<float>("longitude", {npathsVar}, float_params).write(lons);
  og_.vars.createWithScales<float>("latitude", {npathsVar}, float_params).write(lats);

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

  initializeObsGroup(nlocs);

  // Set float_params that add safety in case of missing values and optimize the
  // potential future use of hdf5's "memory file" as the ObsGroup backend.
  const float floatMissing = util::missingValue<float>();
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;
  float_params.compressWithGZIP();
  float_params.setFillValue<float>(floatMissing);

  const ioda::Variable npathsVar = og_.vars["npaths"];
  std::vector<float> buffer(nlocs);
  lons_.resize(nlocs);
  lats_.resize(nlocs);

  obspace.get_db("MetaData", "longitude", buffer);
  for (size_t jj = 0; jj < nlocs; ++jj) lons_[jj] = buffer[jj];
  og_.vars.createWithScales<float>("longitude", {npathsVar}, float_params).write(buffer);

  obspace.get_db("MetaData", "latitude", buffer);
  for (size_t jj = 0; jj < nlocs; ++jj) lats_[jj] = buffer[jj];
  og_.vars.createWithScales<float>("latitude", {npathsVar}, float_params).write(buffer);

  times_.resize(nlocs);
  obspace.get_db("MetaData", "dateTime", times_);
}

// -------------------------------------------------------------------------------------------------

SampledLocations & SampledLocations::operator+=(
    const SampledLocations & other) {
  // Resize ObsGroup to new total size
  const ioda::Variable npathsVar = og_.vars["npaths"];
  const size_t npaths = npathsVar.getDimensions().dimsCur[0];
  const size_t other_npaths = other.og_.vars["npaths"].getDimensions().dimsCur[0];
  const ioda::Dimensions_t total_npaths = npaths + other_npaths;

  og_.resize({std::make_pair(npathsVar, total_npaths)});

  // Append variables from other's ObsGroup to end of local ObsGroup variables
  const std::vector<ioda::Dimensions_t> start(1, npaths);
  const std::vector<ioda::Dimensions_t> other_start(1, 0);
  const std::vector<ioda::Dimensions_t> counts(1, other_npaths);

  ioda::Selection feSelect;
  feSelect.extent({total_npaths}).select({ioda::SelectionOperator::SET, other_start, counts});
  ioda::Selection beSelect;
  beSelect.select({ioda::SelectionOperator::SET, start, counts});

  std::vector<float> buffer(other_npaths);
  other.og_.vars["longitude"].read<float>(gsl::make_span(buffer));
  og_.vars["longitude"].write<float>(buffer, feSelect, beSelect);

  other.og_.vars["latitude"].read<float>(gsl::make_span(buffer));
  og_.vars["latitude"].write<float>(buffer, feSelect, beSelect);

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
  ASSERT(og_.vars.exists("longitude"));
  const size_t npaths = size();
  std::vector<float> lons(npaths);
  og_.vars["longitude"].read<float>(gsl::make_span(lons));
  return lons;
}

// -------------------------------------------------------------------------------------------------

std::vector<float> SampledLocations::lats() const {
  ASSERT(og_.vars.exists("latitude"));
  const size_t npaths = size();
  std::vector<float> lats(npaths);
  og_.vars["latitude"].read<float>(gsl::make_span(lats));
  return lats;
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

void SampledLocations::initializeObsGroup(const size_t npaths) {
  ioda::Engines::BackendCreationParameters params;
  ioda::Group g = constructBackend(ioda::Engines::BackendNames::ObsStore, params);
  const ioda::NewDimensionScales_t dimScales{{
      ioda::NewDimensionScale<int>("npaths", npaths, ioda::Unlimited, npaths)}};
  og_ = ioda::ObsGroup::generate(g, dimScales);
}

// -------------------------------------------------------------------------------------------------

void SampledLocations::print(std::ostream & os) const {
  os << "Lat/lon/time paths: " << size() << " paths on this task " << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ufo
