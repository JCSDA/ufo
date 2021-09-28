/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/Locations.h"

#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "ioda/Engines/Factory.h"
#include "ioda/Group.h"
#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

Locations::Locations(const std::vector<float> & lons, const std::vector<float> & lats,
                     const std::vector<util::DateTime> & times,
                     std::shared_ptr<const ioda::Distribution> dist)
  : dist_(std::move(dist)), times_(std::move(times)) {
  const size_t nlocs = times_.size();
  ASSERT(nlocs == lons.size());
  ASSERT(nlocs == lats.size());

  initializeObsGroup(nlocs);

  // Set float_params that add safety in case of missing values and optimize the
  // potential future use of hdf5's "memory file" as the ObsGroup backend.
  const float floatMissing = util::missingValue(floatMissing);
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;
  float_params.compressWithGZIP();
  float_params.setFillValue<float>(floatMissing);

  const ioda::Variable nlocsVar = og_.vars["nlocs"];
  og_.vars.createWithScales<float>("longitude", {nlocsVar}, float_params).write(lons);
  og_.vars.createWithScales<float>("latitude", {nlocsVar}, float_params).write(lats);

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

  const ioda::ObsSpace obspace(obsparams, comm, bgn, end, oops::mpi::myself());
  const size_t nlocs = obspace.nlocs();
  dist_ = obspace.distribution();

  initializeObsGroup(nlocs);

  // Set float_params that add safety in case of missing values and optimize the
  // potential future use of hdf5's "memory file" as the ObsGroup backend.
  const float floatMissing = util::missingValue(floatMissing);
  ioda::VariableCreationParameters float_params;
  float_params.chunk = true;
  float_params.compressWithGZIP();
  float_params.setFillValue<float>(floatMissing);

  const ioda::Variable nlocsVar = og_.vars["nlocs"];
  std::vector<float> buffer(nlocs);
  obspace.get_db("MetaData", "longitude", buffer);
  og_.vars.createWithScales<float>("longitude", {nlocsVar}, float_params).write(buffer);
  obspace.get_db("MetaData", "latitude", buffer);
  og_.vars.createWithScales<float>("latitude", {nlocsVar}, float_params).write(buffer);

  times_.resize(nlocs);
  obspace.get_db("MetaData", "datetime", times_);
}

// -------------------------------------------------------------------------------------------------
Locations & Locations::operator+=(const Locations & other) {
  // Resize ObsGroup to new total size
  const ioda::Variable nlocsVar = og_.vars["nlocs"];
  const size_t nlocs = nlocsVar.getDimensions().dimsCur[0];
  const size_t other_nlocs = other.og_.vars["nlocs"].getDimensions().dimsCur[0];
  const ioda::Dimensions_t total_nlocs = nlocs + other_nlocs;

  og_.resize({std::make_pair(nlocsVar, total_nlocs)});

  // Append variables from other's ObsGroup to end of local ObsGroup variables
  const std::vector<ioda::Dimensions_t> start(1, nlocs);
  const std::vector<ioda::Dimensions_t> other_start(1, 0);
  const std::vector<ioda::Dimensions_t> counts(1, other_nlocs);

  ioda::Selection feSelect;
  feSelect.extent({total_nlocs}).select({ioda::SelectionOperator::SET, other_start, counts});
  ioda::Selection beSelect;
  beSelect.select({ioda::SelectionOperator::SET, start, counts});

  std::vector<float> buffer(other_nlocs);
  other.og_.vars["longitude"].read<float>(gsl::make_span(buffer));
  og_.vars["longitude"].write<float>(buffer, feSelect, beSelect);

  other.og_.vars["latitude"].read<float>(gsl::make_span(buffer));
  og_.vars["latitude"].write<float>(buffer, feSelect, beSelect);

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
  return times_.size();
}

// -------------------------------------------------------------------------------------------------

std::vector<float> Locations::lons() const {
  ASSERT(og_.vars.exists("longitude"));
  const size_t nlocs = size();
  std::vector<float> lons(nlocs);
  og_.vars["longitude"].read<float>(gsl::make_span(lons));
  return lons;
}

// -------------------------------------------------------------------------------------------------

std::vector<float> Locations::lats() const {
  ASSERT(og_.vars.exists("latitude"));
  const size_t nlocs = size();
  std::vector<float> lats(nlocs);
  og_.vars["latitude"].read<float>(gsl::make_span(lats));
  return lats;
}

// -------------------------------------------------------------------------------------------------

void Locations::initializeObsGroup(const size_t nlocs) {
  ioda::Engines::BackendCreationParameters params;
  ioda::Group g = constructBackend(ioda::Engines::BackendNames::ObsStore, params);
  const ioda::NewDimensionScales_t dimScales{{
      ioda::NewDimensionScale<int>("nlocs", nlocs, ioda::Unlimited, nlocs)}};
  og_ = ioda::ObsGroup::generate(g, dimScales);
}

// -------------------------------------------------------------------------------------------------

void Locations::print(std::ostream & os) const {
  os << "Lat/lon/time locations: " << size() << " locations on this task " << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ufo
