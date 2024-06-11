/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ProfileAverageObsPressure.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/profile/SlantPathLocations.h"

namespace ufo {

static ObsFunctionMaker<ProfileAverageObsPressure>
                       makerProfileAverageObsPressure_("ProfileAverageObsPressure");

// -----------------------------------------------------------------------------

ProfileAverageObsPressure::ProfileAverageObsPressure(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  // Observed air pressure
  invars_ += Variable(options_.observation_vertical_coordinate.value());

  // Air pressure GeoVaLs
  invars_ += Variable(std::string("GeoVaLs/") + options_.model_vertical_coordinate.value());
}

// -----------------------------------------------------------------------------

ProfileAverageObsPressure::~ProfileAverageObsPressure() {}

// -----------------------------------------------------------------------------

void ProfileAverageObsPressure::compute(const ObsFilterData & in,
                                        ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ProfileAverageObsPressure::compute started" << std::endl;

  // ObsSpace.
  const ioda::ObsSpace & obsdb = in.obsspace();

  // Ensure observations have been grouped into profiles.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!obsdb.has("MetaData", "extendedObsSpace"))
    throw eckit::UserError("The extended obs space has not been produced", Here());

  // GeoVaLs.
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Number of locations.
  const size_t nlocs = obsdb.nlocs();

  // Correspondence between record numbers and indices in the data sample.
  const std::vector<std::size_t> &recnums = obsdb.recidx_all_recnums();

  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Get observed vertical coordinate.
  std::vector <float> vert_coord_obs(nlocs);
  in.get(Variable(options_.observation_vertical_coordinate.value()), vert_coord_obs);

  // Name of model vertical coordinate.
  const oops::Variable model_vertical_coordinate{
    options_.model_vertical_coordinate.value()};

  // Vector holding model vertical coordinate at different locations.
  std::vector<double> var_gv(gv->nlevs(model_vertical_coordinate));

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
    oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

    // Get locations of profile in the original ObsSpace and
    // the corresponding profile in the extended ObsSpace.
    // Assuming the extended ObsSpace has been configured correctly, which is
    // checked above, the profile in the extended ObsSpace is always located
    // nprofs positions further on than the profile in the original ObsSpace.
    const std::vector<std::size_t> &locsOriginal = obsdb.recidx_vector(recnums[jprof]);
    const std::vector<std::size_t> &locsExtended = obsdb.recidx_vector(recnums[jprof + nprofs]);

    // Retrieve slant path locations.
    // This function ensures the GeoVaLs are in the correct order.
    const std::vector<std::size_t> slant_path_location =
      ufo::getSlantPathLocations(obsdb,
                                 *gv,
                                 locsOriginal,
                                 options_.observation_vertical_coordinate,
                                 oops::Variable{options_.model_vertical_coordinate},
                                 options_.numIntersectionIterations.value() - 1);

    // Write out values to output vector
    // Pressures in original profile.
    for (size_t loc : locsOriginal)
      out[0][loc] = vert_coord_obs[loc];

    // Pressures in averaged profile.
    for (size_t idx = 0; idx < locsExtended.size() && idx < var_gv.size(); ++idx) {
      // GeoVaL vector for this variable at the slant location.
      gv->getAtLocation(var_gv, model_vertical_coordinate, slant_path_location[idx]);
      // Transfer value into averaged profile.
      out[0][locsExtended[idx]] = var_gv[var_gv.size() - 1 - idx];
    }
  }

  oops::Log::trace() << "ProfileAverageObsPressure::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ProfileAverageObsPressure::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
