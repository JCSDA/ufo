/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/FillAveragedProfileData.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/GeoVaLs.h"
#include "ufo/profile/SlantPathLocations.h"

namespace ufo {

static ObsFunctionMaker<FillAveragedProfileData<float>> floatMaker("FillAveragedProfileData");
static ObsFunctionMaker<FillAveragedProfileData<int>> intMaker("FillAveragedProfileData");
static ObsFunctionMaker<FillAveragedProfileData<std::string>>
stringMaker("FillAveragedProfileData");
static ObsFunctionMaker<FillAveragedProfileData<util::DateTime>>
dateTimeMaker("FillAveragedProfileData");

// -----------------------------------------------------------------------------

template <typename FunctionValue>
FillAveragedProfileData<FunctionValue>::FillAveragedProfileData
(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  // Variable to copy from the original to the averaged profile.
  invars_ += Variable(options_.variable_to_copy);

  // Air pressure GeoVaLs.
  invars_ += Variable(std::string("GeoVaLs/") + options_.model_vertical_coordinate.value());
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
void FillAveragedProfileData<FunctionValue>::compute
(const ObsFilterData & in,
 ioda::ObsDataVector<FunctionValue> & out) const {
  oops::Log::trace() << "FillAveragedProfileData::compute started" << std::endl;

  // Fill values of \p variable_to_copy in the averaged profiles.
  fillAverageProfile(in, out);

  oops::Log::trace() << "FillAveragedProfileData::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
const ufo::Variables & FillAveragedProfileData<FunctionValue>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
void FillAveragedProfileData<FunctionValue>::fillAverageProfile
(const ObsFilterData & in,
 ioda::ObsDataVector<FunctionValue> & out) const {
  const FunctionValue missing = util::missingValue<FunctionValue>();

  // ObsSpace.
  ioda::ObsSpace & obsdb = in.obsspace();

  // Ensure observations have been grouped into profiles.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!obsdb.has("MetaData", "extendedObsSpace"))
    throw eckit::UserError("The extended obs space has not been produced", Here());

  // Model vertical coordinate name.
  const std::string model_vertical_coordinate =
    options_.model_vertical_coordinate.value();

  // This coordinate must be a pressure for the slant path location algorithm to work correctly.
  if (model_vertical_coordinate.find("pressure") == model_vertical_coordinate.npos)
    throw eckit::UserError("Model vertical coordinate must be a pressure", Here());

  // GeoVaLs.
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Number of locations.
  const size_t nlocs = obsdb.nlocs();

  // Correspondence between record numbers and indices in the data sample.
  const std::vector<std::size_t> &recnums = obsdb.recidx_all_recnums();

  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Get variable whose values will be copied to the averaged profile.
  std::vector <FunctionValue> variable_to_copy(nlocs);
  in.get(Variable(options_.variable_to_copy.value()), variable_to_copy);

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
    const std::vector<std::size_t> slant_path_location =
      ufo::getSlantPathLocations(obsdb,
                                 *gv,
                                 locsOriginal,
                                 options_.observation_vertical_coordinate,
                                 oops::Variable{options_.model_vertical_coordinate},
                                 options_.numIntersectionIterations.value() - 1);

    // Output values in original profile.
    for (size_t loc : locsOriginal)
      out[0][loc] = variable_to_copy[loc];

    // Output values in averaged profile.
    for (size_t mlev = 0; mlev < locsExtended.size(); ++mlev)
      out[0][locsExtended[mlev]] = variable_to_copy[slant_path_location[mlev]];
  }
}

}  // namespace ufo
