/* * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProfileAverageObsToModLevels.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/profile/ProfileAverageUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// ProfileAverageObsToModLevels: average the observation values on to model levels.

ProfileAverageObsToModLevels::ProfileAverageObsToModLevels(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ProfileAverageObsToModLevels constructor" << std::endl;

  // Ensure observations have been grouped into profiles.
  if (obsdb_.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!obsdb_.has("MetaData", "extendedObsSpace"))
    throw eckit::UserError("The extended obs space has not been produced", Here());

  // Get parameters from options
  allvars_ += parameters_.observation_vertical_coordinate;
  allvars_ += parameters_.model_vertical_coordinate;
}

// -----------------------------------------------------------------------------

ProfileAverageObsToModLevels::~ProfileAverageObsToModLevels() {
  oops::Log::trace() << "ProfileAverageObsToModLevels destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the Average Observations To Model Levels filter.

void ProfileAverageObsToModLevels::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ProfileAverageObsToModLevels postFilter" << std::endl;
  oops::Log::trace() << "ProfileAverageObsToModLevels obserr: " << *obserr_ << std::endl;

  const oops::ObsVariables observed = obsdb_.obsvariables();
  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");

  // Number of locations (including extended space)
  const size_t nlocs = obsdb_.nlocs();

  // Get model vertical coordinate (in the extended space)
  std::vector<float> model_vert_coord(nlocs);
  data_.get(parameters_.model_vertical_coordinate, model_vert_coord);

  // Get obs vertical coordinate
  std::vector<float> obs_vert_coord(nlocs);
  data_.get(parameters_.observation_vertical_coordinate, obs_vert_coord);

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> &unique_recnums = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;
  // Number of profiles in the original ObsSpace.
  const size_t nprofs = unique_recnums.size() / 2;
  Variables varhofx(filtervars_, "HofX");

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    const std::string varname = filtervars.variable(jv).variable();
    oops::Log::debug() << "varname: " << varname << std::endl;
    const size_t iv = observed.find(varname);
    // H(x)
    std::vector<float> hofx(nlocs);
    data_.get(varhofx.variable(jv), hofx);

    // Loop over profiles.
    for (size_t jprof = 0; jprof < nprofs; ++jprof) {
      oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;
      // Get locations of profile in the original ObsSpace and
      // the corresponding profile in the extended ObsSpace.
      // Assuming the extended ObsSpace has been configured correctly, which is
      // checked above, the profile in the extended ObsSpace is always located
      // nprofs positions further on than the profile in the original ObsSpace.
      const std::vector<size_t> &locsOriginal = obsdb_.recidx_vector(unique_recnums[jprof]);
      const std::vector<size_t> &locsExt = obsdb_.recidx_vector(unique_recnums[jprof + nprofs]);
      oops::Log::debug() << "locsOriginal: " << locsOriginal << std::endl;
      oops::Log::debug() << "locsExt: " << locsExt << std::endl;

      // Create a vector mapping each obs index to a model index (many-to-1):
      const std::vector<size_t> link_obs_and_mod_inds = ufo::linkObsAndModLevIndices
                                                            (locsOriginal,
                                                             locsExt,
                                                             obs_vert_coord,
                                                             model_vert_coord,
                                                             nlocs,
                                                             false);  // link mod->obs AND obs->mod
      // Given a flag on some obs level, set that flag at the corresponding model level:
      ufo::setModLevelFlags(locsExt,
                            link_obs_and_mod_inds,
                            flagged[jv],
                            (*flags_)[iv]);
      // Compute increments at each obs level, take their average for each model level,
      //  and add the average to the HofX at that model level:
      ufo::averageObsToModLevels(locsOriginal,
                                 locsExt,
                                 link_obs_and_mod_inds,
                                 apply,
                                 (*flags_)[iv],
                                 hofx,
                                 obs[jv]);
    }  // profile jprof
    // Save new obs to obsSpace:
    obsdb_.put_db("DerivedObsValue", varname, obs[jv]);
  }  // filter variable jv
}

// -----------------------------------------------------------------------------

void ProfileAverageObsToModLevels::print(std::ostream & os) const {
  os << "ProfileAverageObsToModLevels: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
