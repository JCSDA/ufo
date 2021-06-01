/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ObsProfileAverage.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsProfileAverage> obsProfileAverageMaker_("ProfileAverage");
// -----------------------------------------------------------------------------

ObsProfileAverage::ObsProfileAverage(const ioda::ObsSpace & odb,
                                     const eckit::Configuration & config)
  : ObsOperatorBase(odb, config), odb_(odb)
{
  oops::Log::trace() << "ObsProfileAverage constructor starting" << std::endl;

  // Ensure observations have been grouped into profiles.
  if (odb_.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!odb_.has("MetaData", "extended_obs_space"))
    throw eckit::UserError("The extended obs space has not been produced", Here());

  // Check the observed pressure is present.
  if (!odb_.has("MetaData", "air_pressure"))
    throw eckit::UserError("air_pressure@MetaData not present", Here());

  // Set up configuration options.
  options_.validateAndDeserialize(config);

  // Add air_pressure_levels to the list of variables used in this operator.
  // todo(ctgh): obtain the required simulated variables in a future PR.
  requiredVars_ += oops::Variables({"air_pressure_levels"});

  // If required, set up vectors for OPS comparison.
  if (options_.compareWithOPS.value())
    setUpAuxiliaryReferenceVariables();

  oops::Log::trace() << "ObsProfileAverage constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsProfileAverage::~ObsProfileAverage() {
  oops::Log::trace() << "ObsProfileAverage destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverage::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                    ObsDiagnostics & ydiags) const {
  oops::Log::trace() << "ObsProfileAverage: simulateObs entered" << std::endl;

  const float missing = util::missingValue(missing);

  // Get variable indicating which section of the ObsSpace a profile is in.
  std::vector<int> extended_obs_space(odb_.nlocs());
  odb_.get_db("MetaData", "extended_obs_space", extended_obs_space);

  // Get observed pressure.
  std::vector<float> pressure_obs(odb_.nlocs());
  odb_.get_db("MetaData", "air_pressure", pressure_obs);

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<std::size_t> &recnums = odb_.recidx_all_recnums();
  // Number of profiles in the original ObsSpace.
  const std::size_t nprofs = recnums.size() / 2;

  // Loop over profiles.
  for (std::size_t jprof = 0; jprof < nprofs; ++jprof) {
    oops::Log::debug() << "Profile " << (jprof + 1) << " / " << nprofs << std::endl;

    // Get locations of profile in the original ObsSpace and
    // the corresponding profile in the extended ObsSpace.
    // Assuming the extended ObsSpace has been configured correctly, which is
    // checked above, the profile in the extended ObsSpace is always located
    // nprofs positions further on than the profile in the original ObsSpace.
    const std::vector<std::size_t> &locsOriginal = odb_.recidx_vector(recnums[jprof]);
    const std::vector<std::size_t> &locsExtended = odb_.recidx_vector(recnums[jprof + nprofs]);

    // Observed pressures for this profile.
    std::vector<float> pObs;
    for (const std::size_t loc : locsOriginal)
      pObs.push_back(pressure_obs[loc]);

    // Set up GeoVaLs and H(x) vectors.
    // Number of levels for air_pressure_levels.
    const std::size_t nlevs = gv.nlevs("air_pressure_levels");
    // Vector storing location for each level along the slant path.
    // Initially the first location in the profile is used everywhere.
    std::vector<std::size_t> slant_path_location(nlevs, 0);
    // Vector used to store different pressure GeoVaLs.
    std::vector <float> pressure_gv(nlevs);
    // Loop over model levels and find intersection of profile with model layer boundary.
    for (std::size_t mlev = 0; mlev < nlevs; ++mlev) {
      for (int iter = 0; iter < options_.numIntersectionIterations.value(); ++iter) {
        for (std::size_t jloc = slant_path_location[mlev]; jloc < pObs.size(); ++jloc) {
          gv.getAtLocation(pressure_gv, "air_pressure_levels", jloc);
          // If pressure has not been recorded, move to the next level.
          if (pObs[jloc] == missing) continue;
          // Break from the loop if the observed pressure is lower than
          // the pressure of this model level.
          if (pObs[jloc] <= pressure_gv[mlev]) break;
          // Record the value of this location at this level and all above.
          // This ensures that missing values are dealt with correctly.
          for (std::size_t mlevcolumn = mlev; mlevcolumn < nlevs; ++mlevcolumn)
            slant_path_location[mlevcolumn] = jloc;
        }
      }
    }

    // Fill slanted pressure vector.
    // todo(ctgh): fill other simulated variables in a future PR.
    std::vector<float> slant_pressure;
    for (std::size_t mlev = 0; mlev < nlevs; ++mlev) {
      const std::size_t jloc = slant_path_location[mlev];
      gv.getAtLocation(pressure_gv, "air_pressure_levels", jloc);
      slant_pressure.push_back(pressure_gv[mlev]);
    }

    // Fill H(x) in the extended ObsSpace.
    // todo(ctgh): do this in a future PR.
    // for (std::size_t jloc = 0; jloc < locsExtended.size(); ++jloc) {
    //  for (std::size_t jvar = 0; jvar < ovec.nvars(); ++jvar) {
    //    ovec[(jloc + locsExtended.front()) * ovec.nvars() + jvar] = slant_variable[jloc];
    //  }
    // }

    // If required, compare slant path locations and slant pressure with OPS.
    if (options_.compareWithOPS.value())
      compareAuxiliaryVariables(locsExtended,
                                slant_path_location,
                                slant_pressure);
  }
  oops::Log::trace() << "ObsProfileAverage: simulateObs exit " <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsProfileAverage::setUpAuxiliaryReferenceVariables() {
  if (!(odb_.has("MetOfficeHofX", "slant_path_location") &&
        odb_.has("MetOfficeHofX", "slant_pressure")))
    throw eckit::UserError("At least one reference variable is not present", Here());
  // Get reference values of the slant path locations and pressures.
  slant_path_location_ref_.resize(odb_.nlocs());
  slant_pressure_ref_.resize(odb_.nlocs());
  odb_.get_db("MetOfficeHofX", "slant_path_location", slant_path_location_ref_);
  odb_.get_db("MetOfficeHofX", "slant_pressure", slant_pressure_ref_);
}

// -----------------------------------------------------------------------------

void ObsProfileAverage::compareAuxiliaryVariables
(const std::vector<std::size_t> &locsExtended,
 const std::vector<std::size_t> &slant_path_location,
 const std::vector<float> &slant_pressure) const {
  std::vector<int> slant_path_location_ref_profile;
  std::vector<float> slant_pressure_ref_profile;
  for (const std::size_t loc : locsExtended) {
    slant_path_location_ref_profile.push_back(slant_path_location_ref_[loc]);
    slant_pressure_ref_profile.push_back(slant_pressure_ref_[loc]);
  }
  for (std::size_t jloccomp = 0; jloccomp < locsExtended.size(); ++jloccomp) {
    if (slant_path_location[jloccomp] != slant_path_location_ref_profile[jloccomp])
      throw eckit::BadValue("Mismatch for slant_path_location, jloccomp = " +
                            jloccomp, Here());
    if (!oops::is_close_relative(slant_pressure[jloccomp],
                                 slant_pressure_ref_profile[jloccomp],
                                 1e-9f))
      throw eckit::BadValue("Mismatch for slant_pressure, jloccomp = " +
                            jloccomp, Here());
  }
}

// -----------------------------------------------------------------------------

void ObsProfileAverage::print(std::ostream & os) const {
  os << "ObsProfileAverage operator" << std::endl;
  os << "config = " << options_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
