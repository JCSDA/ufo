/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FloatCompare.h"
#include "oops/util/missingValues.h"

#include "ufo/profile/ObsProfileAverageData.h"

namespace ufo {

  ObsProfileAverageData::ObsProfileAverageData(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
    : odb_(odb)
  {
    // Ensure observations have been grouped into profiles.
    if (odb_.obs_group_vars().empty())
      throw eckit::UserError("Group variables configuration is empty", Here());

    // Ensure observations have been sorted by air pressure in descending order.
    if (odb_.obs_sort_var() != "air_pressure")
      throw eckit::UserError("Sort variable must be air_pressure", Here());
    if (odb_.obs_sort_order() != "descending")
      throw eckit::UserError("Profiles must be sorted in descending order", Here());

    // Check the ObsSpace has been extended. If this is not the case
    // then it will not be possible to access profiles in the original and
    // extended sections of the ObsSpace.
    if (!odb_.has("MetaData", "extended_obs_space"))
      throw eckit::UserError("The extended obs space has not been produced", Here());

    // Check the observed pressure is present. This is necessary in order to
    // determine the slant path locations.
    if (!odb_.has("MetaData", "air_pressure"))
      throw eckit::UserError("air_pressure@MetaData not present", Here());

    // Add air_pressure_levels to the list of variables used in this operator.
    // This GeoVaL is used to determine the slant path locations.
    requiredVars_ += oops::Variables({"air_pressure_levels"});

    // Add any simulated variables to the list of variables used in this operator.
    requiredVars_ += odb_.obsvariables();

    // Set up configuration options.
    options_.validateAndDeserialize(config);

    // If required, set up vectors for OPS comparison.
    if (options_.compareWithOPS.value())
      this->setUpAuxiliaryReferenceVariables();
  }

  const oops::Variables & ObsProfileAverageData::requiredVars() const {
    return requiredVars_;
  }

  std::vector<std::size_t> ObsProfileAverageData::getSlantPathLocations
  (const std::vector<std::size_t> & locsOriginal,
   const std::vector<std::size_t> & locsExtended,
   const GeoVaLs & gv) const
  {
    const float missing = util::missingValue(missing);

    // Get observed pressure.
    std::vector<float> pressure_obs(odb_.nlocs());
    odb_.get_db("MetaData", "air_pressure", pressure_obs);

    // Set up GeoVaLs and H(x) vectors.
    // Number of levels for air_pressure_levels.
    const std::size_t nlevs_p = gv.nlevs("air_pressure_levels");
    // Vector storing location for each level along the slant path.
    // Initially the first location in the profile is used everywhere.
    std::vector<std::size_t> slant_path_location(nlevs_p, locsOriginal.front());
    // Vector of slanted pressures, used for comparisons with OPS.
    std::vector<float> slant_pressure;
    // Vector used to store different pressure GeoVaLs.
    std::vector <float> pressure_gv(nlevs_p);

    // Loop over model levels and find intersection of profile with model layer boundary.
    // This can be performed multiple times in order to account for slanted model levels.
    std::size_t idxstart = 0;  // Starting index for loop over levels.
    // This counter records locations in the slanted profile.
    // It is initialised at the first location in the original profile.
    std::size_t jlocslant = locsOriginal.front();
    // Maximum iteration value.
    const int itermax = options_.numIntersectionIterations.value() - 1;
    // Loop over each model level in turn.
    for (std::size_t mlev = 0; mlev < nlevs_p; ++mlev) {
      for (int iter = 0; iter <= itermax; ++iter) {
        // Get the GeoVaL that corresponds to the current slanted profile location.
        gv.getAtLocation(pressure_gv, "air_pressure_levels", jlocslant);
        // Define an iteration-specific location that is initialised to the
        // current slanted profile location.
        std::size_t jlociter = jlocslant;
        // Determine the intersection of the observed profile with the current model level.
        // The intersection is taken to be the location with the largest observed pressure
        // that is less than or equal to the model pressure at this level.
        for (std::size_t idx = idxstart; idx < locsOriginal.size(); ++idx) {
          // Intersection location.
          const std::size_t jlocintersect = locsOriginal[idx];
          // If pressure has not been recorded, move to the next level.
          if (pressure_obs[jlocintersect] == missing) continue;
          // Break from the loop if the observed pressure is lower than
          // the pressure of this model level.
          if (pressure_obs[jlocintersect] <= pressure_gv[mlev]) break;
          // Update the iteration-specific location with the new intersection location.
          jlociter = jlocintersect;
          // Update the loop starting index when the last iteration is reached.
          if (iter == itermax) idxstart = idx;
        }
        // Modify the slanted location in the original profile.
        jlocslant = jlociter;
        if (iter == itermax) {
          // Record the value of the slant path location at this model level and all above.
          // This ensures that missing values are dealt with correctly.
          for (std::size_t mlevcolumn = mlev; mlevcolumn < nlevs_p; ++mlevcolumn)
           slant_path_location[mlevcolumn] = jlocslant;
          // Fill the slant pressure if the comparison with OPS will be performed.
          if (options_.compareWithOPS.value())
            slant_pressure.push_back(pressure_gv[mlev]);
        }
      }
    }

    // If required, compare slant path locations and slant pressure with OPS output.
    if (options_.compareWithOPS.value())
      this->compareAuxiliaryReferenceVariables(locsExtended,
                                               slant_path_location,
                                               slant_pressure);

    return slant_path_location;
  }

  void ObsProfileAverageData::setUpAuxiliaryReferenceVariables() {
    if (!(odb_.has("MetOfficeHofX", "slant_path_location") &&
          odb_.has("MetOfficeHofX", "slant_pressure")))
      throw eckit::UserError("At least one reference variable is not present", Here());
    // Get reference values of the slant path locations and pressures.
    slant_path_location_ref_.resize(odb_.nlocs());
    slant_pressure_ref_.resize(odb_.nlocs());
    odb_.get_db("MetOfficeHofX", "slant_path_location", slant_path_location_ref_);
    odb_.get_db("MetOfficeHofX", "slant_pressure", slant_pressure_ref_);
  }

  void ObsProfileAverageData::compareAuxiliaryReferenceVariables
  (const std::vector<std::size_t> & locsExtended,
   const std::vector<std::size_t> & slant_path_location,
   const std::vector<float> & slant_pressure) const {
    std::vector<int> slant_path_location_ref_profile;
    std::vector<float> slant_pressure_ref_profile;
    for (const std::size_t loc : locsExtended) {
      slant_path_location_ref_profile.push_back(slant_path_location_ref_[loc]);
      slant_pressure_ref_profile.push_back(slant_pressure_ref_[loc]);
    }
    std::stringstream errmsg;
    for (std::size_t jloccomp = 0; jloccomp < locsExtended.size(); ++jloccomp) {
      if (slant_path_location[jloccomp] != slant_path_location_ref_profile[jloccomp]) {
        errmsg << "Mismatch for slant_path_location, location = " << jloccomp
               << " (this code, OPS): "
               << slant_path_location[jloccomp] << ", "
               << slant_path_location_ref_profile[jloccomp];
        throw eckit::BadValue(errmsg.str(), Here());
      }
      if (!oops::is_close_relative(slant_pressure[jloccomp],
                                   slant_pressure_ref_profile[jloccomp],
                                   1e-5f)) {
        errmsg << "Mismatch for slant_pressure, location = " << jloccomp
               << " (this code, OPS): "
               << slant_pressure[jloccomp] << ", "
               << slant_pressure_ref_profile[jloccomp];
        throw eckit::BadValue(errmsg.str(), Here());
      }
    }
  }

  void ObsProfileAverageData::print(std::ostream & os) const {
    os << "ObsProfileAverage operator" << std::endl;
    os << "config = " << options_ << std::endl;
  }
}  // namespace ufo
