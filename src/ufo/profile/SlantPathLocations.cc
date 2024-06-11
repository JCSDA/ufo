/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "ioda/ObsSpace.h"

#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/profile/SlantPathLocations.h"
#include "ufo/utils/StringUtils.h"  // for splitVarGroup

namespace ufo {

  std::vector<std::size_t> getSlantPathLocations(const ioda::ObsSpace & odb,
                                                 const GeoVaLs & gv,
                                                 const std::vector<std::size_t> & locs,
                                                 const std::string & obsVerticalCoord,
                                                 const oops::Variable & modelVerticalCoord,
                                                 const int itermax) {
    const float missing = util::missingValue<float>();

    // Get observed pressure.
    std::vector<float> pressure_obs(odb.nlocs());
    std::string obsVar, obsGroup;
    splitVarGroup(obsVerticalCoord, obsVar, obsGroup);
    odb.get_db(obsGroup, obsVar, pressure_obs);
    // Number of levels for model pressure.
    const std::size_t nlevs_p = gv.nlevs(modelVerticalCoord);
    // Vector storing location for each level along the slant path.
    // Initially the first location in the profile is used everywhere.
    std::vector<std::size_t> slant_path_location(nlevs_p, locs.front());
    // Vector used to store different pressure GeoVaLs.
    std::vector <float> pressure_gv(nlevs_p);

    // Loop over model levels and find intersection of profile with model layer boundary.
    // This can be performed multiple times in order to account for slanted model levels.
    std::size_t idxstart = 0;  // Starting index for loop over levels.
    // This counter records locations in the slanted profile.
    // It is initialised at the first location in the original profile.
    std::size_t jlocslant = locs.front();
    // Loop over each model level in turn.
    for (int mlev = nlevs_p - 1; mlev >= 0; --mlev) {
      for (int iter = 0; iter <= itermax; ++iter) {
        // Get the GeoVaL that corresponds to the current slanted profile location.
        gv.getAtLocation(pressure_gv, modelVerticalCoord, jlocslant);
        // The GeoVaLs must be ordered from top to bottom for this algorithm to work.
        if (pressure_gv.front() > pressure_gv.back())
          throw eckit::BadValue("Pressure GeoVaLs are in the wrong order.", Here());
        // Define an iteration-specific location that is initialised to the
        // current slanted profile location.
        std::size_t jlociter = jlocslant;
        // Determine the intersection of the observed profile with the current model level.
        // The intersection is taken to be the location with the largest observed pressure
        // that is less than or equal to the model pressure at this level.
        for (std::size_t idx = idxstart; idx < locs.size(); ++idx) {
          // Intersection location.
          const std::size_t jlocintersect = locs[idx];
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
          for (int mlevcolumn = nlevs_p - 1 - mlev; mlevcolumn < nlevs_p; ++mlevcolumn)
           slant_path_location[mlevcolumn] = jlocslant;
        }
      }
    }

    return slant_path_location;
  }
}  // namespace ufo
