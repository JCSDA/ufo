/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <set>
#include <utility>

#include "ufo/filters/QCflags.h"

#include "ufo/profile/ProfileAverageUtils.h"

namespace ufo {

  void ProfileAverageUtils::passNonMissingAveragedObservations
  (ProfileDataHandler & profileDataHandler,
   std::vector <ProfileDataHolder> & profiles,
   const std::string & flag_name,
   const std::string & obs_name)
  {
    const float missing = util::missingValue<float>();
    const size_t halfnprofs = profileDataHandler.getObsdb().nrecs() / 2;

    // Retrieve filter QC flags.
    ioda::ObsDataVector<int> & filterFlagsAll = profileDataHandler.getFilterFlags();
    ioda::ObsDataRow<int> & filterFlags = filterFlagsAll[flag_name];

    // The averaging routines fill in previously missing values, so the associated flag is updated.
    const std::vector<size_t> &recnums = profileDataHandler.getObsdb().recidx_all_recnums();
    for (std::size_t jprof = halfnprofs; jprof < 2 * halfnprofs; ++jprof) {
      const std::vector<std::size_t> &locsExtended =
        profileDataHandler.getObsdb().recidx_vector(recnums[jprof]);
      auto& profileExtended = profiles[jprof];
      const std::vector <float> &obs = profileExtended.get<float>(obs_name);
      for (size_t idx = 0; idx < locsExtended.size(); ++idx) {
        if (obs[idx] != missing)
          filterFlags[locsExtended[idx]] = QCflags::pass;
      }
    }
  }

  void ProfileAverageUtils::fillValidationData
  (ProfileDataHolder & profile,
   bool extended_obs_space,
   const std::string & average_name,
   const std::string & qcflags_name,
   const oops::Variable & geovals_testreference_name,
   const oops::Variable & geovals_qcflags_name)
  {
    const float missing = util::missingValue<float>();

    // Retrieve, then save, the OPS versions of the variable averaged onto model levels,
    // and the QC flags associated with the averaging process.
    // The quantities retrieved depend on whether the profile lies in the original
    // or averaged sections of the ObsSpace.
    if (extended_obs_space) {
      std::vector<float>& averaged_values =
        profile.getGeoVaLVector(geovals_testreference_name);
      // Ensure all vectors are the correct size to be saved to the ObsSpace.
      const size_t numModelLevels = profile.getNumProfileLevels();
      averaged_values.resize(numModelLevels, missing);
      profile.set<float>(ufo::addOPSPrefix(average_name), std::move(averaged_values));
      // The QC flags are stored as floats but are converted to integers here.
      // Due to the loss of precision, 5 must be added to the missing value.
      const std::vector <float>& average_qcflags_float =
        profile.getGeoVaLVector(geovals_qcflags_name);
      std::vector <int> average_qcflags_int
        (average_qcflags_float.begin(),
         average_qcflags_float.end());
      std::replace(average_qcflags_int.begin(),
                   average_qcflags_int.end(),
                   -2147483648L,
                   -2147483643L);
      // Ensure all vectors are the correct size to be saved to the ObsSpace.
      average_qcflags_int.resize(numModelLevels, 0);
      profile.set<int>(ufo::addOPSPrefix(qcflags_name), std::move(average_qcflags_int));
    } else {
      // Create a copy here because the vector will be used later in the routine.
      std::vector <float> avg = profile.get<float>(average_name);
      profile.set<float>(ufo::addOPSPrefix(average_name), std::move(avg));
    }
  }

  // -----------------------------------------------------------------------------

  // Create vector linking observation levels and model levels.
  std::vector<size_t> linkObsAndModLevIndices(const std::vector<size_t> &locsOriginal,
                                        const std::vector<size_t> &locsExt,
                                        const std::vector<float> &obs_vert_coord,
                                        const std::vector<float> &model_vert_coord,
                                        const size_t &nlocs,
                                        const bool &obs_to_mod_only) {
    const int64_t missingValInt = util::missingValue<int64_t>();  // size_t type
    const float missingValFloat = util::missingValue<float>();
    const size_t nlocs_obs = locsOriginal.size();
    const size_t nlocs_ext = locsExt.size();
    const bool vert_coord_increasing = validateVertCoords(locsOriginal,
                                                          locsExt,
                                                          obs_vert_coord,
                                                          model_vert_coord);

    std::vector<size_t> vector_linking_indices(nlocs, missingValInt);
    // Find mid-levels (halfway between each model level)
    std::vector<float> mid_levels(nlocs_ext);
    for (size_t mlev = 0; mlev < nlocs_ext-1; ++mlev) {
      mid_levels[mlev] = 0.5*(model_vert_coord[locsExt[mlev]]+
                              model_vert_coord[locsExt[mlev+1]]);
      oops::Log::debug() << "mid_levels[mlev]: " << mid_levels[mlev] << std::endl;
      oops::Log::debug() << "model_vert_coord[mlev]: " <<
                            model_vert_coord[locsExt[mlev]] << std::endl;
    }
    if (nlocs_ext > 1) {
      // set final (bottom-most) mid-level:
      mid_levels[nlocs_ext-1] = mid_levels[nlocs_ext-2] +
                              2*(model_vert_coord[locsExt[nlocs_ext-1]]-mid_levels[nlocs_ext-2]);
    } else {  // single model level
      mid_levels[0] = 2*model_vert_coord[locsExt[0]];
    }
    // Fill the vector indices mapping obs levels to model levels:
    for (size_t jlev = 0; jlev < nlocs_obs; ++jlev) {
      for (size_t mlev = 0; mlev < nlocs_ext; ++mlev) {
        if ( (obs_vert_coord[locsOriginal[jlev]] != missingValFloat) &&
              ( (vert_coord_increasing &&
                mid_levels[mlev] > obs_vert_coord[locsOriginal[jlev]]) ||
              (!vert_coord_increasing &&
                mid_levels[mlev] < obs_vert_coord[locsOriginal[jlev]]) ) ) {
          vector_linking_indices[locsOriginal[jlev]] = locsExt[mlev];
          break;
        }
      }  // loop over model levels
    }  // loop over obs levels
    if (!obs_to_mod_only) {  // link model levels to obs levels as well
      for (size_t mlev = 0; mlev < nlocs_ext; ++mlev) {
        for (size_t jlev = 0; jlev < nlocs_obs; ++jlev) {
          if ( (vert_coord_increasing &&
                obs_vert_coord[locsOriginal[jlev]] >= model_vert_coord[locsExt[mlev]]) ||
              (!vert_coord_increasing &&
                obs_vert_coord[locsOriginal[jlev]] <= model_vert_coord[locsExt[mlev]]) ) {
            vector_linking_indices[locsExt[mlev]] = locsOriginal[jlev];
            break;
          }
        }  // loop over obs levels
      }  // loop over model levels
    }
    return vector_linking_indices;
  }

  // -----------------------------------------------------------------------------

  // Copy values of `flagged` and `flags`, from model levels to their corresponding obs
  //   levels.
  void setModLevelFlags(const std::vector<size_t> &locsExt,
                        const std::vector<size_t> &vector_linking_indices,
                        std::vector<bool> &flagged,
                        std::vector<int> &flags) {
    const size_t nlocs_ext = locsExt.size();
    const int64_t missingValInt = util::missingValue<int64_t>();  // size_t type

    // Set flags
    for (size_t mlev = 0; mlev < nlocs_ext; ++mlev) {
      const size_t obsInd = vector_linking_indices[locsExt[mlev]];
      if (obsInd != missingValInt) {
        flags[locsExt[mlev]] = flags[obsInd];
        flagged[locsExt[mlev]] = flagged[obsInd];
      }
    }
  }

  // -----------------------------------------------------------------------------

  // Compute average of observation values associated with each model level, and write
  //   them into the extended space.
  void averageObsToModLevels(const std::vector<size_t> &locsOriginal,
                            const std::vector<size_t> &locsExt,
                            const std::vector<size_t> &vector_linking_indices,
                            const std::vector<bool> &apply,
                            const std::vector<int> &flags,
                            const std::vector<float> &hofx,
                            std::vector<float> &obs) {
    const int64_t missingValInt = util::missingValue<int64_t>();  // size_t type
    const float missingValFloat = util::missingValue<float>();
    const size_t nlocs_obs = locsOriginal.size();
    const size_t nlocs_ext = locsExt.size();
    const size_t nlocs = hofx.size();

    // Compute increments (obs minus BG)
    std::vector<float> increments(nlocs, missingValFloat);
    for (size_t jlev : locsOriginal) {
      if (apply[jlev] && flags[jlev] == QCflags::pass &&
          obs[jlev] != missingValFloat &&
          hofx[jlev] != missingValFloat) {
        increments[jlev] = obs[jlev] - hofx[jlev];
      }
    }

    // Average increments and add to BG
    for (size_t mlev = 0; mlev < nlocs_ext; ++mlev) {
      float sum_increments = 0.0;
      size_t count_increments = 0;
      for (size_t jlev : locsOriginal) {
        if (vector_linking_indices[jlev] == locsExt[mlev] &&
            increments[jlev] != missingValFloat) {
          sum_increments += increments[jlev];
          count_increments++;
          oops::Log::debug() << "vector_linking_indices[jlev]: " << vector_linking_indices[jlev] <<
                                "[" << jlev << "]" << std::endl;
        }
      }
      oops::Log::debug() << "Increments[mlev]: " << sum_increments <<
                            "[" << mlev << "]" << std::endl;
      if (count_increments > 0) {
        obs[locsExt[mlev]] = hofx[locsExt[mlev]] + sum_increments / count_increments;
      }
    }
  }

  // -----------------------------------------------------------------------------

  // Check model vertical coordinate non-zero and goes in same direction as observation
  //  vertical coordinate. Returns true if the vertical coordinates are increasing.
  bool validateVertCoords(const std::vector<size_t> &locsOriginal,
                          const std::vector<size_t> &locsExt,
                          const std::vector<float> &obs_vert_coord,
                          const std::vector<float> &model_vert_coord) {
    const size_t nlocs_obs = locsOriginal.size();
    const size_t nlocs_ext = locsExt.size();

    // Stop if no model vertical coordinate:
    std::vector<size_t> model_vert_coord_ext(nlocs_ext);
    for (size_t mlev = 0; mlev < nlocs_ext; ++mlev) {
      model_vert_coord_ext[mlev] = model_vert_coord[locsExt[mlev]];
      oops::Log::debug() << "model_vert_coord[" << mlev << "]: " <<
                            model_vert_coord[locsExt[mlev]] << std::endl;
    }
    if (std::all_of(model_vert_coord_ext.begin(), model_vert_coord_ext.end(),
                    [](float i) { return (i == 0); })) {
      throw eckit::UserError(": The model vertical coordinate extended space is all zeros. "
                              "Did you remember to include the vertical coordinate variable "
                              "when applying the ProfileAverage obsOperator?", Here());
    }

    // Stop if model vertical coordinate not in same direction as obs vertical coordinate:

    // if no clear case is found, all levels are equal and default to increasing model levels
    bool mod_vert_coord_increasing = true;
    bool mod_vert_coord_equal = true;
    // check every model level until a clear non-equal-levels case is found
    for (size_t mLev = 0; mLev < nlocs_ext-1; ++mLev) {
      if (model_vert_coord[locsExt[mLev+1]] > model_vert_coord[locsExt[mLev]]) {
        mod_vert_coord_increasing = true;
        mod_vert_coord_equal = false;
        break;
      } else if (model_vert_coord[locsExt[mLev+1]] < model_vert_coord[locsExt[mLev]]) {
        mod_vert_coord_increasing = false;
        mod_vert_coord_equal = false;
        break;
      }
    }

    if (nlocs_obs <= 1) {
      return mod_vert_coord_increasing;
    }

    // if no clear case is found, all levels are equal and default to increasing obs levels
    bool obs_vert_coord_increasing = true;
    bool obs_vert_coord_equal = true;
    // check every obs level until a clear non-equal-levels case is found
    for (size_t oLev = 0; oLev < nlocs_obs-1; ++oLev) {
      if (obs_vert_coord[locsOriginal[oLev+1]] > obs_vert_coord[locsOriginal[oLev]]) {
        obs_vert_coord_increasing = true;
        obs_vert_coord_equal = false;
        break;
      } else if (obs_vert_coord[locsOriginal[oLev+1]] < obs_vert_coord[locsOriginal[oLev]]) {
        obs_vert_coord_increasing = false;
        obs_vert_coord_equal = false;
        break;
      }
    }

    if (obs_vert_coord_equal) {
      oops::Log::debug() << " obs_vert_coord_equal " << obs_vert_coord[locsOriginal[1]] << " == "
                         << obs_vert_coord[locsOriginal[0]] << std::endl;
      return mod_vert_coord_increasing;
    } else if (mod_vert_coord_equal) {
      oops::Log::debug() << " mod_vert_coord_equal " << model_vert_coord[locsExt[1]] << " == "
                         << model_vert_coord[locsExt[0]] << std::endl;
      return obs_vert_coord_increasing;
    }
    if (obs_vert_coord_increasing && !mod_vert_coord_increasing) {
      throw eckit::UserError(": The model vertical coordinate is decreasing, but the observation "
            "vertical coordinate is increasing. They must go in the same direction.", Here());
    } else if (!obs_vert_coord_increasing && mod_vert_coord_increasing) {
      throw eckit::UserError(": The model vertical coordinate is increasing, but the observation "
            "vertical coordinate is decreasing. They must go in the same direction.", Here());
    } else {
      return obs_vert_coord_increasing;
    }
  }

  std::string addOPSPrefix(const std::string & fullname)
  {
    std::string varname;
    std::string groupname;
    ufo::splitVarGroup(fullname, varname, groupname);
    return groupname + std::string("/OPS_") + varname;
  }
}  // namespace ufo
