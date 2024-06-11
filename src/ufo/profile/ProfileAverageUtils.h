/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGEUTILS_H_
#define UFO_PROFILE_PROFILEAVERAGEUTILS_H_

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oops/util/missingValues.h"

#include "ufo/profile/ProfileDataHolder.h"

namespace ufo {
  class ProfileDataHandler;
}

namespace ufo {

  class ProfileAverageUtils {
   public:
    /// Modify filter QC flags based on values of the averaged data.
    static void passNonMissingAveragedObservations
      (ProfileDataHandler & profileDataHandler,
       std::vector <ProfileDataHolder> & profiles,
       const std::string & flag_name,
       const std::string & obs_name);

    /// Fill validation data for use in comparisons with OPS.
    static void fillValidationData
      (ProfileDataHolder & profile,
       bool extended_obs_space,
       const std::string & average_name,
       const std::string & qcflags_name,
       const oops::Variable & geovals_testreference_name,
       const oops::Variable & geovals_qcflags_name);

    /// Set values in a profile to missing.
    template <typename T>
    static void setProfileMissing
      (ProfileDataHolder & profile,
       const std::vector <std::string> & variableNames)
    {
      const T missing = util::missingValue<T>();
      for (const std::string & variableName : variableNames) {
        profile.set<T>(variableName,
                       std::move(std::vector<T>(profile.getNumProfileLevels(), missing)));
      }
    }

    /// Transfer one variable in a profile to another variable in the same profile.
    template <typename T>
      static void copyProfileValues(ProfileDataHolder & profile,
                                    const std::string & sourceVariableName,
                                    const std::string & destinationVariableName)
    {
      profile.set<T>(destinationVariableName,
                     std::move(profile.get<T>(sourceVariableName)));
    }
  };

  /// \brief Create vector mapping observation levels to model levels.
  /// \details Returns a vector `nlocs` long (original+extended), such that each element links
  ///  an original-space index (observation level) to an extended-space index (model level)
  ///  or vice versa, i.e.
  ///   `vector_linking_indices[locsOriginal[ind]]` is the extended-space index that corresponds
  ///   to `locsOriginal[ind]`;
  ///   `vector_linking_indices[locsExt[ind]]` is the original-space index that corresponds to
  ///   `locsExt[ind]`.
  ///
  ///  Many obs levels can link to the same model level: by finding the mid-levels halfway
  ///  between adjacent model levels, all the obs levels that fall between the two mid-levels
  ///  bounding a model level, are linked to that model level. (Some model levels may have 0
  ///  obs levels linking to them, but every non-missing obs level links to a model level,
  ///  if within the full vertical range of model levels.)
  ///
  ///  While the mapping from obs levels to model levels is many-to-one, the reverse mapping
  ///  from model levels to obs levels is one-to-one. Each model level has one corresponding
  ///  obs level, that being the one closest to it that equals or exceeds it in height (or
  ///  depth, in the oceans).
  ///
  /// \param[in] locsOriginal: indices of original space profile
  /// \param[in] locsExt: indices of the same profile in extended space
  /// \param[in] obs_vert_coord: observation vertical coordinate values
  /// \param[in] model_vert_coord: model vertical coordinate values (only the extended space
  ///            is used, so checks that they aren't all 0)
  /// \param[in] nlocs: number of original plus extended space levels in the profile
  /// \param[in] obs_to_mod_only: if false, fill both original- and extended-space portions
  ///            of the output `vector_linking_indices`; if true, only fill the original-
  ///            space portion (and save some computation)
  std::vector<size_t> linkObsAndModLevIndices(const std::vector<size_t> &locsOriginal,
                                              const std::vector<size_t> &locsExt,
                                              const std::vector<float> &obs_vert_coord,
                                              const std::vector<float> &model_vert_coord,
                                              const size_t &nlocs,
                                              const bool &obs_to_mod_only);

  /// \brief Copy values of `flagged` and `flags`, from model levels to their corresponding obs
  ///   levels.
  /// \details Using the result of `ufo::linkObsAndModLevelIndices`, set each model level's
  ///  `flagged` and `flags` to its corresponding obs level's values of `flagged` and `flags`.
  ///
  /// \param[in] locsExt: indices of the profile in extended space
  /// \param[in] vector_linking_indices: vector linking model levels and obs levels
  /// \param[inout] flagged: true if location rejected, false if accepted
  /// \param[inout] flags: JEDI QC flags (indicating rejection reason)
  void setModLevelFlags(const std::vector<size_t> &locsExt,
                        const std::vector<size_t> &vector_linking_indices,
                        std::vector<bool> &flagged,
                        std::vector<int> &flags);

  /// \brief Compute average of observation values associated with each model level, and write
  ///   them into the extended space.
  /// \details The increments (observation minus background) are calculated at each obs level,
  ///   then their average is taken for each model level, using the mapping in
  ///   `vector_linking_indices`. The average includes only levels that are where-included,
  ///   QC passing, and both obs value and H(x) not missing. These average increments are then
  ///   added to the H(x) value at their corresponding model levels, and the results are written
  ///   to the extended space of `obs`.
  ///
  /// \param[in] locsOriginal: indices of original space profile
  /// \param[in] locsExt: indices of the same profile in extended space
  /// \param[in] vector_linking_indices: vector mapping obs levels to model levels
  /// \param[in] apply: true if location where-included, false if where-excluded
  /// \param[in] flags: JEDI QC flags (indicating rejection reason)
  /// \param[in] hofx: model values at obs locations
  /// \param[inout] obs: observation values
  void averageObsToModLevels(const std::vector<size_t> &locsOriginal,
                             const std::vector<size_t> &locsExt,
                             const std::vector<size_t> &vector_linking_indices,
                             const std::vector<bool> &apply,
                             const std::vector<int> &flags,
                             const std::vector<float> &hofx,
                             std::vector<float> &obs);

  /// \brief Check model vertical coordinate non-zero and goes in same direction as observation
  ///   vertical coordinate.
  /// \details Returns a bool true if the vertical coordinates are increasing (e.g. ocean depth
  ///   going surface to depth) and false if decreasing (e.g. air pressure going from surface
  ///   to top of atmosphere). Throws an exception if the model vertical coordinate is all zeros
  ///   in the extended space, or if the model and observation vertical coordinates do not both
  ///   increase or both decrease.
  ///
  /// \param[in] locsOriginal: indices of original space profile
  /// \param[in] locsExt: indices of the same profile in extended space
  /// \param[in] obs_vert_coord: observation vertical coordinate values
  /// \param[in] model_vert_coord: model vertical coordinate values (only the extended space
  ///            is used, so checks that they aren't all 0).
  bool validateVertCoords(const std::vector<size_t> &locsOriginal,
                          const std::vector<size_t> &locsExt,
                          const std::vector<float> &obs_vert_coord,
                          const std::vector<float> &model_vert_coord);

  std::string addOPSPrefix(const std::string & fullname);
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGEUTILS_H_
