/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATAHANDLER_H_
#define UFO_PROFILE_PROFILEDATAHANDLER_H_

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "boost/variant.hpp"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/CompareNVectors.h"

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ProfileIndices.h"

#include "ufo/utils/Flags.h"
#include "ufo/utils/StringUtils.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

  /// \brief Retrieve and store data for individual profiles.
  /// To do this, first the vector of values in the entire data sample is retrieved
  /// then the relevant data corresponding to this profile are extracted.
  class ProfileDataHandler {
   public:
    ProfileDataHandler(ioda::ObsSpace &obsdb,
                       const ProfileConsistencyCheckParameters &options,
                       EntireSampleDataHandler &entireSampleDataHandler,
                       const ProfileIndices &profileIndices);

    /// Retrieve a vector containing the requested variable for the current profile.
    ///    -# If the variable has previously been placed in a vector, return the vector.
    ///    -# Otherwise obtain the vector from the entire data sample, as long as the entire sample
    ///       is not empty.
    /// Also store the name of the variable, enabling it to be retrieved later.
    template <typename T>
      std::vector<T>& get(const std::string &fullname)
      {
        // Determine variable and group names, optional, and number of entries per profile.
        std::string varname;
        std::string groupname;
        ufo::splitVarGroup(fullname, varname, groupname);
        bool optional = options_.getOptional(groupname);
        size_t entriesPerProfile = options_.getEntriesPerProfile(groupname);

        std::vector <T> vec_prof;  // Vector storing data for current profile.
        if (profileData_.find(fullname) != profileData_.end()) {
          // If the vector is already present, return it.
          // If the type T is incorrect then boost::get will return an exception.
          // Provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (profileData_[fullname]);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get"
                                      " probably has the wrong type", Here());
          }
        } else {
          // Retrieve variable vector from entire sample.
          const std::vector <T> &vec_all = entireSampleDataHandler_.get<T>(fullname);
          // Only proceed if the vector is not empty.
          if (vec_all.size() > 0) {
            // If the number of entries per profile was not specified use the indices
            // that were obtained by sorting and grouping the record numbers.
            if (entriesPerProfile == -1) {
              for (const auto& profileIndex : profileIndices_.getProfileIndices()) {
                vec_prof.emplace_back(vec_all[profileIndex]);
              }
            } else {
              // Otherwise, loop over the relevant portion of the entire sample.
              size_t profileNumCurrent = profileIndices_.getProfileNumCurrent();
              for (size_t profileIndex = profileNumCurrent * entriesPerProfile;
                   profileIndex < (profileNumCurrent + 1) * entriesPerProfile;
                   ++profileIndex) {
                vec_prof.emplace_back(vec_all[profileIndex]);
              }
            }
          }
        }
        // Add vector to map (even if it is empty).
        profileData_.emplace(fullname, std::move(vec_prof));
        return boost::get <std::vector<T>> (profileData_[fullname]);
      }

    /// Directly set a vector for the current profile.
    /// Typically used to store variables that are used locally in checks
    /// (e.g. intermediate values).
    template <typename T>
      void set(const std::string &fullname, std::vector<T> &&vec_in)
      {
        // Check whether vector is already in map.
        if (profileData_.find(fullname) != profileData_.end()) {
          // Replace vector in map.
          profileData_[fullname] = vec_in;
        } else {
          // Add vector to map.
          profileData_.emplace(fullname, vec_in);
        }
      }

    /// Reset profile information (vectors and corresponding names).
    /// This should be called every time a new profile will be retrieved.
    void reset();

    /// Transfer values from one vector to another (as long as neither is empty).
    template <typename T>
      void updateValueIfPresent(const std::vector <T> &vecIn, const size_t &idxIn,
                                std::vector <T> &vecOut, const size_t &idxOut)
      {
        // Ensure neither vector is empty.
        if (oops::anyVectorEmpty(vecIn, vecOut)) return;
        vecOut[idxOut] = vecIn[idxIn];
      }

    /// If any variables in the current profile were modified by the checks,
    /// the equivalent variables in the entire sample are set to the modified values.
    /// The variables that are (potentially) modified are hardcoded but this could be
    /// changed to a configurable list if requred.
    void updateEntireSampleData();

    /// Set final report flags based on the NumAnyErrors counter.
    void setFinalReportFlags();

    /// Update the 'flagged' vector based on any flags that may have changed during the checks.
    /// The QC flag group is hardocded but this could be changed to
    /// a configurable value if required.
    void setFlagged(const size_t nvars, std::vector<std::vector<bool>> &flagged);

   private:
    /// Container of each variable in the current profile.
    std::unordered_map <std::string, boost::variant
      <std::vector <int>, std::vector <float>, std::vector <std::string>>> profileData_;

    /// Observation database.
    ioda::ObsSpace &obsdb_;

    /// Configurable parameters.
    const ProfileConsistencyCheckParameters &options_;

    /// Class that handles the entire data sample.
    EntireSampleDataHandler &entireSampleDataHandler_;

    /// Indices in entire sample of observations in a particular profile.
    const ProfileIndices &profileIndices_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATAHANDLER_H_
