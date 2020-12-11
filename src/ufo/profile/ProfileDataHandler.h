/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATAHANDLER_H_
#define UFO_PROFILE_PROFILEDATAHANDLER_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "boost/variant.hpp"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/CompareNVectors.h"
#include "oops/util/missingValues.h"

#include "ufo/profile/DataHandlerParameters.h"
#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ProfileIndices.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"
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
                       const DataHandlerParameters &options,
                       const std::vector <bool> &apply);

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
        const bool optional = options_.getOptional(groupname);
        const size_t entriesPerProfile = options_.getEntriesPerProfile(groupname);

        if (profileData_.find(fullname) != profileData_.end()) {
          // If the vector is already present, return it.
          // If the type T is incorrect then boost::get will return an exception;
          // provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (profileData_[fullname]);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get"
                                      " probably has the wrong type", Here());
          }
        } else {
          std::vector <T> vec_prof;  // Vector storing data for current profile.
          // Retrieve variable vector from entire sample.
          const std::vector <T> &vec_all = entireSampleDataHandler_->get<T>(fullname);
          // Only proceed if the vector is not empty.
          if (!vec_all.empty()) {
            getProfileIndicesInEntireSample(groupname);
            for (const auto& profileIndex : profileIndicesInEntireSample_)
              vec_prof.emplace_back(vec_all[profileIndex]);
          }
          // Add vector to map (even if it is empty).
          profileData_.emplace(fullname, std::move(vec_prof));
          return boost::get<std::vector<T>> (profileData_[fullname]);
        }
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

    /// Initialise the next profile prior to applying checks.
    /// Clears \p profileData_ and determines the \p profileIndices_ for the next profile.
    void initialiseNextProfile();

    /// Update information for this profile.
    /// This function calls three other functions which take the following actions:
    /// 1. Set final report flags in this profile,
    /// 2. Modify 'flagged' vector for each filter variable based on check results,
    /// 3. If any variables in the current profile were modified by the checks,
    ///    the equivalent variables in the entire sample are set to the modified values.
    void updateProfileInformation(const size_t nvars, std::vector<std::vector<bool>> &flagged);

    /// Write various quantities to the obsdb so they can be used in future QC checks.
    /// Use the method in EntireSampleDataHandler to do this.
    void writeQuantitiesToObsdb();

    /// Return obsdb
    ioda::ObsSpace &getObsdb() {return obsdb_;}

    /// Return number of levels to which QC checks should be applied.
    int getNumProfileLevels() const {return profileIndices_->getNumProfileLevels();}

   private:  // functions
    /// Reset profile information (vectors and corresponding names).
    /// This should be called every time a new profile will be retrieved.
    void resetProfileInformation();

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

    /// Transfer values from one vector to another (as long as neither is empty).
    template <typename T>
      void updateValueIfPresent(const std::vector <T> &vecIn, const size_t &idxIn,
                                std::vector <T> &vecOut, const size_t &idxOut)
      {
        // Ensure neither vector is empty.
        if (oops::anyVectorEmpty(vecIn, vecOut)) return;
        vecOut[idxOut] = vecIn[idxIn];
      }

    /// Get indices in entire sample corresponding to current profile.
    void getProfileIndicesInEntireSample(const std::string& groupname);

   private:  // members
    /// Container of each variable in the current profile.
    std::unordered_map <std::string, boost::variant
      <std::vector <int>, std::vector <float>, std::vector <std::string>>> profileData_;

    /// Observation database.
    ioda::ObsSpace &obsdb_;

    /// Configurable parameters.
    const DataHandlerParameters &options_;

    /// Class that handles the entire data sample.
    std::unique_ptr <EntireSampleDataHandler> entireSampleDataHandler_;

    /// Class that handles profile indices.
    std::unique_ptr <ProfileIndices> profileIndices_;

    /// Indices in the entire data sample that correspond to the current profile.
    std::vector <size_t> profileIndicesInEntireSample_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATAHANDLER_H_
