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

#include "oops/base/Variable.h"
#include "oops/util/CompareNVectors.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"

#include "ufo/profile/DataHandlerParameters.h"
#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ProfileIndices.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"
#include "ufo/utils/StringUtils.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;
  class ProfileDataHolder;
}

namespace ufo {

  /// \brief Retrieve and store data for individual profiles.
  /// To do this, first the vector of values in the entire data sample is retrieved
  /// then the relevant data corresponding to this profile are extracted.
  class ProfileDataHandler {
   public:
    ProfileDataHandler(const ObsFilterData &data,
                       ioda::ObsDataVector<int> &flags,
                       const DataHandlerParameters &options,
                       const std::vector <bool> &apply,
                       const Variables &filtervars,
                       std::vector<std::vector<bool>> &flagged);

    /// Retrieve a vector containing the requested variable for the current profile.
    ///    -# If the variable has previously been placed in a vector, return the vector.
    ///    -# Otherwise obtain the vector from the entire data sample, as long as the entire sample
    ///       is not empty.
    /// Also store the name of the variable, enabling it to be retrieved later.
    template <typename T>
      std::vector<T>& get(const std::string &fullname)
      {
        // Determine variable and group names
        std::string varname;
        std::string groupname;
        ufo::splitVarGroup(fullname, varname, groupname);

        auto it_profileData = profileData_.find(fullname);
        if (it_profileData != profileData_.end()) {
          // If the vector is already present, return it.
          // If the type T is incorrect then boost::get will return an exception;
          // provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (it_profileData->second);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get for " +
                                      fullname + " probably has the wrong type", Here());
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
    /// Also initialise a vector in the entire sample, allowing the data to
    /// be stored between checks.
    template <typename T>
      void set(const std::string &fullname, std::vector<T> &&vec_in)
      {
        // Determine variable and group names
        std::string varname;
        std::string groupname;
        ufo::splitVarGroup(fullname, varname, groupname);
        // Check whether vector is already in map.
        auto it_profileData = profileData_.find(fullname);
        if (it_profileData != profileData_.end()) {
          // Replace vector in map.
          it_profileData->second = std::move(vec_in);
        } else {
          // Add vector to map.
          profileData_.emplace(fullname, std::move(vec_in));
        }
        entireSampleDataHandler_->initialiseVector<T>(fullname);
        // Transfer this profile's data into the entire sample.
        getProfileIndicesInEntireSample(groupname);
        std::vector <T>& entireSampleData = entireSampleDataHandler_->get<T>(fullname);
        const std::vector <T>& profileData = this->get<T>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
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
    void updateProfileInformation();

    /// Write various quantities to the obsdb so they can be used in future QC checks.
    /// Use the method in EntireSampleDataHandler to do this.
    void writeQuantitiesToObsdb();

    /// Return obsdb
    ioda::ObsSpace &getObsdb() {return obsdb_;}

    /// Return number of levels to which QC checks should be applied.
    int getNumProfileLevels() const {return profileIndices_->getNumProfileLevels();}

    /// Get GeoVaLs for a particular profile.
    std::vector <float>& getGeoVaLVector(const oops::Variable &variable);

    /// Get filter flags
    ioda::ObsDataVector<int>& getFilterFlags() const {return flags_;}

    /// Reset profile indices (required if it is desired to loop through
    /// the entire sample again).
    void resetProfileIndices() {profileIndices_->reset();}

    /// Produce a vector of all profiles, loading the requested variables into each one.
    std::vector <ProfileDataHolder> produceProfileVector
      (const std::vector <std::string> &variableNamesInt,
       const std::vector <std::string> &variableNamesFloat,
       const std::vector <std::string> &variableNamesString,
       const oops::Variables &variableNamesGeoVaLs);

    /// Read values from a collection of profiles and update information related to each one.
    void updateAllProfiles(std::vector <ProfileDataHolder> &profiles);

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
    void setFlagged();

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

    /// Get the name of the vertical coordinate that is used to determine the slant path
    /// locations for the variable \p variableName.
    oops::Variable getAssociatedVerticalCoordinate(const oops::Variable & variable) const;

   private:  // members
    /// Container of each variable in the current profile.
    std::unordered_map <std::string, boost::variant
                        <std::vector <int>,
                         std::vector <float>,
                         std::vector <std::string>,
                         std::vector <bool>>> profileData_;

    /// Container of GeoVaLs in the current profile.
    std::unordered_map <oops::Variable, std::vector <float>> GeoVaLData_;

    /// Observation database.
    ioda::ObsSpace &obsdb_;

    /// GeoVaLs.
    std::unique_ptr<GeoVaLs> geovals_;

    /// Filter flags
    ioda::ObsDataVector<int> &flags_;

    /// Configurable parameters.
    const DataHandlerParameters &options_;

    /// Filter variables
    const Variables &filtervars_;

    /// Flagged values
    std::vector<std::vector<bool>> &flagged_;

    /// Class that handles the entire data sample.
    std::unique_ptr <EntireSampleDataHandler> entireSampleDataHandler_;

    /// Class that handles profile indices.
    std::unique_ptr <ProfileIndices> profileIndices_;

    /// Indices in the entire data sample that correspond to the current profile.
    std::vector <size_t> profileIndicesInEntireSample_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATAHANDLER_H_
