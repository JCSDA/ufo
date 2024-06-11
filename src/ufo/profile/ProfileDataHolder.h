/*
 * (C) Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATAHOLDER_H_
#define UFO_PROFILE_PROFILEDATAHOLDER_H_

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "boost/variant.hpp"

#include "oops/base/Variables.h"
#include "oops/util/missingValues.h"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  enum class ObsSpaceSection {Original, Extended};

  /// \brief Profile data holder class.
  ///
  /// \details Stores data for a single profile.
  /// Each profile is filled with a user-specified list of variables that are transferred
  /// from an associated ProfileDataHandler.
  /// A collection of ProfileDataHolders can be used in checks that act on the entire
  /// sample of profiles at once (unlike checks that act on each profile individually,
  /// in which case the ProfileDataHandler can be used).
  /// This class can also be used to modify the values in the associated ProfileDataHandler,
  /// which is necessary when updating values in the entire sample.
  class ProfileDataHolder {
   public:
    explicit ProfileDataHolder(ProfileDataHandler &profileDataHandler);

    /// Fill profile with data.
    void fill(const std::vector <std::string> &variableNamesInt,
              const std::vector <std::string> &variableNamesFloat,
              const std::vector <std::string> &variableNamesString,
              const oops::Variables &variableNamesGeoVaLs);

    /// Retrieve a vector if it is present. If not, throw an exception.
    template <typename T>
      std::vector <T>& get(const std::string &fullname)
      {
        auto it_profileData = profileData_.find(fullname);
        if (it_profileData != profileData_.end()) {
          // If the type T is incorrect then boost::get will return an exception;
          // provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (it_profileData->second);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get for " +
                                      fullname + " probably has the wrong type", Here());
          }
        } else {
          throw eckit::BadValue("Variable " + fullname + " not present in profile. "
                                "Please add it to the relevant argument in the call "
                                "to produceProfileVector()", Here());
        }
      }

    /// Retrieve a GeoVaL vector if it is present. If not, throw an exception.
    std::vector <float>& getGeoVaLVector(const oops::Variable& fullname);

    /// Set values in a vector.
    template <typename T>
      void set(const std::string &fullname, std::vector<T> &&vec_in)
      {
        // Ensure vector has expected length for this profile;
        // this may not be the case for GeoVaLs.
        vec_in.resize(this->getNumProfileLevels(), util::missingValue<T>());

        // Check whether vector is already in map.
        auto it_profileData = profileData_.find(fullname);
        if (it_profileData != profileData_.end()) {
          // Replace vector in map.
          it_profileData->second = std::move(vec_in);
        } else {
          // Add vector to map.
          profileData_.emplace(fullname, std::move(vec_in));
        }
      }

    /// Get number of profile levels for this profile.
    int getNumProfileLevels() const {return static_cast<int>(numProfileLevels_);}

    /// Move all values to the associated ProfileDataHandler.
    void moveValuesToHandler();

    /// Check this profile is in the expected ObsSpace section (original or extended).
    void checkObsSpaceSection(ufo::ObsSpaceSection section);

   private:
    /// Number of profile levels
    std::size_t numProfileLevels_;

    /// Container of each variable in the current profile.
    std::unordered_map <std::string, boost::variant
                        <std::vector <int>,
                         std::vector <float>,
                         std::vector <std::string>,
                         std::vector <bool>>> profileData_;

    /// Container of GeoVaLs in the current profile.
    std::unordered_map <oops::Variable, std::vector <float>> profileGeoVaLs_;

    /// Profile data handler
    ProfileDataHandler &profileDataHandler_;

    /// Names of int variables
    std::vector <std::string> variableNamesInt_;

    /// Names of float variables
    std::vector <std::string> variableNamesFloat_;

    /// Names of string variables
    std::vector <std::string> variableNamesString_;

    /// Names of GeoVaLs
    oops::Variables variableNamesGeoVaLs_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATAHOLDER_H_
