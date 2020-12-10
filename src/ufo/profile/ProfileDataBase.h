/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATABASE_H_
#define UFO_PROFILE_PROFILEDATABASE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "boost/variant.hpp"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  /// \brief Profile data base class
  class ProfileDataBase {
   public:
    explicit ProfileDataBase(ProfileDataHandler &profileDataHandler);
    virtual ~ProfileDataBase() {}

    /// Fill profile with data
    virtual void fill() = 0;

    /// Retrieve a vector if it is present. If not, throw an exception;
    /// the variable should be requested in the fill() routine.
    template <typename T>
      std::vector <T>& get(const std::string &fullname)
      {
        if (profileData_.find(fullname) != profileData_.end()) {
          // If the type T is incorrect then boost::get will return an exception;
          // provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (profileData_[fullname]);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get"
                                      " probably has the wrong type", Here());
          }
        } else {
          throw eckit::BadValue("Variable not present in profile. "
                                "Please add it to the fill() routine", Here());
        }
      }

   protected:
    /// Container of each variable in the current profile.
    std::unordered_map <std::string, boost::variant
      <std::vector <int>, std::vector <float>, std::vector <std::string>>> profileData_;

    /// Profile data handler
    ProfileDataHandler &profileDataHandler_;

    // to add - list of indices for each profile
    // (might need them when writing back to the entire sample)
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATABASE_H_
