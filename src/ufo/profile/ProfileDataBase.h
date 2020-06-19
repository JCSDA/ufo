/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATABASE_H_
#define UFO_PROFILE_PROFILEDATABASE_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/CompareNVectors.h"

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/profile/ProfileIndices.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

  /// \brief Store data and manipulate into individual profiles
  class ProfileDataBase {
   public:
    ProfileDataBase(ioda::ObsSpace &obsdb,
                    const ProfileConsistencyCheckParameters &options,
                    const ProfileIndices &profileIndices);
    virtual ~ProfileDataBase() {}

    /// Set number of a particular profile
    void setProfileNum(const size_t jprof) {jprof_ = jprof;}

   protected:  // functions
    /// Fill profile vector for a particular variable
    template <typename T>
      void fillProfileData(const std::vector <T> &vec_all,
                           std::vector <T> &vec_prof)
      {
        vec_prof.clear();
        if (vec_all.size() > 0) {
          for (auto profileIndex : profileIndices_.getProfileIndices()) {
            vec_prof.emplace_back(vec_all[profileIndex]);
          }
        }
      }

    /// Retrieve numerical data from obsdb
    template <typename T>
      void retrieveDataVector(const std::string &varname,
                              const std::string &groupname,
                              std::vector <T> &datavec,
                              bool optional = false)
      {
        if (obsdb_.has(groupname, varname)) {
          datavec.assign(obsdb_.nlocs(), 0);
          obsdb_.get_db(groupname, varname, datavec);
        } else if (optional) {
          datavec.assign(obsdb_.nlocs(), 0);
        }
      }

    /// Retrieve string data from obsdb
    void retrieveDataVector(const std::string &varname,
                            const std::string &groupname,
                            std::vector <std::string> &datavec,
                            bool optional = false)
    {
      if (obsdb_.has(groupname, varname)) {
        datavec.assign(obsdb_.nlocs(), "");
        obsdb_.get_db(groupname, varname, datavec);
      } else if (optional) {
        datavec.assign(obsdb_.nlocs(), "");
      }
    }

    /// Retrieve counter from obsdb
    template <typename T>
      void retrieveCounterVector(const std::string &varname,
                                 const std::string &groupname,
                                 std::vector <T> &countervec)
      {
        countervec.assign(obsdb_.nrecs(), 0);
        if (obsdb_.has(groupname, varname)) {
          obsdb_.get_db(groupname, varname, countervec);
        }
      }

    /// Put data on obsdb
    template <typename T>
      void putDataVector(const std::string &varname,
                         const std::string &groupname,
                         const std::vector <T> &datavec) const
      {
        if (datavec.empty()) return;
        obsdb_.put_db(groupname, varname, datavec);
      }

    /// Get values of all variables for entire sample
    virtual void retrieveAllData() = 0;

    /// Get values of all variables for one profile
    virtual void fillProfileValues() = 0;

   protected:  // members
    /// Observation database
    ioda::ObsSpace &obsdb_;

    /// Configurable parameters
    const ProfileConsistencyCheckParameters &options_;

    /// Indices in entire sample of observations in a particular profile
    const ProfileIndices &profileIndices_;

    /// Profile number
    size_t jprof_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATABASE_H_
