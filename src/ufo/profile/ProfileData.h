/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEDATA_H_
#define UFO_PROFILE_PROFILEDATA_H_

#include <string>
#include <vector>

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/profile/ProfileDataBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ProfileIndices;
}

namespace ufo {

  /// \brief Store data and manipulate into individual profiles
  class ProfileData : public ProfileDataBase {
   public:
    ProfileData(ioda::ObsSpace &obsdb,
                const ProfileConsistencyCheckParameters &options,
                const ProfileIndices &profileIndices);

    /// Fill values of all variables for one profile
    void fillProfileValues() override;

    /// Return profile pressures
    std::vector <float> getPressures() const {return pressures_prof_;}

    /// Return profile tObs
    std::vector <float> gettObs() const {return tObs_prof_;}

    /// Return profile tBkg
    std::vector <float> gettBkg() const {return tBkg_prof_;}

    /// Return profile zObs
    std::vector <float> getzObs() const {return zObs_prof_;}

    /// Return profile zBkg
    std::vector <float> getzBkg() const {return zBkg_prof_;}

    /// Return profile PstarBackgr
    std::vector <float> getPstarBackgr() const {return PstarBackgr_prof_;}

    /// Return profile station ID
    std::string getStationID() const
      {
        if (stationID_prof_.size() > 0)
          return stationID_prof_[0];
        else
          return "";
      }

   private:
    /// Retrieve values of all variables for entire sample
    void retrieveAllData() override;

    //=== Entire sample values ===//

    /// Entire sample pressures
    std::vector <float> pressures_;

    /// Entire sample tObs
    std::vector <float> tObs_;

    /// Entire sample tBkg
    std::vector <float> tBkg_;

    /// Entire sample zObs
    std::vector <float> zObs_;

    /// Entire sample zBkg
    std::vector <float> zBkg_;

    /// Entire sample PstarBackgr
    std::vector <float> PstarBackgr_;

    /// Entire sample station ID
    std::vector <std::string> stationID_;

    //=== Individual profile values ===//

    /// Individual profile pressures
    std::vector <float> pressures_prof_;

    /// Individual profile tObs
    std::vector <float> tObs_prof_;

    /// Individual profile tBkg
    std::vector <float> tBkg_prof_;

    /// Individual profile zObs
    std::vector <float> zObs_prof_;

    /// Individual profile zBkg
    std::vector <float> zBkg_prof_;

    /// Individual profile PstarBackgr
    std::vector <float> PstarBackgr_prof_;

    /// Individual profile station ID
    std::vector <std::string> stationID_prof_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEDATA_H_
