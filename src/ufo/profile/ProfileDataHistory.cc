/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/profile/ProfileDataHistory.h"

namespace ufo {
  ProfileDataHistory::ProfileDataHistory(ProfileDataHandler &profileDataHandler)
    : ProfileDataBase(profileDataHandler)
  {}

  void ProfileDataHistory::fill()
  {
    // Fill int variables
    for (const auto& variable : {ufo::VariableNames::qcflags_air_temperature,
          ufo::VariableNames::qcflags_relative_humidity,
          ufo::VariableNames::qcflags_eastward_wind,
          ufo::VariableNames::qcflags_observation_report,
          ufo::VariableNames::ObsType})
      profileData_.emplace(variable, profileDataHandler_.get<int>(variable));

    // Fill float variables
    for (const auto& variable : {ufo::VariableNames::obs_air_temperature,
          ufo::VariableNames::obs_relative_humidity,
          ufo::VariableNames::obs_eastward_wind,
          ufo::VariableNames::Latitude,
          ufo::VariableNames::Longitude,
          ufo::VariableNames::Time})
      profileData_.emplace(variable, profileDataHandler_.get<float>(variable));

    // Fill string variables
    profileData_.emplace(ufo::VariableNames::station_ID,
                         profileDataHandler_.get<std::string>(ufo::VariableNames::station_ID));
  }
}  // namespace ufo
