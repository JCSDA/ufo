/*
 * (C) Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <set>
#include <utility>

#include "ufo/profile/ProfileDataHolder.h"

namespace ufo {
  ProfileDataHolder::ProfileDataHolder(ProfileDataHandler &profileDataHandler)
    : profileDataHandler_(profileDataHandler)
  {
    numProfileLevels_ = profileDataHandler_.getNumProfileLevels();
  }

  void ProfileDataHolder::fill(const std::vector <std::string> &variableNamesInt,
                               const std::vector <std::string> &variableNamesFloat,
                               const std::vector <std::string> &variableNamesString,
                               const oops::Variables &variableNamesGeoVaLs)
  {
    variableNamesInt_ = variableNamesInt;
    variableNamesFloat_ = variableNamesFloat;
    variableNamesString_ = variableNamesString;
    variableNamesGeoVaLs_ = variableNamesGeoVaLs;

    for (const auto& variable : variableNamesInt_)
      profileData_.emplace(variable, profileDataHandler_.get<int>(variable));
    for (const auto& variable : variableNamesFloat_)
      profileData_.emplace(variable, profileDataHandler_.get<float>(variable));
    for (const auto& variable : variableNamesString_)
      profileData_.emplace(variable, profileDataHandler_.get<std::string>(variable));
    for (const auto& variable : variableNamesGeoVaLs_)
      profileGeoVaLs_.emplace(variable, profileDataHandler_.getGeoVaLVector(variable));
  }

  std::vector <float>& ProfileDataHolder::getGeoVaLVector(const oops::Variable &var)
  {
    auto it_profileGeoVaLs = profileGeoVaLs_.find(var);
    if (it_profileGeoVaLs != profileGeoVaLs_.end()) {
      return it_profileGeoVaLs->second;
    } else {
      throw eckit::BadValue("GeoVaL " + var.name() + " not present in profile. "
                            "Please add it to the relevant argument in the call "
                            "to produceProfileVector()", Here());
    }
  }

  void ProfileDataHolder::moveValuesToHandler()
  {
    oops::Log::debug() << "  Moving values to ProfileDataHandler" << std::endl;
    for (const auto& variable : variableNamesInt_) {
      oops::Log::debug() << "   " << variable << std::endl;
      profileDataHandler_.set<int>(variable, std::move(this->get<int>(variable)));
    }
    for (const auto& variable : variableNamesFloat_) {
      oops::Log::debug() << "   " << variable << std::endl;
      profileDataHandler_.set<float>(variable, std::move(this->get<float>(variable)));
    }
    for (const auto& variable : variableNamesString_) {
      oops::Log::debug() << "   " << variable << std::endl;
      profileDataHandler_.set<std::string>(variable, std::move(this->get<std::string>(variable)));
    }
    oops::Log::debug() << std::endl;
    profileData_.clear();
  }

  void ProfileDataHolder::checkObsSpaceSection(ufo::ObsSpaceSection section)
  {
    // If extended_obs_space is not present this will throw an exception.
    const auto &extended_obs_space = this->get<int>(ufo::VariableNames::extended_obs_space);
    if (section == ufo::ObsSpaceSection::Original &&
        std::find(extended_obs_space.begin(),
                  extended_obs_space.end(), 1) != extended_obs_space.end())
      throw eckit::BadValue("This profile is expected to be in the original ObsSpace "
                            "but has been labelled as being in the extended ObsSpace.", Here());
    if (section == ufo::ObsSpaceSection::Extended &&
        std::find(extended_obs_space.begin(),
                  extended_obs_space.end(), 0) != extended_obs_space.end())
      throw eckit::BadValue("This profile is expected to be in the extended ObsSpace "
                            "but has been labelled as being in the original ObsSpace.", Here());
  }
}  // namespace ufo
