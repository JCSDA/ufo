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
                               const std::vector <std::string> &variableNamesGeoVaLs,
                               const std::vector <std::string> &variableNamesObsDiags)
  {
    variableNamesInt_ = variableNamesInt;
    variableNamesFloat_ = variableNamesFloat;
    variableNamesString_ = variableNamesString;
    variableNamesGeoVaLs_ = variableNamesGeoVaLs;
    variableNamesObsDiags_ = variableNamesObsDiags;

    for (const auto& variable : variableNamesInt_)
      profileData_.emplace(variable, profileDataHandler_.get<int>(variable));
    for (const auto& variable : variableNamesFloat_)
      profileData_.emplace(variable, profileDataHandler_.get<float>(variable));
    for (const auto& variable : variableNamesString_)
      profileData_.emplace(variable, profileDataHandler_.get<std::string>(variable));
    for (const auto& variable : variableNamesGeoVaLs_)
      profileGeoVaLs_.emplace(variable, profileDataHandler_.getGeoVaLVector(variable));
    for (const auto& variable : variableNamesObsDiags_)
      profileObsDiags_.emplace(variable, profileDataHandler_.getObsDiag(variable));
  }

  std::vector <float>& ProfileDataHolder::getGeoVaLVector(const std::string &fullname)
  {
    auto it_profileGeoVaLs = profileGeoVaLs_.find(fullname);
    if (it_profileGeoVaLs != profileGeoVaLs_.end()) {
      return it_profileGeoVaLs->second;
    } else {
      throw eckit::BadValue("GeoVaL " + fullname + " not present in profile. "
                            "Please add it to the relevant argument in the call "
                            "to produceProfileVector()", Here());
    }
  }

  std::vector <float>& ProfileDataHolder::getObsDiagVector(const std::string &fullname)
  {
    auto it_profileObsDiags = profileObsDiags_.find(fullname);
    if (it_profileObsDiags != profileObsDiags_.end()) {
      return it_profileObsDiags->second;
    } else {
      throw eckit::BadValue("ObsDiag " + fullname + " not present in profile. "
                            "Please add it to the relevant argument in the call "
                            "to produceProfileVector()", Here());
    }
  }

  void ProfileDataHolder::moveValuesToHandler()
  {
    for (const auto& variable : variableNamesInt_)
      profileDataHandler_.set<int>(variable, std::move(this->get<int>(variable)));
    for (const auto& variable : variableNamesFloat_)
      profileDataHandler_.set<float>(variable, std::move(this->get<float>(variable)));
    for (const auto& variable : variableNamesString_)
      profileDataHandler_.set<std::string>(variable, std::move(this->get<std::string>(variable)));
    profileData_.clear();
  }
}  // namespace ufo
