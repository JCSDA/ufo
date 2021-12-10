/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/utils/StringTools.h"

#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"

#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"
#include "ufo/profile/SlantPathLocations.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {
  ProfileDataHandler::ProfileDataHandler(const ObsFilterData &data,
                                         ioda::ObsDataVector<int> &flags,
                                         const DataHandlerParameters &options,
                                         const std::vector <bool> &apply,
                                         const Variables &filtervars,
                                         std::vector<std::vector<bool>> &flagged)
    : obsdb_(data.obsspace()),
      flags_(flags),
      options_(options),
      filtervars_(filtervars),
      flagged_(flagged)
  {
    if (data.getGeoVaLs() && data.getGeoVaLs()->nlocs() > 0) {
      geovals_.reset(new GeoVaLs(*(data.getGeoVaLs())));
      if (geovals_->has(ufo::VariableNames::geovals_pressure))
        geovals_->reorderzdir(ufo::VariableNames::geovals_pressure, "bottom2top");
      else if (geovals_->has(ufo::VariableNames::geovals_pressure_rho))
        geovals_->reorderzdir(ufo::VariableNames::geovals_pressure_rho, "bottom2top");
      else
        throw eckit::BadValue("GeoVaLs must contain a pressure coordinate");
    }
    profileIndices_.reset(new ProfileIndices(obsdb_, options, apply));
    entireSampleDataHandler_.reset(new EntireSampleDataHandler(data, options));
  }

  void ProfileDataHandler::resetProfileInformation()
  {
    profileData_.clear();
    GeoVaLData_.clear();
  }

  void ProfileDataHandler::initialiseNextProfile()
  {
    resetProfileInformation();
    profileIndices_->updateNextProfileIndices();
  }

  void ProfileDataHandler::updateProfileInformation()
  {
    // Set final report flags in this profile.
    setFinalReportFlags();

    // Modify 'flagged' vector for each filter variable based on check results.
    setFlagged();

    // If any variables in the current profile were modified by the checks,
    // the equivalent variables in the entire sample are set to the modified values.
    updateEntireSampleData();
  }

  void ProfileDataHandler::writeQuantitiesToObsdb()
  {
    entireSampleDataHandler_->writeQuantitiesToObsdb();
  }

  void ProfileDataHandler::getProfileIndicesInEntireSample(const std::string& groupname)
  {
    profileIndicesInEntireSample_.clear();
    const size_t entriesPerProfile = options_.getEntriesPerProfile(groupname);
    // If the number of entries per profile was not specified, use the indices
    // that were obtained by sorting and grouping the record numbers.
    if (entriesPerProfile == 0) {
      profileIndicesInEntireSample_ = profileIndices_->getProfileIndices();
    } else {
      // Otherwise increment the indices sequentially, starting at the
      // relevant position.
      profileIndicesInEntireSample_.resize(entriesPerProfile);
      std::iota(profileIndicesInEntireSample_.begin(),
                profileIndicesInEntireSample_.end(),
                profileIndices_->getProfileNumCurrent() * entriesPerProfile);
    }
  }

  void ProfileDataHandler::updateEntireSampleData()
  {
    for (const auto &it_profile : profileData_) {
      std::string fullname = it_profile.first;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);

      if (groupname == "QCFlags" ||
          groupname == "Counters") {
        const std::vector <int>& profileData = get<int>(fullname);
        getProfileIndicesInEntireSample(groupname);
        std::vector <int>& entireSampleData = entireSampleDataHandler_->get<int>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      } else if (groupname == "Corrections" ||
                 groupname == "DerivedObsValue" ||
                 groupname == "DerivedMetaData" ||
                 groupname == "DerivedModelValue" ||
                 groupname == "GrossErrorProbability" ||
                 fullname == ufo::VariableNames::obs_air_pressure) {
        const std::vector <float>& profileData = get<float>(fullname);
        getProfileIndicesInEntireSample(groupname);
        std::vector <float>& entireSampleData = entireSampleDataHandler_->get<float>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      } else if (eckit::StringTools::startsWith(groupname, "DiagnosticFlags")) {
        const std::vector <bool>& profileData = get<bool>(fullname);
        getProfileIndicesInEntireSample(groupname);
        std::vector <bool>& entireSampleData = entireSampleDataHandler_->get<bool>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      }
    }
  }

  void ProfileDataHandler::setFinalReportFlags()
  {
    std::vector <int> &ReportFlags = get<int>(ufo::VariableNames::qcflags_observation_report);
    const std::vector <int> &NumAnyErrors = get<int>(ufo::VariableNames::counter_NumAnyErrors);
    if (!NumAnyErrors.empty() && NumAnyErrors[0] > options_.nErrorsFail.value()) {
      oops::Log::debug() << " " << NumAnyErrors[0]
                         << " errors detected, whole profile rejected" << std::endl;
      for (size_t jlev = 0; jlev < ReportFlags.size(); ++jlev) {
        ReportFlags[jlev] |= ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport;
      }
    }
  }

  void ProfileDataHandler::setFlagged()
  {
    oops::Log::debug() << "   Flagging observations" << std::endl;

    for (const auto& it_profile : profileData_) {
      std::string fullname = it_profile.first;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);

      if (groupname == "QCFlags") {
        oops::Log::debug() << "    " << fullname << std::endl;

        // Obtain QC flags
        const std::vector <int> &Flags = get<int>(fullname);
        if (Flags.empty()) continue;
        getProfileIndicesInEntireSample(groupname);

        // Determine index of varname in the filter variables.
        // If it is not present then the variable will not be flagged individually.
        size_t idxvar = filtervars_.size();
        for (size_t idx = 0; idx < idxvar; ++idx) {
          if (filtervars_[idx].variable() == varname) {
            idxvar = idx;
            break;
          }
        }

        // If varname is observation_report then all filter variables will be rejected.
        bool isObservationReport = varname == "observation_report";

        // Index of elements in this profile.
        size_t idxprof = 0;
        // Loop over indices of elements in entire profile sample.
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          // Flag all filter variables if the whole observation has been rejected.
          if (isObservationReport &&
              Flags[idxprof] & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport) {
            oops::Log::debug() << "     Reject all variables, index " << profileIndex << std::endl;
            for (size_t jvar = 0; jvar < filtervars_.size(); ++jvar)
              flagged_[jvar][profileIndex] = true;
            // Move to next element in the profile.
            idxprof++;
            continue;
          }
          // Flag variable if its specific value has been rejected.
          if (idxvar < filtervars_.size() &&
              Flags[idxprof] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag) {
            oops::Log::debug() << "     Reject " << varname
                               << ", index " << profileIndex << std::endl;
            flagged_[idxvar][profileIndex] = true;
          }
          idxprof++;
        }
      }
    }
  }

  std::string ProfileDataHandler::getAssociatedVerticalCoordinate
  (const std::string & variableName) const
  {
    // Obtain the map between non-default vertical coordinates and variable names.
    const auto & alternativeVerticalCoordinate = options_.alternativeVerticalCoordinate.value();
    auto it_altCoord = alternativeVerticalCoordinate.find(variableName);
    if (it_altCoord != alternativeVerticalCoordinate.end()) {
      // This variable has an associated alternative vertical coordinate.
      return it_altCoord->second;
    } else {
      // This variable uses the default vertical coordinate.
      return options_.defaultVerticalCoordinate.value();
    }
  }

  std::vector <float>& ProfileDataHandler::getGeoVaLVector(const std::string &variableName)
  {
    auto it_GeoVaLData = GeoVaLData_.find(variableName);
    if (it_GeoVaLData != GeoVaLData_.end()) {
      // If the GeoVaL vector is already present, return it.
      return it_GeoVaLData->second;
    } else {
      std::vector <float> vec_GeoVaL_column;
      // Only fill the GeoVaL vector if the required GeoVaLs are present
      // and there is at least one observation location.
      if (geovals_ &&
          obsdb_.nlocs() > 0 &&
          geovals_->has(variableName)) {
        // Locations at which to retrieve the GeoVaL.
        const std::vector<std::size_t> slant_path_location =
          ufo::getSlantPathLocations(obsdb_,
                                     *geovals_,
                                     profileIndices_->getProfileIndices(),
                                     this->getAssociatedVerticalCoordinate(variableName));
        // Vector storing GeoVaL data for current profile.
        vec_GeoVaL_column.assign(geovals_->nlevs(variableName), 0.0);
        // Check the number of entries in the slant path location vector is equal
        // to the number of entries in the GeoVaL for this variable.
        // If not, the GeoVaL at the first location in the profile is used;
        // in other words, drift is not accounted for.
        // todo(ctgh): revisit this choice in a future PR.
        if (slant_path_location.size() == vec_GeoVaL_column.size()) {
          std::vector<float> vec_GeoVaL_loc(geovals_->nlevs(variableName));
          // Take the GeoVaL at each slant path location and copy the relevant
          // value from each GeoVaL into the output vector.
          for (std::size_t mlev = 0; mlev < slant_path_location.size(); ++mlev) {
            const std::size_t jloc = slant_path_location[mlev];
            geovals_->getAtLocation(vec_GeoVaL_loc, variableName, jloc);
            vec_GeoVaL_column[mlev] = vec_GeoVaL_loc[mlev];
          }
        } else {
          // Take the GeoVaL at the first location.
          const std::size_t jloc = profileIndices_->getProfileIndices()[0];
          geovals_->getAtLocation(vec_GeoVaL_column, variableName, jloc);
        }
      }
      // Add GeoVaL vector to map (even if it is empty).
      GeoVaLData_.emplace(variableName, std::move(vec_GeoVaL_column));
      return GeoVaLData_[variableName];
    }
  }

  std::vector <ProfileDataHolder> ProfileDataHandler::produceProfileVector
  (const std::vector <std::string> &variableNamesInt,
   const std::vector <std::string> &variableNamesFloat,
   const std::vector <std::string> &variableNamesString,
   const std::vector <std::string> &variableNamesGeoVaLs)
  {
    profileIndices_->reset();
    std::vector <ProfileDataHolder> profiles;
    oops::Log::debug() << "  Filling vector of profiles" << std::endl;
    for (size_t jprof = 0; jprof < obsdb_.nrecs(); ++jprof) {
      initialiseNextProfile();
      ProfileDataHolder profile(*this);
      profile.fill(variableNamesInt,
                   variableNamesFloat,
                   variableNamesString,
                   variableNamesGeoVaLs);
      profiles.emplace_back(profile);
    }
    return profiles;
  }

  void ProfileDataHandler::updateAllProfiles(std::vector <ProfileDataHolder> &profiles)
  {
    this->resetProfileIndices();
    for (size_t jprof = 0; jprof < obsdb_.nrecs(); ++jprof) {
      this->initialiseNextProfile();
      auto& profile = profiles[jprof];
      // Move values from profile to this object.
      profile.moveValuesToHandler();
      // Update information, including the 'flagged' vector, for this profile.
      this->updateProfileInformation();
    }
  }
}  // namespace ufo
