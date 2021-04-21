/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {
  ProfileDataHandler::ProfileDataHandler(const ObsFilterData &data,
                                         const DataHandlerParameters &options,
                                         const std::vector <bool> &apply,
                                         const Variables &filtervars,
                                         std::vector<std::vector<bool>> &flagged)
    : obsdb_(data.obsspace()),
      geovals_(data.getGeoVaLs()),
      obsdiags_(data.getObsDiags()),
      options_(options),
      filtervars_(filtervars),
      flagged_(flagged)
  {
    profileIndices_.reset(new ProfileIndices(obsdb_, options, apply));
    entireSampleDataHandler_.reset(new EntireSampleDataHandler(obsdb_, options));
  }

  void ProfileDataHandler::resetProfileInformation()
  {
    profileData_.clear();
    GeoVaLData_.clear();
    obsDiagData_.clear();
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
          groupname == "ModelLevelsFlags" ||
          groupname == "ModelRhoLevelsFlags" ||
          fullname == ufo::VariableNames::counter_NumAnyErrors) {
        const std::vector <int>& profileData = get<int>(fullname);
        getProfileIndicesInEntireSample(groupname);
        std::vector <int>& entireSampleData = entireSampleDataHandler_->get<int>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      } else if (groupname == "Corrections" ||
                 groupname == "DerivedValue" ||
                 groupname == "ModelLevelsDerivedValue" ||
                 groupname == "ModelRhoLevelsDerivedValue") {
        const std::vector <float>& profileData = get<float>(fullname);
        getProfileIndicesInEntireSample(groupname);
        std::vector <float>& entireSampleData = entireSampleDataHandler_->get<float>(fullname);
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
    oops::Log::debug() << "Flagging observations" << std::endl;

    for (const auto& it_profile : profileData_) {
      std::string fullname = it_profile.first;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);

      if (groupname == "QCFlags") {
        oops::Log::debug() << " " << fullname << std::endl;

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
            oops::Log::debug() << "  Reject all variables, index " << profileIndex << std::endl;
            for (size_t jvar = 0; jvar < filtervars_.size(); ++jvar)
              flagged_[jvar][profileIndex] = true;
            // Move to next element in the profile.
            idxprof++;
            continue;
          }
          // Flag variable if its specific value has been rejected.
          if (idxvar < filtervars_.size() &&
              Flags[idxprof] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag) {
            oops::Log::debug() << "  Reject " << varname
                               << ", index " << profileIndex << std::endl;
            flagged_[idxvar][profileIndex] = true;
          }
          idxprof++;
        }
      }
    }
  }

  std::vector <float>& ProfileDataHandler::getGeoVaLVector(const std::string &variableName)
  {
    if (GeoVaLData_.find(variableName) != GeoVaLData_.end()) {
      // If the GeoVaL vector is already present, return it.
      return GeoVaLData_[variableName];
    } else {
      std::vector <float> vec_GeoVaL_column;
      // Only fill the GeoVaL vector if the required GeoVaLs are present
      // and there is at least one observation location.
      if (geovals_ &&
          obsdb_.nlocs() > 0 &&
          geovals_->has(variableName)) {
        // Location at which to retrieve the GeoVaL.
        // This assumes each model column for each observation in a profile is identical
        // so takes the first entry in each case.
        // todo(ctgh): this is an approximation that should be revisited
        // when considering horizontal drift.
        const size_t jloc = profileIndices_->getProfileIndices()[0];
        // Vector storing GeoVaL data for current profile.
        vec_GeoVaL_column.assign(geovals_->nlevs(variableName), 0.0);
        // Get GeoVaLs at the specified location.
        geovals_->getAtLocation(vec_GeoVaL_column, variableName, jloc);
      }
      // Add GeoVaL vector to map (even if it is empty).
      GeoVaLData_.emplace(variableName, std::move(vec_GeoVaL_column));
      return GeoVaLData_[variableName];
    }
  }

  std::vector <float>& ProfileDataHandler::getObsDiag(const std::string &fullname)
  {
    if (obsDiagData_.find(fullname) != obsDiagData_.end()) {
      // If the ObsDiag vector is already present, return it.
      return obsDiagData_[fullname];
    } else {
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);
      std::vector <float> vec_ObsDiag;
      // Attempt to retrieve variable vector from entire sample.
      // If it is not present, the vector will remain empty.
      std::vector <float> &vec_all = entireSampleDataHandler_->get<float>(fullname);
      // If the vector is empty, attempt to fill it from the ObsDiags
      // (if they are present and have the required variable).
      if (vec_all.empty() &&
          obsdiags_ &&
          obsdb_.nlocs() > 0 &&
          obsdiags_->has(varname)) {
        vec_all.assign(obsdb_.nlocs(), util::missingValue(1.0f));
        obsdiags_->get(vec_all, varname);
      }
      // If the ObsDiags vector for the entire sample is not empty,
      // fill the values for this profile.
      if (!vec_all.empty()) {
        getProfileIndicesInEntireSample(groupname);
        for (const auto& profileIndex : profileIndicesInEntireSample_)
          vec_ObsDiag.emplace_back(vec_all[profileIndex]);
      }
      // Add ObsDiag vector to map (even if it is empty).
      obsDiagData_.emplace(fullname, std::move(vec_ObsDiag));
      return obsDiagData_[fullname];
    }
  }
}  // namespace ufo
