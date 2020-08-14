/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {
  ProfileDataHandler::ProfileDataHandler(ioda::ObsSpace &obsdb,
                                         const ProfileConsistencyCheckParameters &options,
                                         EntireSampleDataHandler &entireSampleDataHandler,
                                         const ProfileIndices &profileIndices)
    : obsdb_(obsdb),
      options_(options),
      entireSampleDataHandler_(entireSampleDataHandler),
      profileIndices_(profileIndices)
  {}

  void ProfileDataHandler::reset()
  {
    profileData_.clear();
  }

  void ProfileDataHandler::updateEntireSampleData()
  {
    for (const auto &it_profile : profileData_) {
      std::string fullname = it_profile.first;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);

      if (groupname == "QCFlags") {
        const std::vector <int>& profileData = get<int>(fullname);
        std::vector <int>& entireSampleData = entireSampleDataHandler_.get<int>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndices_.getProfileIndices()) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      }

      if (groupname == "Corrections") {
        const std::vector <float>& profileData = get<float>(fullname);
        std::vector <float>& entireSampleData = entireSampleDataHandler_.get<float>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndices_.getProfileIndices()) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      }

      if (fullname == ufo::VariableNames::name_counter_NumAnyErrors) {
        const std::vector <int>& profileData = get<int>(fullname);
        std::vector <int>& entireSampleData = entireSampleDataHandler_.get<int>(fullname);
        const size_t profileNumCurrent = profileIndices_.getProfileNumCurrent();
        const size_t entriesPerProfile = profileData.size();
        size_t idx = 0;
        for (size_t profileIndex = profileNumCurrent * entriesPerProfile;
             profileIndex < (profileNumCurrent + 1) * entriesPerProfile;
             ++profileIndex) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
        }
      }
    }
  }

  void ProfileDataHandler::setFinalReportFlags()
  {
    std::vector <int> &ReportFlags = get<int>(ufo::VariableNames::name_qc_ReportFlags);
    const std::vector <int> &NumAnyErrors = get<int>(ufo::VariableNames::name_counter_NumAnyErrors);
    if (NumAnyErrors[0] > options_.nErrorsFail.value()) {
      oops::Log::debug() << " " << NumAnyErrors[0]
                         << " errors detected, whole profile rejected" << std::endl;
      for (size_t jlev = 0; jlev < ReportFlags.size(); ++jlev) {
        ReportFlags[jlev] |= ufo::FlagsWholeObReport::FinalRejectReport;
      }
    }
  }

  void ProfileDataHandler::setFlagged(const size_t nvars,
                                      std::vector<std::vector<bool>> &flagged)
  {
    oops::Log::debug() << "Flagging observations" << std::endl;

    for (const auto& it_profile : profileData_) {
      std::string fullname = it_profile.first;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);

      if (groupname == "QCFlags") {
        const std::vector <int> &Flags = get<int>(fullname);
        size_t idx = 0;
        for (const auto& jloc : profileIndices_.getProfileIndices()) {
          // Please note this concise code relies on both FlagsElem::FinalRejectFlag
          // and FlagsWholeObReport::FinalRejectReport being equal to the same value
          // (as is the case in OPS).
          // If one or both values change, and clash with another flag in the enum,
          // this will have to be rewritten.
          if (!Flags.empty() &&
              (Flags[idx] & ufo::FlagsElem::FinalRejectFlag ||
               Flags[idx] & ufo::FlagsWholeObReport::FinalRejectReport)) {
            oops::Log::debug() << " " << jloc << std::endl;
            // Flag all variables
            for (size_t jv = 0; jv < nvars; ++jv) {
              flagged[jv][jloc] = true;
            }
          }
          idx++;
        }
      }
    }
  }
}  // namespace ufo


