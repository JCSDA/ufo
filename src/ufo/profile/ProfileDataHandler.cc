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
                                         const DataHandlerParameters &options,
                                         const std::vector <bool> &apply)
    : obsdb_(obsdb),
      options_(options)
  {
    profileIndices_.reset(new ProfileIndices(obsdb, options, apply));
    entireSampleDataHandler_.reset(new EntireSampleDataHandler(obsdb, options));
  }

  void ProfileDataHandler::resetProfileInformation()
  {
    profileData_.clear();
  }

  void ProfileDataHandler::initialiseNextProfile()
  {
    resetProfileInformation();
    profileIndices_->updateNextProfileIndices();
  }

  void ProfileDataHandler::updateProfileInformation(const size_t nvars,
                                                    std::vector<std::vector<bool>> &flagged)
  {
    // Set final report flags in this profile.
    setFinalReportFlags();

    // Modify 'flagged' vector for each filter variable based on check results.
    setFlagged(nvars, flagged);

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

      if (groupname == "QCFlags" || fullname == ufo::VariableNames::counter_NumAnyErrors) {
        getProfileIndicesInEntireSample(groupname);
        const std::vector <int>& profileData = get<int>(fullname);
        std::vector <int>& entireSampleData = entireSampleDataHandler_->get<int>(fullname);
        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          updateValueIfPresent(profileData, idx, entireSampleData, profileIndex);
          idx++;
        }
      } else if (groupname == "Corrections") {
        getProfileIndicesInEntireSample(groupname);
        const std::vector <float>& profileData = get<float>(fullname);
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
        oops::Log::debug() << " " << fullname << std::endl;

        const std::vector <int> &Flags = get<int>(fullname);
        getProfileIndicesInEntireSample(groupname);

        size_t idx = 0;
        for (const auto& profileIndex : profileIndicesInEntireSample_) {
          // Please note this concise code relies on both FlagsElem::FinalRejectFlag
          // and FlagsWholeObReport::FinalRejectReport being equal to the same value
          // (as is the case in OPS).
          // If one or both values change, and clash with another flag in the enum,
          // this will have to be rewritten.
          if (!Flags.empty() &&
              (Flags[idx] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag ||
               Flags[idx] & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport)) {
            oops::Log::debug() << "  " << profileIndex << std::endl;
            // Flag all variables
            for (size_t jv = 0; jv < nvars; ++jv)
              flagged[jv][profileIndex] = true;
          }
          idx++;
        }
      }
    }
  }
}  // namespace ufo


