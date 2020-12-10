/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckPermanentReject.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckPermanentReject>
  makerProfileCheckPermanentReject_("PermanentReject");

  ProfileCheckPermanentReject::ProfileCheckPermanentReject
  (const ProfileConsistencyCheckParameters &options,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckPermanentReject::runCheck()
  {
    oops::Log::debug() << " Permanent rejection check" << std::endl;

    const size_t numProfileLevels = profileDataHandler_.getNumProfileLevels();
    const bool ModelLevels = options_.modellevels.value();
    std::vector <int> &tFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &rhFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    std::vector <int> &uFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &vFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_northward_wind);
    std::vector <int> &zFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_geopotential_height);
    std::vector <int> &ReportFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_observation_report);

    if (ReportFlags.empty()) {
      oops::Log::debug() << "ReportFlags vector is empty. "
                         << "Permanent rejection check will not be performed." << std::endl;
      return;
    }

    // Set PermRejectFlag on individual elements if whole report has PermReject.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (ReportFlags[jlev] & ufo::MetOfficeQCFlags::WholeObReport::PermRejectReport) {
        SetQCFlag(ufo::MetOfficeQCFlags::Elem::PermRejectFlag, jlev,
                  tFlags, rhFlags, uFlags, vFlags, zFlags);
      }
    }

    // Set FinalRejectFlag on individual elements if a variety of criteria
    // are met on model-level data.
    if (ModelLevels) {
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if ((ReportFlags[jlev] & ufo::MetOfficeQCFlags::WholeObReport::PermRejectReport) ||
            (ReportFlags[jlev] & ufo::MetOfficeQCFlags::WholeObReport::TrackRejectReport) ||
            (ReportFlags[jlev] & ufo::MetOfficeQCFlags::WholeObReport::SurplusReport) ||
            (ReportFlags[jlev] & ufo::MetOfficeQCFlags::WholeObReport::OutOfAreaReport)) {
          ReportFlags[jlev] |= ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport;
          SetQCFlag(ufo::MetOfficeQCFlags::Elem::FinalRejectFlag, jlev,
                    tFlags, rhFlags, uFlags, vFlags, zFlags);
        }
      }
    }
  }
}  // namespace ufo
