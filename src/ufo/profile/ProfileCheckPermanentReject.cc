/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/profile/ProfileCheckPermanentReject.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckPermanentReject>
  makerProfileCheckPermanentReject_("PermanentReject");

  ProfileCheckPermanentReject::ProfileCheckPermanentReject
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckPermanentReject::doPrepareOutputColumns(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Permanent rejection check" << std::endl;

    const std::vector<std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_final_reject_report,
       ufo::ProfileVariableNames::diagflags_perm_reject_report,
       ufo::ProfileVariableNames::diagflags_track_reject_report,
       ufo::ProfileVariableNames::diagflags_surplus_report,
       ufo::ProfileVariableNames::diagflags_out_of_area_report,
       ufo::ProfileVariableNames::diagflags_perm_reject_t,
       ufo::ProfileVariableNames::diagflags_perm_reject_u,
       ufo::ProfileVariableNames::diagflags_perm_reject_v,
       ufo::ProfileVariableNames::diagflags_perm_reject_rh,
       ufo::ProfileVariableNames::diagflags_perm_reject_z,
       ufo::ProfileVariableNames::diagflags_final_reject_t,
       ufo::ProfileVariableNames::diagflags_final_reject_u,
       ufo::ProfileVariableNames::diagflags_final_reject_v,
       ufo::ProfileVariableNames::diagflags_final_reject_rh,
       ufo::ProfileVariableNames::diagflags_final_reject_z};

    for (const auto &variableBool : variableNamesBool) {
      profileDataHandler.get<bool>(variableBool);
    }
  }

  void ProfileCheckPermanentReject::doRunCheck(ProfileDataHandler &profileDataHandler)
  {
    std::vector <bool> &diagFlagsReportFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_report);
    const std::vector <bool> &diagFlagsReportPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_report);
    const std::vector <bool> &diagFlagsReportTrackReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_track_reject_report);
    const std::vector <bool> &diagFlagsReportSurplus =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_surplus_report);
    const std::vector <bool> &diagFlagsReportOutOfArea =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_out_of_area_report);
    std::vector <bool> &diagFlagsTPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_t);
    std::vector <bool> &diagFlagsUPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_u);
    std::vector <bool> &diagFlagsVPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_v);
    std::vector <bool> &diagFlagsRHPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_rh);
    std::vector <bool> &diagFlagsZPermReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_perm_reject_z);
    std::vector <bool> &diagFlagsTFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_t);
    std::vector <bool> &diagFlagsUFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_u);
    std::vector <bool> &diagFlagsVFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_v);
    std::vector <bool> &diagFlagsRHFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_rh);
    std::vector <bool> &diagFlagsZFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_z);

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();

    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(diagFlagsReportPermReject,
                                         diagFlagsReportTrackReject,
                                         diagFlagsReportSurplus,
                                         diagFlagsReportOutOfArea,
                                         diagFlagsReportFinalReject)) {
      oops::Log::debug() << "At least one set of observation report diagnostic flags is empty. "
                         << "Permanent rejection check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(diagFlagsReportPermReject,
                                                    diagFlagsReportTrackReject,
                                                    diagFlagsReportSurplus,
                                                    diagFlagsReportOutOfArea,
                                                    diagFlagsReportFinalReject)
                         << std::endl;
      return;
    }

    // Set PermRejectFlag on individual elements if whole report has PermReject.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (diagFlagsReportPermReject[jlev]) {
        SetDiagnosticFlag(jlev,
                          diagFlagsTPermReject,
                          diagFlagsUPermReject,
                          diagFlagsVPermReject,
                          diagFlagsRHPermReject,
                          diagFlagsZPermReject);
      }
    }

    // Set FinalRejectFlag on individual elements if a variety of criteria
    // are met on model-level data.
    if (ModelLevels) {
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (diagFlagsReportPermReject[jlev] ||
            diagFlagsReportTrackReject[jlev] ||
            diagFlagsReportSurplus[jlev] ||
            diagFlagsReportOutOfArea[jlev]) {
          diagFlagsReportFinalReject[jlev] = true;
          SetDiagnosticFlag(jlev,
                            diagFlagsTFinalReject,
                            diagFlagsUFinalReject,
                            diagFlagsVFinalReject,
                            diagFlagsRHFinalReject,
                            diagFlagsZFinalReject);
        }
      }
    }
  }
}  // namespace ufo
