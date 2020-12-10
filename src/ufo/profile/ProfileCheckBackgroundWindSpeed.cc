/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBackgroundWindSpeed.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundWindSpeed>
  makerProfileCheckBackgroundWindSpeed_("BackgroundWindSpeed");

  ProfileCheckBackgroundWindSpeed::ProfileCheckBackgroundWindSpeed
  (const ProfileConsistencyCheckParameters &options,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckBackgroundWindSpeed::runCheck()
  {
    oops::Log::debug() << " Background check for wind velocity" << std::endl;

    const size_t numProfileLevels = profileDataHandler_.getNumProfileLevels();
    const bool ModelLevels = options_.modellevels.value();
    const std::vector <float> &uObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_eastward_wind);
    const std::vector <float> &uObsErr =
       profileDataHandler_.get<float>(ufo::VariableNames::obserr_eastward_wind);
    const std::vector <float> &uBkg =
      profileDataHandler_.get<float>(ufo::VariableNames::hofx_eastward_wind);
    const std::vector <float> &uBkgErr =
      profileDataHandler_.get<float>(ufo::VariableNames::bkgerr_eastward_wind);
    std::vector <float> &uPGE =
      profileDataHandler_.get<float>(ufo::VariableNames::pge_eastward_wind);
    std::vector <float> &uPGEBd =
      profileDataHandler_.get<float>(ufo::VariableNames::pgebd_eastward_wind);
    std::vector <int> &uFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    const std::vector <float> &vObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_northward_wind);
    const std::vector <float> &vObsErr =
       profileDataHandler_.get<float>(ufo::VariableNames::obserr_northward_wind);
    const std::vector <float> &vBkg =
      profileDataHandler_.get<float>(ufo::VariableNames::hofx_northward_wind);
    const std::vector <float> &vBkgErr =
      profileDataHandler_.get<float>(ufo::VariableNames::bkgerr_northward_wind);
    std::vector <float> &vPGE =
      profileDataHandler_.get<float>(ufo::VariableNames::pge_northward_wind);
    std::vector <float> &vPGEBd =
      profileDataHandler_.get<float>(ufo::VariableNames::pgebd_northward_wind);
    std::vector <int> &vFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_northward_wind);
    const std::vector <int> &timeFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_time);

    if (!oops::allVectorsSameNonZeroSize(uObs, uObsErr, uBkg, uBkgErr,
                                         uPGE, uFlags,
                                         vObs, vObsErr, vBkg, vBkgErr,
                                         vPGE, vFlags, timeFlags)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(uObs, uObsErr, uBkg, uBkgErr,
                                                      uPGE, uFlags,
                                                      vObs, vObsErr, vBkg, vBkgErr,
                                                      vPGE, vFlags, timeFlags)
                           << std::endl;
      return;
    }

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, options_.BkCheck_PdBad_uv.value());

    // Modify observation PGE if certain flags have been set.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (uFlags[jlev] & ufo::MetOfficeQCFlags::Profile::InterpolationFlag)
        uPGE[jlev] = 0.5 + 0.5 * uPGE[jlev];
      if (timeFlags[jlev])
        uFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           uObs,
                           uObsErr,
                           uBkg,
                           uBkgErr,
                           PdBad,
                           ModelLevels,
                           uFlags,
                           uPGE,
                           uPGEBd,
                           -1,
                           &vObs,
                           &vBkg);

    // Update v PGE and flags.
    vPGE = uPGE;
    vFlags = uFlags;
    vPGEBd = uPGEBd;
  }
}  // namespace ufo
