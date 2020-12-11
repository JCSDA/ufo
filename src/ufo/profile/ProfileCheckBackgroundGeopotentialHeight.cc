/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBackgroundGeopotentialHeight.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundGeopotentialHeight>
  makerProfileCheckBackgroundGeopotentialHeight_("BackgroundGeopotentialHeight");

  ProfileCheckBackgroundGeopotentialHeight::ProfileCheckBackgroundGeopotentialHeight
  (const ProfileConsistencyCheckParameters &options,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckBackgroundGeopotentialHeight::runCheck()
  {
    oops::Log::debug() << " Background check for geopotential height" << std::endl;

    const size_t numProfileLevels = profileDataHandler_.getNumProfileLevels();
    const bool ModelLevels = options_.modellevels.value();
    const std::vector <float> &Zstation =
      profileDataHandler_.get<float>(ufo::VariableNames::Zstation);
    const std::vector <float> &pressures =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &zObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_geopotential_height);
    const std::vector <float> &zObsErr =
       profileDataHandler_.get<float>(ufo::VariableNames::obserr_geopotential_height);
    const std::vector <float> &zBkg =
      profileDataHandler_.get<float>(ufo::VariableNames::hofx_geopotential_height);
    std::vector <float> &zBkgErr =
      profileDataHandler_.get<float>(ufo::VariableNames::bkgerr_geopotential_height);
    std::vector <float> &zPGE =
      profileDataHandler_.get<float>(ufo::VariableNames::pge_geopotential_height);
    std::vector <float> &zPGEBd =
      profileDataHandler_.get<float>(ufo::VariableNames::pgebd_geopotential_height);
    std::vector <int> &zFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_geopotential_height);
    const std::vector <int> &tFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_air_temperature);
    const std::vector <float> &zObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::obscorrection_geopotential_height);
    const std::vector <int> &timeFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_time);

    if (!oops::allVectorsSameNonZeroSize(Zstation, pressures,
                                         zObs, zObsErr, zBkg,
                                         zPGE, zFlags, zObsCorrection,
                                         tFlags, timeFlags)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(Zstation, pressures,
                                                      zObs, zObsErr, zBkg,
                                                      zPGE, zFlags, zObsCorrection,
                                                      tFlags, timeFlags)
                           << std::endl;
      return;
    }

    std::vector <float> zObsFinal;
    correctVector(zObs, zObsCorrection, zObsFinal);

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, 0);
    // The z background error may not have been set before this point.
    if (zBkgErr.empty())
      zBkgErr.assign(numProfileLevels, missingValueFloat);
    // Background error estimates are taken from ECMWF Research Manual 1
    // (ECMWF Data Assimilation Scientific Documentation, 3/92, 3rd edition), table 2.1.
    // They are then multiplied by 1.2.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (zObsFinal[jlev] != missingValueFloat) {
        // Permanently reject any levels at/below surface
        if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag ||
            zObsFinal[jlev] <= Zstation[jlev]) {
          zFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
        }
        const float Plevel = pressures[jlev] / 100.0;  // hPa
        const std::vector<float> zBkgErrs = options_.BkCheck_zBkgErrs.value();
        const std::vector<float> zBadPGEs = options_.BkCheck_zBadPGEs.value();
        size_t idx = 0;
        for (const auto& PlevelThreshold : options_.BkCheck_PlevelThresholds.value()) {
          if (Plevel > PlevelThreshold) {
            zBkgErr[jlev] = 1.2 * zBkgErrs[idx];
            PdBad[jlev] = zBadPGEs[idx];
            break;
          }
          idx++;
        }
      }
    }

    // Modify observation PGE if certain flags have been set.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (zFlags[jlev] & ufo::MetOfficeQCFlags::Profile::InterpolationFlag)
        zPGE[jlev] = 0.5 + 0.5 * zPGE[jlev];
      if (zFlags[jlev] & ufo::MetOfficeQCFlags::Profile::HydrostaticFlag)
        zPGE[jlev] = 0.5 + 0.5 * zPGE[jlev];
      if (timeFlags[jlev])
        zFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           zObsFinal,
                           zObsErr,
                           zBkg,
                           zBkgErr,
                           PdBad,
                           ModelLevels,
                           zFlags,
                           zPGE,
                           zPGEBd);
  }
}  // namespace ufo
