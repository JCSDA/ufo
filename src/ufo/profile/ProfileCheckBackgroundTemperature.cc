/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBackgroundTemperature.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundTemperature>
  makerProfileCheckBackgroundTemperature_("BackgroundTemperature");

  ProfileCheckBackgroundTemperature::ProfileCheckBackgroundTemperature
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckBackgroundTemperature::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Background check for temperature" << std::endl;

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &Latitude =
      profileDataHandler.get<float>(ufo::VariableNames::Latitude);
    const std::vector <float> &pressures =
      profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
      profileDataHandler.get<float>(ufo::VariableNames::obs_air_temperature);
    const std::vector <float> &tObsErr =
      profileDataHandler.get<float>(ufo::VariableNames::obserr_air_temperature);
    const std::vector <float> &tBkg =
      profileDataHandler.get<float>(ufo::VariableNames::hofx_air_temperature);
    const std::vector <float> &tBkgErr =
      profileDataHandler.get<float>
      (options_.bkgErrGroup.value() + "/" + options_.bkgErrName_air_temperature.value());
    std::vector <float> &tPGE =
      profileDataHandler.get<float>(ufo::VariableNames::pge_air_temperature);
    std::vector <int> &tFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_air_temperature);
    const std::vector <float> &tObsCorrection =
       profileDataHandler.get<float>(ufo::VariableNames::obscorrection_air_temperature);
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::VariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(Latitude, pressures,
                                         tObs, tObsErr, tBkg, tBkgErr,
                                         tPGE, tFlags, tObsCorrection)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(Latitude, pressures,
                                                    tObs, tObsErr, tBkg, tBkgErr,
                                                    tPGE, tFlags, tObsCorrection)
                         << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, options_.BkCheck_PdBad_t.value());
    // Local version of temperature background error.
    std::vector <float> BackgrErrT(numProfileLevels, 0);
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      BackgrErrT[jlev] = tBkgErr[jlev];
    }
    // Extra representivity error for data on reported levels.
    if (!ModelLevels) {
      const float Psplit =
        std::fabs(Latitude[0]) < options_.BkCheck_Psplit_latitude_tropics ?
                                 options_.BkCheck_Psplit_tropics.value() :
                                 options_.BkCheck_Psplit_extratropics.value();
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (pressures[jlev] <= Psplit) {
          BackgrErrT[jlev] = tBkgErr[jlev] == missingValueFloat ?
            missingValueFloat :
            tBkgErr[jlev] * options_.BkCheck_ErrorInflationBelowPsplit.value();
        } else {
          BackgrErrT[jlev] = tBkgErr[jlev] == missingValueFloat ?
            missingValueFloat :
            tBkgErr[jlev] * options_.BkCheck_ErrorInflationAbovePsplit.value();
        }
      }
    }

    // Modify observation PGE if certain flags have been set.
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SuperadiabatFlag)
        tPGE[jlev] = 0.5 + 0.5 * tPGE[jlev];
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::InterpolationFlag)
        tPGE[jlev] = 0.5 + 0.5 * tPGE[jlev];
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::HydrostaticFlag)
        tPGE[jlev] = 0.5 + 0.5 * tPGE[jlev];
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           tObsFinal,
                           tObsErr,
                           tBkg,
                           BackgrErrT,  // Used instead of tBkgErr.
                           PdBad,
                           ModelLevels,
                           tFlags,
                           tPGE);
  }
}  // namespace ufo
