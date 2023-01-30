/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBackgroundRelativeHumidity.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckBackgroundRelativeHumidity>
  makerProfileCheckBackgroundRelativeHumidity_("BackgroundRelativeHumidity");

  ProfileCheckBackgroundRelativeHumidity::ProfileCheckBackgroundRelativeHumidity
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileCheckBackgroundRelativeHumidity::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Background check for relative humidity" << std::endl;

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &rhObs =
      profileDataHandler.get<float>(ufo::VariableNames::obs_relative_humidity);
    const std::vector <float> &rhObsErr =
      profileDataHandler.get<float>(ufo::VariableNames::obserr_relative_humidity);
    const std::vector <float> &rhBkg =
      profileDataHandler.get<float>(ufo::VariableNames::hofx_relative_humidity);
    const std::vector <float> &rhBkgErr =
      profileDataHandler.get<float>
      (options_.bkgErrGroup.value() + "/" + options_.bkgErrName_relative_humidity.value());
    std::vector <float> &rhPGE =
      profileDataHandler.get<float>(ufo::VariableNames::pge_relative_humidity);
    std::vector <int> &rhFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    const std::vector <int> &extended_obs_space =
      profileDataHandler.get<int>(ufo::VariableNames::extended_obs_space);
    const bool ModelLevels = std::find(extended_obs_space.begin(), extended_obs_space.end(), 1)
      != extended_obs_space.end();

    if (!oops::allVectorsSameNonZeroSize(rhObs, rhObsErr, rhBkg, rhBkgErr,
                                         rhPGE, rhFlags)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(rhObs, rhObsErr, rhBkg, rhBkgErr,
                                                    rhPGE, rhFlags)
                         << std::endl;
      return;
    }

    // Probability density of 'bad' observations.
    std::vector <float> PdBad(numProfileLevels, options_.BkCheck_PdBad_rh.value());
    // Local version of relative humidity background error.
    std::vector <float> BackgrErrRH(numProfileLevels, 0.0);
    // Local version of relative humidity observation error.
    std::vector <float> ObErrRH(numProfileLevels, 0.0);

    // Relax QC to take account of long-tailed error distributions.
    const float sqrt2 = std::sqrt(2.0);
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      BackgrErrRH[jlev] = missingValueFloat;
      ObErrRH[jlev] = missingValueFloat;
      if (rhBkgErr[jlev] != missingValueFloat)
        BackgrErrRH[jlev] = sqrt2 * rhBkgErr[jlev];
      if (rhObsErr[jlev] != missingValueFloat)
        ObErrRH[jlev] = sqrt2 * rhObsErr[jlev];
    }

    // Calculate probability of gross error.
    ufo::BayesianPGEUpdate(options_.PGEParameters,
                           rhObs,
                           ObErrRH,  // Used instead of rhObsErr.
                           rhBkg,
                           BackgrErrRH,  // Used instead of rhBkgErr.
                           PdBad,
                           ModelLevels,
                           rhFlags,
                           rhPGE,
                           options_.BkCheck_ErrVarMax_rh.value());
  }
}  // namespace ufo
