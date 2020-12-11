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
  (const ProfileConsistencyCheckParameters &options,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckBackgroundRelativeHumidity::runCheck()
  {
    oops::Log::debug() << " Background check for relative humidity" << std::endl;

    const size_t numProfileLevels = profileDataHandler_.getNumProfileLevels();
    const bool ModelLevels = options_.modellevels.value();
    const std::vector <float> &rhObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_relative_humidity);
    const std::vector <float> &rhObsErr =
       profileDataHandler_.get<float>(ufo::VariableNames::obserr_relative_humidity);
    const std::vector <float> &rhBkg =
      profileDataHandler_.get<float>(ufo::VariableNames::hofx_relative_humidity);
    const std::vector <float> &rhBkgErr =
      profileDataHandler_.get<float>(ufo::VariableNames::bkgerr_relative_humidity);
    std::vector <float> &rhPGE =
      profileDataHandler_.get<float>(ufo::VariableNames::pge_relative_humidity);
    std::vector <float> &rhPGEBd =
      profileDataHandler_.get<float>(ufo::VariableNames::pgebd_relative_humidity);
    std::vector <int> &rhFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    const std::vector <int> &timeFlags =
      profileDataHandler_.get<int>(ufo::VariableNames::qcflags_time);

    if (!oops::allVectorsSameNonZeroSize(rhObs, rhObsErr, rhBkg, rhBkgErr,
                                         rhPGE, rhFlags, timeFlags)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(rhObs, rhObsErr, rhBkg, rhBkgErr,
                                                      rhPGE, rhFlags, timeFlags)
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
      if (timeFlags[jlev])
        rhFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
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
                           rhPGEBd,
                           options_.BkCheck_ErrVarMax_rh.value());
  }
}  // namespace ufo
