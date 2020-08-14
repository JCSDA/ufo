/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckUnstableLayer.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckUnstableLayer>
  makerProfileCheckUnstableLayer_("UnstableLayer");

  ProfileCheckUnstableLayer::ProfileCheckUnstableLayer
  (const ProfileConsistencyCheckParameters &options,
   const ProfileIndices &profileIndices,
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckUnstableLayer::runCheck()
  {
    oops::Log::debug() << " Unstable layer/superadiabat check" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();

    const std::vector <float> &pressures =
       profileDataHandler_.get<float>(ufo::VariableNames::name_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler_.get<float>(ufo::VariableNames::name_obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::name_hofx_air_temperature);
    std::vector <int> &tFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::name_qc_tFlags);
    std::vector <int> &NumAnyErrors =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_NumAnyErrors);
    std::vector <int> &NumSuperadiabat =
       profileDataHandler_.get<int>(ufo::VariableNames::name_counter_NumSuperadiabat);
    const std::vector <float> &tObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::name_tObsCorrection);

    if (oops::anyVectorEmpty(pressures, tObs, tBkg, tFlags, tObsCorrection)) {
      oops::Log::warning() << "At least one vector is empty. "
                           << "Check will not be performed." << std::endl;
      return;
    }
    if (!oops::allVectorsSameSize(pressures, tObs, tBkg, tFlags, tObsCorrection)) {
      oops::Log::warning() << "Not all vectors have the same size. "
                           << "Check will not be performed." << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    PBottom_ = 0.0;

    int jlevprev = 0;
    for (int jlev = 0; jlev < numLevelsToCheck; ++jlev) {
      if (tFlags[jlev] & ufo::FlagsElem::FinalRejectFlag) continue;  // Ignore this level
      if (tObsFinal[jlev] != missingValueFloat &&
          pressures[jlev] > options_.ULCheck_MinP.value()) {
        if (PBottom_ == 0.0) {
          PBottom_ = pressures[jlev];
        } else {
          //  Temperature calculated adiabatically from prev level
          const float Tadiabat = tObsFinal[jlevprev] *
            std::pow(pressures[jlev] / pressures[jlevprev], ufo::Constants::rd_over_cp);
          if (tObsFinal[jlev] - Tadiabat <= options_.ULCheck_SuperadiabatTol.value() &&
              pressures[jlevprev] <= PBottom_ - options_.ULCheck_PBThresh.value()) {
            NumAnyErrors[0]++;
            NumSuperadiabat[0]++;
            tFlags[jlevprev] |= ufo::FlagsProfile::SuperadiabatFlag;
            tFlags[jlev]     |= ufo::FlagsProfile::SuperadiabatFlag;

            oops::Log::debug() << " -> Failed unstable layer/superadiabat check for levels "
                               << jlevprev << " and " << jlev << std::endl;
            oops::Log::debug() << " -> Tadiabat = " << Tadiabat - ufo::Constants::t0c << "C"
                               << std::endl;
            oops::Log::debug() << " -> Level " << jlevprev << ": "
                               << "P = " << pressures[jlevprev] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlevprev] - ufo::Constants::t0c << "C, tBkg = "
                               << tBkg[jlevprev] - ufo::Constants::t0c << "C, "
                               << "tObs - Tadiabat = " << tObsFinal[jlev] - Tadiabat
                               << "C" << std::endl;
            oops::Log::debug() << " -> Level " << jlev << ": "
                               << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlev] - ufo::Constants::t0c << "C, tBkg = "
                               << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;
          }
        }
        jlevprev = jlev;
      }
    }
  }

  void ProfileCheckUnstableLayer::fillValidator()
  {
    std::vector <float> PBottom(profileIndices_.getNumLevelsToCheck(), PBottom_);
    profileDataHandler_.set(ufo::VariableNames::name_PBottom, std::move(PBottom));
  }
}  // namespace ufo
