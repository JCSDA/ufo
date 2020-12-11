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
   ProfileDataHandler &profileDataHandler,
   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckUnstableLayer::runCheck()
  {
    oops::Log::debug() << " Unstable layer/superadiabat check" << std::endl;

    const int numProfileLevels = profileDataHandler_.getNumProfileLevels();

    const std::vector <float> &pressures =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::hofx_air_temperature);
    std::vector <int> &tFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &NumAnyErrors =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_NumAnyErrors);
    std::vector <int> &NumSuperadiabat =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_NumSuperadiabat);
    const std::vector <float> &tObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::obscorrection_air_temperature);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg, tFlags, tObsCorrection)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(pressures, tObs, tBkg, tFlags, tObsCorrection)
                           << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    PBottom_ = 0.0;

    int jlevprev = 0;
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      // Ignore this level if it has been flagged as rejected.
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag) continue;
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
            tFlags[jlevprev] |= ufo::MetOfficeQCFlags::Profile::SuperadiabatFlag;
            tFlags[jlev]     |= ufo::MetOfficeQCFlags::Profile::SuperadiabatFlag;

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
    std::vector <float> PBottom(profileDataHandler_.getNumProfileLevels(), PBottom_);
    profileDataHandler_.set(ufo::VariableNames::PBottom, std::move(PBottom));
  }
}  // namespace ufo
