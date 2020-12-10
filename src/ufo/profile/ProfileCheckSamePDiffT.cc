/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckSamePDiffT.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckSamePDiffT> makerProfileCheckSamePDiffT_("SamePDiffT");

  ProfileCheckSamePDiffT::ProfileCheckSamePDiffT(const ProfileConsistencyCheckParameters &options,
                                                 ProfileDataHandler &profileDataHandler,
                                                 ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator)
  {}

  void ProfileCheckSamePDiffT::runCheck()
  {
    oops::Log::debug() << " Test for same pressure and different temperature" << std::endl;
    int jlevprev = -1;
    int NumErrors = 0;

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
    std::vector <int> &NumSamePErrObs =
      profileDataHandler_.get<int>(ufo::VariableNames::counter_NumSamePErrObs);
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

    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (tObs[jlev] == missingValueFloat) continue;

      if (jlevprev == -1) {
        jlevprev = jlev;
        continue;
      }

      if (pressures[jlev] == pressures[jlevprev]) {
        int jlevuse = jlevprev;
        if (std::abs(tObsFinal[jlev] - tObsFinal[jlevprev]) > options_.SPDTCheck_TThresh.value()) {
          NumErrors++;
          NumAnyErrors[0]++;

          // Choose which level to flag
          if (std::abs(tObsFinal[jlev] - tBkg[jlev]) <=
              std::abs(tObsFinal[jlevprev] - tBkg[jlevprev])) {
            tFlags[jlevprev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
            tFlags[jlev]     |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
            jlevuse = jlev;
          } else {
            tFlags[jlevprev] |= ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
            tFlags[jlev]     |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
          }

          oops::Log::debug() << " -> Failed same P/different T check for levels "
                             << jlevprev << " and " << jlev << std::endl;
          oops::Log::debug() << " -> Level " << jlevprev << ": "
                             << "P = " << pressures[jlevprev] * 0.01 << "hPa, tObs = "
                             << tObsFinal[jlevprev] - ufo::Constants::t0c << "C, "
                             << "tBkg = " << tBkg[jlevprev] - ufo::Constants::t0c
                             << "C" << std::endl;
          oops::Log::debug() << " -> Level " << jlev << ": "
                             << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                             << tObsFinal[jlev] - ufo::Constants::t0c << "C, "
                             << "tBkg = " << tBkg[jlev] - ufo::Constants::t0c
                             << "C" << std::endl;
          oops::Log::debug() << " -> tObs difference: " << tObsFinal[jlev] - tObsFinal[jlevprev]
                             << std::endl;
          oops::Log::debug() << " -> Use level " << jlevuse << std::endl;
        }
        jlevprev = jlevuse;
      } else {  // Distinct pressures
        jlevprev = jlev;
      }
    }
    if (NumErrors > 0) NumSamePErrObs[0]++;
  }
}  // namespace ufo

