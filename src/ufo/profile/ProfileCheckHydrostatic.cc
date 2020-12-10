/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckHydrostatic.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileCheckHydrostatic> makerProfileCheckHydrostatic_("Hydrostatic");

  ProfileCheckHydrostatic::ProfileCheckHydrostatic(const ProfileConsistencyCheckParameters &options,
                                                   ProfileDataHandler &profileDataHandler,
                                                   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileDataHandler, profileCheckValidator),
    ProfileStandardLevels(options)
  {}

  void ProfileCheckHydrostatic::runCheck()
  {
    oops::Log::debug() << " Hydrostatic check" << std::endl;

    const int numProfileLevels = profileDataHandler_.getNumProfileLevels();

    const std::vector <float> &pressures =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <float> &tObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_air_temperature);
    const std::vector <float> &tBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::hofx_air_temperature);
    const std::vector <float> &zObs =
       profileDataHandler_.get<float>(ufo::VariableNames::obs_geopotential_height);
    const std::vector <float> &zBkg =
       profileDataHandler_.get<float>(ufo::VariableNames::hofx_geopotential_height);
    std::vector <int> &tFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &zFlags =
       profileDataHandler_.get<int>(ufo::VariableNames::qcflags_geopotential_height);
    std::vector <int> &NumAnyErrors =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_NumAnyErrors);
    std::vector <int> &Num925Miss =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_Num925Miss);
    std::vector <int> &Num100Miss =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_Num100Miss);
    std::vector <int> &NumStdMiss =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_NumStdMiss);
    std::vector <int> &NumHydErrObs =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_NumHydErrObs);
    std::vector <int> &NumIntHydErrors =
       profileDataHandler_.get<int>(ufo::VariableNames::counter_NumIntHydErrors);
    const std::vector <float> &tObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::obscorrection_air_temperature);
    std::vector <float> &zObsCorrection =
       profileDataHandler_.get<float>(ufo::VariableNames::obscorrection_geopotential_height);

    if (!oops::allVectorsSameNonZeroSize(pressures, tObs, tBkg, zObs, zBkg, tFlags, zFlags,
                                         tObsCorrection, zObsCorrection)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Check will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(pressures, tObs, tBkg, zObs, zBkg, tFlags,
                                                      zFlags, tObsCorrection, zObsCorrection)
                           << std::endl;
      return;
    }

    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    calcStdLevels(numProfileLevels, pressures, tObsFinal, tFlags);
    findHCheckStdLevs();

    HydDesc_ = options_.HydDesc.value();
    DC_.assign(numProfileLevels, missingValueFloat);
    ETol_.assign(numProfileLevels, missingValueFloat);
    D_.assign(numProfileLevels, missingValueFloat);
    E_.assign(numProfileLevels + 1, missingValueFloat);
    HydError_.assign(numProfileLevels, 0);

    int NumErrors = 0;
    // Find large thickness residuals
    for (int jlevstd = 1; jlevstd < NumStd_; ++jlevstd) {
      int jlev = StdLev_[jlevstd];  // Standard level
      int jlevB = StdLev_[jlevstd - 1];  // Standard level below this one
      if (zObs[jlev] == missingValueFloat ||
          zObs[jlevB] == missingValueFloat) continue;
      if (IndStd_[jlevstd - 1] == -1) {  // Surface
        if (std::fabs(pressures[jlevB] - pressures[jlev]) >
            options_.HCheck_SurfacePThresh.value()) continue;
      } else if (IndStd_[jlevstd] == Ind925_ + 1 &&
                 IndStd_[jlevstd - 1] == Ind925_ - 1) {  // Missed 925 hPa
        Num925Miss[0]++;
      } else if (IndStd_[jlevstd] - IndStd_[jlevstd - 1] != 1) {
        if (IndStd_[jlevstd - 1] < Ind925_ &&
            IndStd_[jlevstd] > Ind925_) {  // Allow for bigger gaps than two standard levels
          Num925Miss[0]++;
        } else if (IndStd_[jlevstd - 1] < Ind100_ &&
                   IndStd_[jlevstd] > Ind100_) {  // Missed 100 hPa
          Num100Miss[0]++;
        } else {  // Missed any other standard level
          NumStdMiss[0]++;
        }

        oops::Log::debug() << " Gap in standard levels" << std::endl;
        oops::Log::debug() << " -> Level " << jlev << ": "
                           << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                           << tObsFinal[jlev] - ufo::Constants::t0c << "C, "
                           << "tBkg = " << tBkg[jlev] - ufo::Constants::t0c << "C" << std::endl;
        oops::Log::debug() << " -> Level " << jlevB << ": "
                           << "P = " << pressures[jlevB] * 0.01 << "hPa, tObs = "
                           << tObsFinal[jlevB] - ufo::Constants::t0c << "C, "
                           << "tBkg = " << tBkg[jlevB] - ufo::Constants::t0c
                           << "C" << std::endl;
        oops::Log::debug() << " -> IndStd[" << jlevstd << "] = "
                           << IndStd_[jlevstd] << ", "
                           << "IndStd[" << jlevstd - 1 << "] = "
                           << IndStd_[jlevstd - 1] << std::endl;
        continue;
      }

      DC_[jlevstd] = 0.5 * ufo::Constants::rd_over_g * (LogP_[jlevB] - LogP_[jlev]);
      D_[jlevstd]  = DC_[jlevstd] *
        (tObsFinal[jlevB] + tObsFinal[jlev]);  // Thickness calculated from temperature

      // For neutral stability T(jlev) = T(jlevB) * TRatio
      float TRatio = std::pow(pressures[jlev] / pressures[jlevB], ufo::Constants::rd_over_cp);
      float DB = DC_[jlevstd] * (1.0 + TRatio) *
        tObsFinal[jlevB];  // Min thickness given T(jlevB)
      float DA = DC_[jlevstd] * (1.0 / TRatio + 1.0) *
        tObsFinal[jlev];  // Max thickness given T(jlev)
      ETol_[jlevstd] = options_.HCheck_ETolMult.value() * (DA - DB);
      float ETolMax = options_.HCheck_ETolMax.value();
      if (pressures[jlevB] <= options_.HCheck_ETolMaxPThresh.value())
        ETolMax = options_.HCheck_ETolMaxLarger.value();
      float ETolMin = options_.HCheck_ETolMin.value();  // ETolMin = 20.0 m in GGDPS
      ETol_[jlevstd] = std::max(std::min(ETol_[jlevstd], ETolMax), ETolMin);
      E_[jlevstd] = zObs[jlev] - zObs[jlevB] - D_[jlevstd];
      if (std::fabs(E_[jlevstd]) > ETol_[jlevstd]) {
        NumErrors++;
        NumAnyErrors[0]++;
        HydError_[jlevstd] = 3;  // T or Z error
      } else {
        HydError_[jlevstd] = 0;  // Probably OK
        if (std::fabs(E_[jlevstd]) <= options_.HCheck_EThresh.value() &&
            std::fabs(E_[jlevstd - 1]) <= options_.HCheck_EThreshB.value() &&
            tFlags[jlevB] & ufo::MetOfficeQCFlags::Profile::InterpolationFlag) {
          tFlags[jlevB] &= ~ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
          oops::Log::debug() << " -> removed interpolation flag on level " << jlevB << std::endl;
        }
      }
    }

    // Hydrostatic decision making algorithm
    if (NumErrors > 0) {
      NumHydErrObs[0]++;

      for (int jlevstd = 2; jlevstd < NumStd_; ++jlevstd) {
        // Check for duplicate std levels
        if (eckit::types::is_approximately_equal(DC_[jlevstd - 1], 0.0f) ||
            eckit::types::is_approximately_equal(DC_[jlevstd], 0.0f)) continue;
        int jlev = StdLev_[jlevstd];  // Standard level
        int jlevB = StdLev_[jlevstd - 1];  // Standard level below
        if (HydError_[jlevstd] == 3 || HydError_[jlevstd - 1] == 3) {
          if (E_[jlevstd] == missingValueFloat) {
            // Checked previous time as top level
            continue;
          }
          if (E_[jlevstd - 1] == missingValueFloat) {
            if (E_[jlevstd + 1] == missingValueFloat) {
              zFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
              tFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
              zFlags[jlev]  |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
              tFlags[jlev]  |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
              oops::Log::debug() << " -> Isolated large residual on levels "
                                 << jlev << " and " << jlevB << std::endl;
            }
            continue;
          }
          // Possible temperature corrections

          // These equations are expressed differently to their equivalents in the OPS code.
          // This avoids a vectorisation problem in the clang compiler:
          // previously, when E_ was divided by DC_, some (extraneous) machine registers
          // were filled with zeros prior to the division occurring.
          // A floating point exception was then thrown due to division by zero even though
          // the registers in question played no further part in the calculation.
          float EDC1 = E_[jlevstd - 1] * DC_[jlevstd];
          float EDC2 = E_[jlevstd] * DC_[jlevstd - 1];
          float MinAbsEDC = std::min(std::fabs(EDC1), std::fabs(EDC2));
          float CorrMinThreshDC = options_.HCheck_CorrMinThresh.value()
            * DC_[jlevstd] * DC_[jlevstd - 1];
          float AbsEDCDiff = std::fabs(EDC1 - EDC2);
          float CorrDiffThreshDC = options_.HCheck_CorrDiffThresh.value()
            * DC_[jlevstd] * DC_[jlevstd - 1];

          float MinAbsE = std::min(std::fabs(E_[jlevstd - 1]), std::fabs(E_[jlevstd]));

          float ENext = missingValueFloat;
          if (jlevstd < NumStd_ - 1) ENext = E_[jlevstd + 1];

          //  Height error
          if ((std::fabs(E_[jlevstd - 1] + E_[jlevstd]) <=
               options_.HCheck_ESumThresh.value() &&
               MinAbsE >= options_.HCheck_MinAbsEThresh.value()) ||
              (std::fabs(E_[jlevstd - 1] + E_[jlevstd]) <=
               options_.HCheck_ESumThreshLarger.value() &&
               MinAbsE >= options_.HCheck_MinAbsEThreshLarger.value())) {
            zFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            oops::Log::debug() << " -> Failed hydrostatic check (height error) on level "
                               << jlevB << std::endl;
            HydError_[jlevstd - 1] = 1;
            HydError_[jlevstd] = 0;
            float Corr = 0.5 * (E_[jlevstd] - E_[jlevstd - 1]);  // Average of adjacent levels
            float CorrApp = 100.0 * std::round(Corr / 100.0);  // Round to nearest 100 m
            if (std::fabs(Corr - CorrApp) > options_.HCheck_CorrThresh.value()) CorrApp = 0.0;
            oops::Log::debug() << " -> P = " << pressures[jlevB] * 0.01
                               << "hPa, zObs = " << zObs[jlevB] << "m, "
                               << "Z Correction? " << Corr << "m"
                               << ", rounded = " << CorrApp << "m" << std::endl;
            if (CorrApp != 0.0) {
              zFlags[jlevB] |= ufo::MetOfficeQCFlags::Elem::DataCorrectFlag;
              if (options_.HCheck_CorrectZ.value()) {
                zObsCorrection[jlevB] = CorrApp;
                oops::Log::debug() << " -> Uncorrected zObs: " << zObs[jlevB] << "m" << std::endl;
                oops::Log::debug() << "    zObs correction: " << CorrApp << "m" << std::endl;
                oops::Log::debug() << "    Corrected zObs: "
                                   << zObs[jlevB] + zObsCorrection[jlevB] << "m" << std::endl;
              } else {
                // Observation is rejected
                zFlags[jlevB] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
              }
            }
            // Height errors in two adjacent levels
          } else if (ENext != missingValueFloat &&
                     std::fabs(E_[jlevstd - 1] + E_[jlevstd] + ENext) <=
                     options_.HCheck_ESumNextThresh.value() &&
                     MinAbsE >= options_.HCheck_MinAbsEThresh.value()) {
            zFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            zFlags[jlev]  |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            HydError_[jlevstd - 1] = 1;
            HydError_[jlevstd] = 1;

            oops::Log::debug() << " -> Failed hydrostatic check (height error) on levels "
                               << jlevB << " and " << jlev << std::endl;

            float Corr = -E_[jlevstd - 1];
            oops::Log::debug() << " -> P = " << pressures[jlevB] * 0.01 << "hPa, zObs = "
                               << zObs[jlevB] << "m, "
                               << "Z Correction? " << Corr << "m" << std::endl;
            Corr = ENext;
            oops::Log::debug() << " -> P = " << pressures[jlev] * 0.01 << "hPa, zObs = "
                               << zObs[jlev] << "m, "
                               << "Z Correction? " << Corr << "m" << std::endl;

            // Temperature error
          } else if (MinAbsE >= options_.HCheck_MinAbsEThreshT.value() &&
                     AbsEDCDiff <= CorrDiffThreshDC &&
                     MinAbsEDC >= CorrMinThreshDC) {
            tFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            HydError_[jlevstd - 1] = 2;
            HydError_[jlevstd] = 0;

            oops::Log::debug() << " -> Failed hydrostatic check (temperature error) on level "
                               << jlevB << std::endl;

            // Potential T correction
            float Corr1 = E_.at(jlevstd - 1) / DC_.at(jlevstd - 1);
            float Corr2 = E_.at(jlevstd) / DC_.at(jlevstd);
            float Corr = 0.5 * (Corr1 + Corr2);
            oops::Log::debug() << " -> P = " << pressures[jlevB] * 0.01 << "hPa, tObs = "
                               << tObsFinal[jlevB] - ufo::Constants::t0c << "C, "
                               << "T Correction? " << Corr << "C, "
                               << " Corr1, Corr2 = "
                               << Corr1 << "C, " << Corr2 << "C , DC[" << jlevstd - 1
                               << "], DC[" << jlevstd << "] = "
                               << DC_[jlevstd - 1] << ", " << DC_[jlevstd]
                               << std::endl;

            if (tFlags[jlevB] & ufo::MetOfficeQCFlags::Profile::InterpolationFlag) {
              int SigB = SigBelow_[jlevstd - 1];
              int SigA = SigAbove_[jlevstd - 1];

              tFlags[SigB] &= ~ufo::MetOfficeQCFlags::Profile::InterpolationFlag;
              tFlags[SigA] &= ~ufo::MetOfficeQCFlags::Profile::InterpolationFlag;

              NumIntHydErrors[0]++;
              oops::Log::debug() << " -> Hyd: remove interpolation flags on levels "
                                 << SigB << " " << SigA << std::endl;
            }

            // Bottom level error in T or Z, usually jlevstd = 3
          } else if (HydError_[jlevstd - 1] == 3 && HydError_[jlevstd] == 0) {
            if (E_[jlevstd - 2] == missingValueFloat) {
              int L1 = StdLev_[jlevstd - 2];

              zFlags[L1] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
              tFlags[L1] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;

              oops::Log::debug() << " -> Failed hydrostatic check "
                                 << "(bottom level error in T or Z) on level " << L1 << std::endl;

              HydError_[jlevstd - 2] = 4;
              HydError_[jlevstd - 1] = 0;

              if (tFlags[L1] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag) {
                oops::Log::debug() << " -> Baseline error for level " << L1
                                   << "? P = " << pressures[L1] * 0.01 << "hPa, zObs = "
                                   << zObs[L1] << "m, zBkg = " << zBkg[L1]
                                   << ", zObs + E = "
                                   << zObs[L1] + E_[jlevstd - 1] << "m" << std::endl;
              }
            } else {
              // Error in all subsequent heights?
              HydError_[jlevstd - 1] = 6;
              zFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;

              oops::Log::debug() << " -> Failed hydrostatic check "
                                 << "(error in all subsequent heights) on level "
                                 << jlevB << std::endl;
            }
          } else if (HydError_[jlevstd - 1] == 3) {  // T and/or Z error
            zFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            tFlags[jlevB] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;

            oops::Log::debug() << " -> Failed hydrostatic check "
                               << "(T and/or Z error) on level " << jlevB << std::endl;
          }

          // Top level error in T or Z
          if (HydError_[jlevstd] == 3 && E_[jlevstd + 1] == missingValueFloat) {
            zFlags[jlev] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            tFlags[jlev] |= ufo::MetOfficeQCFlags::Profile::HydrostaticFlag;
            HydError_[jlevstd] = 5;

            oops::Log::debug() << " -> Failed hydrostatic check "
                               << "(top level error in T or Z) on level " << jlev << std::endl;
          }
        }
      }

      for (int jlevstd = 0; jlevstd < NumStd_; ++jlevstd) {
        int HydType = HydError_[jlevstd];
        int jlev = StdLev_[jlevstd];  // Standard level

        oops::Log::debug() << " -> Level " << jlev << ": "
                           << "P = " << pressures[jlev] * 0.01 << "hPa, tObs = "
                           << tObsFinal[jlev] - ufo::Constants::t0c << "C, "
                           << "tBkg = " << tBkg[jlev] - ufo::Constants::t0c << "C, "
                           << "zObs = " << zObs[jlev] << "m, zBkg = " << zBkg[jlev] << "m, "
                           << "D = " << D_[jlevstd] << ", E = " << E_[jlevstd]
                           << ", ETol = " << ETol_[jlevstd] << ", DC = " << DC_[jlevstd]
                           << ", HydDesc = " << HydDesc_[HydType] << " " << std::endl;
      }
    }
  }

  void ProfileCheckHydrostatic::fillValidator()
  {
    profileDataHandler_.set(ufo::VariableNames::DC, std::move(DC_));
    profileDataHandler_.set(ufo::VariableNames::ETol, std::move(ETol_));
    profileDataHandler_.set(ufo::VariableNames::D, std::move(D_));
    profileDataHandler_.set(ufo::VariableNames::E, std::move(E_));
    profileDataHandler_.set(ufo::VariableNames::HydError, std::move(HydError_));
  }
}  // namespace ufo

