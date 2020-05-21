/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckHydrostatic.h"

namespace ufo {
  ProfileCheckHydrostatic::ProfileCheckHydrostatic(const ProfileConsistencyCheckParameters &options,
                                                   const ProfileIndices &profileIndices,
                                                   const ProfileData &profileData,
                                                   ProfileFlags &profileFlags,
                                                   ProfileCheckValidator &profileCheckValidator)
    : ProfileCheckBase(options, profileIndices, profileData, profileFlags, profileCheckValidator),
    ProfileStandardLevels(options)
  {}

  void ProfileCheckHydrostatic::runCheck()
  {
    oops::Log::debug() << " Hydrostatic check" << std::endl;

    const int numLevelsToCheck = profileIndices_.getNumLevelsToCheck();
    const std::vector <float> &pressures = profileData_.getPressures();
    const std::vector <float> &tObs = profileData_.gettObs();
    const std::vector <float> &tBkg = profileData_.gettBkg();
    const std::vector <float> &zObs = profileData_.getzObs();
    const std::vector <float> &zBkg = profileData_.getzBkg();
    std::vector <int> &tFlags = profileFlags_.gettFlags();
    std::vector <int> &zFlags = profileFlags_.getzFlags();
    const std::vector <float> &tObsCorrection = profileFlags_.gettObsCorrection();
    std::vector <float> &zObsCorrection =
      profileFlags_.getzObsCorrection();  // Potentially modified here
    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);

    calcStdLevels(numLevelsToCheck, pressures, tObsFinal, tFlags);
    findHCheckStdLevs();

    HydDesc_ = options_.HydDesc.value();
    DC_.assign(numLevelsToCheck, missingValueFloat);
    ETol_.assign(numLevelsToCheck, missingValueFloat);
    D_.assign(numLevelsToCheck, missingValueFloat);
    E_.assign(numLevelsToCheck + 1, missingValueFloat);
    HydError_.assign(numLevelsToCheck, 0);

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
        profileFlags_.incrementCounterCumul("Num925Miss");
      } else if (IndStd_[jlevstd] - IndStd_[jlevstd - 1] != 1) {
        if (IndStd_[jlevstd - 1] < Ind925_ &&
            IndStd_[jlevstd] > Ind925_) {  // Allow for bigger gaps than two standard levels
          profileFlags_.incrementCounterCumul("Num925Miss");
        } else if (IndStd_[jlevstd - 1] < Ind100_ &&
                   IndStd_[jlevstd] > Ind100_) {  // Missed 100 hPa
          profileFlags_.incrementCounterCumul("Num100Miss");
        } else {
          profileFlags_.incrementCounterCumul("NumStdMiss");  // Missed any other standard level
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
        profileFlags_.incrementCounter("NumAnyErrors");
        HydError_[jlevstd] = 3;  // T or Z error
      } else {
        HydError_[jlevstd] = 0;  // Probably OK
        if (std::fabs(E_[jlevstd]) <= options_.HCheck_EThresh.value() &&
            std::fabs(E_[jlevstd - 1]) <= options_.HCheck_EThreshB.value() &&
            tFlags[jlevB] & ufo::FlagsProfile::InterpolationFlag) {  // use 0.5*ETol(jlevstd) ?
          tFlags[jlevB] &= ~ufo::FlagsProfile::InterpolationFlag;
          oops::Log::debug() << " -> removed interpolation flag on level " << jlevB << std::endl;
        }
      }
    }

    // Hydrostatic decision making algorithm
    if (NumErrors > 0) {
      profileFlags_.incrementCounterCumul("NumHydErrObs");

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
              zFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;
              tFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;
              zFlags[jlev]  |= ufo::FlagsProfile::HydrostaticFlag;
              tFlags[jlev]  |= ufo::FlagsProfile::HydrostaticFlag;
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
            zFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;
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
              zFlags[jlevB] |= ufo::FlagsElem::DataCorrectFlag;
              if (options_.HCheck_CorrectZ.value()) {
                zObsCorrection[jlevB] = CorrApp;
                oops::Log::debug() << " -> Uncorrected zObs: " << zObs[jlevB] << "m" << std::endl;
                oops::Log::debug() << "    zObs correction: " << CorrApp << "m" << std::endl;
                oops::Log::debug() << "    Corrected zObs: "
                                   << zObs[jlevB] + zObsCorrection[jlevB] << "m" << std::endl;
              } else {
                // Observation is rejected
                zFlags[jlevB] |= ufo::FlagsElem::FinalRejectFlag;
              }
            }
            // Height errors in two adjacent levels
          } else if (ENext != missingValueFloat &&
                     std::fabs(E_[jlevstd - 1] + E_[jlevstd] + ENext) <=
                     options_.HCheck_ESumNextThresh.value() &&
                     MinAbsE >= options_.HCheck_MinAbsEThresh.value()) {
            zFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;
            zFlags[jlev]  |= ufo::FlagsProfile::HydrostaticFlag;
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
            tFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;
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

            if (tFlags[jlevB] & ufo::FlagsProfile::InterpolationFlag) {
              int SigB = SigBelow_[jlevstd - 1];
              int SigA = SigAbove_[jlevstd - 1];

              tFlags[SigB] &= ~ufo::FlagsProfile::InterpolationFlag;
              tFlags[SigA] &= ~ufo::FlagsProfile::InterpolationFlag;

              profileFlags_.incrementCounterCumul("NumIntHydErrors");
              oops::Log::debug() << " -> Hyd: remove interpolation flags on levels "
                                 << SigB << " " << SigA << std::endl;
            }

            // Bottom level error in T or Z, usually jlevstd = 3
          } else if (HydError_[jlevstd - 1] == 3 && HydError_[jlevstd] == 0) {
            if (E_[jlevstd - 2] == missingValueFloat) {
              int L1 = StdLev_[jlevstd - 2];

              zFlags[L1] |= ufo::FlagsProfile::HydrostaticFlag;
              tFlags[L1] |= ufo::FlagsProfile::HydrostaticFlag;

              oops::Log::debug() << " -> Failed hydrostatic check "
                                 << "(bottom level error in T or Z) on level " << L1 << std::endl;

              HydError_[jlevstd - 2] = 4;
              HydError_[jlevstd - 1] = 0;

              if (tFlags[L1] & ufo::FlagsProfile::SurfaceLevelFlag) {
                oops::Log::debug() << " -> Baseline error for level " << L1
                                   << "? P = " << pressures[L1] * 0.01 << "hPa, zObs = "
                                   << zObs[L1] << "m, zBkg = " << zBkg[L1]
                                   << ", zObs + E = "
                                   << zObs[L1] + E_[jlevstd - 1] << "m" << std::endl;
              }
            } else {
              // Error in all subsequent heights?
              HydError_[jlevstd - 1] = 6;
              zFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;

              oops::Log::debug() << " -> Failed hydrostatic check "
                                 << "(error in all subsequent heights) on level "
                                 << jlevB << std::endl;
            }
          } else if (HydError_[jlevstd - 1] == 3) {  // T and/or Z error
            zFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;
            tFlags[jlevB] |= ufo::FlagsProfile::HydrostaticFlag;

            oops::Log::debug() << " -> Failed hydrostatic check "
                               << "(T and/or Z error) on level " << jlevB << std::endl;
          }

          // Top level error in T or Z
          if (HydError_[jlevstd] == 3 && E_[jlevstd + 1] == missingValueFloat) {
            zFlags[jlev] |= ufo::FlagsProfile::HydrostaticFlag;
            tFlags[jlev] |= ufo::FlagsProfile::HydrostaticFlag;
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
    profileCheckValidator_.settFlags(profileFlags_.gettFlags());
    profileCheckValidator_.setzFlags(profileFlags_.getzFlags());
    profileCheckValidator_.setNumAnyErrors(profileFlags_.getCounter("NumAnyErrors"));
    profileCheckValidator_.setNum925Miss(profileFlags_.getCounter("Num925Miss"));
    profileCheckValidator_.setNum100Miss(profileFlags_.getCounter("Num100Miss"));
    profileCheckValidator_.setNumStdMiss(profileFlags_.getCounter("NumStdMiss"));
    profileCheckValidator_.setNumHydErrObs(profileFlags_.getCounter("NumHydErrObs"));
    profileCheckValidator_.setNumIntHydErrors(profileFlags_.getCounter("NumIntHydErrors"));
    profileCheckValidator_.setDC(DC_);
    profileCheckValidator_.setETol(ETol_);
    profileCheckValidator_.setD(D_);
    profileCheckValidator_.setE(E_);
    profileCheckValidator_.setHydError(HydError_);
  }
}  // namespace ufo

