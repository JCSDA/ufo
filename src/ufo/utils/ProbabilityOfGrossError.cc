/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/utils/ProbabilityOfGrossError.h"

namespace ufo {
  void BayesianPGEUpdate(const ProbabilityOfGrossErrorParameters &options,
                         const std::vector<float> &obsVal,
                         const std::vector<float> &obsErr,
                         const std::vector<float> &bkgVal,
                         const std::vector<float> &bkgErr,
                         const std::vector<float> &PdBad,
                         const bool ModelLevels,
                         std::vector<int> &flags,
                         std::vector<float> &PGE,
                         std::vector<float> &PGEBd,
                         float ErrVarMax,
                         const std::vector<float> *obsVal2,
                         const std::vector<float> *bkgVal2)
  {
    const float missingValueFloat = util::missingValue(1.0f);
    // PGE multiplication factor used to store PGE values for later use.
    const double PGEMult = 1000.0;
    // Missing data indicator for stored PGEs.
    const double PGEMDI = 1.111;
    // Maximum value of exponent in background QC.
    const double ExpArgMax = options.PGE_ExpArgMax.value();
    // PGE rejection limit.
    const double PGECrit = options.PGE_PGECrit.value();
    // Multiplication factor for observation errors.
    const float ObErrMult = options.PGE_ObErrMult.value();
    // Multiplication factor for background errors.
    const float BkgErrMult = options.PGE_BkgErrMult.value();
    // Critical value for squared difference from background / ErrVar.
    const double SDiffCrit = obsVal2 && bkgVal2 ?
      options.PGE_SDiffCrit.value() * 2.0 :
      options.PGE_SDiffCrit.value();
    // Number of levels in profile.
    const size_t numProfileLevels = obsVal.size();

    // Initialise buddy check PGE to missing data indicator.
    PGEBd.assign(obsVal.size(), PGEMDI);

    // Combined (obs and bkg) error variance.
    float ErrVar = 0.0;
    // Squared difference from background / ErrVar.
    double SDiff = 0.0;
    // Probability density of good observation.
    double PdGood = 0.0;
    // PGE after background check.
    double PGEBk = 0.0;

    const bool obsErrEmpty = obsErr.empty();
    const bool bkgErrEmpty = bkgErr.empty();

    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
      // Calculate combined error variance.
      if (!obsErrEmpty && !bkgErrEmpty &&
          obsErr[jlev] >= 0 && bkgErr[jlev] >= 0) {
        ErrVar = std::pow(ObErrMult * obsErr[jlev], 2) +
          std::pow(BkgErrMult * bkgErr[jlev], 2);
      } else {
        ErrVar = missingValueFloat;
      }
      // Set combined error variance to maximum value (if defined).
      if (ErrVarMax > 0.0) {
        ErrVar = std::min(ErrVar, ErrVarMax);
      }

      // Update PGE (if all of the required values are present).
      if (obsVal[jlev] != missingValueFloat &&
          bkgVal[jlev] != missingValueFloat &&
          ErrVar != missingValueFloat) {
        if (obsVal2 && bkgVal2 &&
            (*obsVal2)[jlev] != missingValueFloat &&
            (*bkgVal2)[jlev] != missingValueFloat) {  // Vector observable.
          SDiff = (std::pow(obsVal[jlev] - bkgVal[jlev], 2) +
                   std::pow((*obsVal2)[jlev] - (*bkgVal2)[jlev], 2)) / static_cast<double>(ErrVar);
          // Bivariate normal distribution; square root does not appear in denominator.
          PdGood = std::exp(-0.5 * std::min(SDiff, 2.0 * ExpArgMax)) / (2.0 * M_PI * ErrVar);
        } else {  // Scalar observable.
          SDiff = std::pow(obsVal[jlev] - bkgVal[jlev], 2) / ErrVar;
          // Univariate normal distribution; square root appears in denominator.
          PdGood = std::exp(-0.5 * std::min(SDiff, 2.0 * ExpArgMax)) /
            std::sqrt(2.0 * M_PI * ErrVar);
        }

        // Calculate PGE after background check, normalising appropriately.
        PGEBk = (PdBad[jlev] * PGE[jlev]) /
          (PdBad[jlev] * PGE[jlev] + PdGood * (1.0 - PGE[jlev]));

        // Set QC flags.
        flags[jlev] |= ufo::MetOfficeQCFlags::Elem::BackPerfFlag;
        if (PGEBk >= PGECrit) {
          flags[jlev] |= ufo::MetOfficeQCFlags::Elem::BackRejectFlag;
        }
      } else {
        // Deal with missing data.
        SDiff = SDiffCrit;
        PGEBk = PGEMDI;
      }

      // Pack PGEs for use in later routines.
      PGE[jlev] = trunc(PGEBk * PGEMult) + PGE[jlev];

      // Model-level data may have additional processing.
      if (ModelLevels &&
          (SDiff >= SDiffCrit ||
           flags[jlev] & ufo::MetOfficeQCFlags::Elem::PermRejectFlag ||
           flags[jlev] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag)) {
        PGEBk = PGEMDI;  // Do not apply buddy check in this case.
        flags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
      }

      // Update PGE for buddy check.
      PGEBd[jlev] = PGEBk;
    }
  }
}  // namespace ufo
