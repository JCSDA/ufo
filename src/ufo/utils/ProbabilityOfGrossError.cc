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
                         const bool PerformSDiffCheck,
                         std::vector<int> &flags,
                         std::vector<float> &PGE,
                         float ErrVarMax,
                         const std::vector<float> *obsVal2,
                         const std::vector<float> *bkgVal2,
                         std::vector<float> *TotalPd)
  {
    const float missingValueFloat = util::missingValue<float>();
    // PGE multiplication factor used to store PGE values for later use.
    const double PGEMult = 1000.0;
    // Missing data indicator for stored PGEs.
    const double PGEMDI = 1.111;
    // PGEMDI * PGEMult to avoid truncation errors
    const double PGEMDIMult = 1111.0;
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
    // Number of levels in profile, or total number of single-level obs.
    const size_t numLocs = obsVal.size();

    // Combined (obs and bkg) error variance.
    float ErrVar = 0.0;
    // Squared difference from background / ErrVar.
    double SDiff = 0.0;
    // Probability density of good observation.
    double PdGood = 0.0;

    const bool obsErrEmpty = obsErr.empty();
    const bool bkgErrEmpty = bkgErr.empty();

    for (size_t jloc = 0; jloc < numLocs; ++jloc) {
      // Calculate combined error variance.
      if (!obsErrEmpty && !bkgErrEmpty &&
          obsErr[jloc] >= 0 && bkgErr[jloc] >= 0) {
        ErrVar = std::pow(ObErrMult * obsErr[jloc], 2) +
          std::pow(BkgErrMult * bkgErr[jloc], 2);
      } else {
        ErrVar = missingValueFloat;
      }
      // Set combined error variance to maximum value (if defined).
      if (ErrVarMax > 0.0) {
        ErrVar = std::min(ErrVar, ErrVarMax);
      }

      // Update PGE (if all of the required values are present).
      if (obsVal[jloc] != missingValueFloat &&
          bkgVal[jloc] != missingValueFloat &&
          ErrVar != missingValueFloat) {
        if (obsVal2 && bkgVal2 &&
            (*obsVal2)[jloc] != missingValueFloat &&
            (*bkgVal2)[jloc] != missingValueFloat) {  // Vector observable.
          SDiff = (std::pow(obsVal[jloc] - bkgVal[jloc], 2) +
                   std::pow((*obsVal2)[jloc] - (*bkgVal2)[jloc], 2)) / static_cast<double>(ErrVar);
          // Bivariate normal distribution; square root does not appear in denominator.
          PdGood = std::exp(-0.5 * std::min(SDiff, 2.0 * ExpArgMax)) / (2.0 * M_PI * ErrVar);
        } else {  // Scalar observable.
          SDiff = std::pow(obsVal[jloc] - bkgVal[jloc], 2) / ErrVar;
          // Univariate normal distribution; square root appears in denominator.
          PdGood = std::exp(-0.5 * std::min(SDiff, 2.0 * ExpArgMax)) /
            std::sqrt(2.0 * M_PI * ErrVar);
        }

        // Calculate PGE after background check, normalising appropriately.
        const double TotalProbDist = PdBad[jloc] * PGE[jloc] + PdGood * (1.0 - PGE[jloc]);
        PGE[jloc] = (PdBad[jloc] * PGE[jloc]) / TotalProbDist;
        if (TotalPd)
            (*TotalPd)[jloc] = TotalProbDist;

        // Set QC flags.
        flags[jloc] |= ufo::MetOfficeQCFlags::Elem::BackPerfFlag;
        if (PGE[jloc] >= PGECrit) {
          flags[jloc] |= ufo::MetOfficeQCFlags::Elem::BackRejectFlag;
        }
      } else {
        // Deal with missing data.
        if (TotalPd) {
          (*TotalPd)[jloc] = PdBad[jloc];
        }
        SDiff = SDiffCrit;
        PGE[jloc] = PGEMDI;
      }

      // Apply squared difference check if required:
      // reject (o - b)**2/errvar >= SDiffCrit
      if (PerformSDiffCheck &&
          (SDiff >= SDiffCrit ||
           flags[jloc] & ufo::MetOfficeQCFlags::Elem::PermRejectFlag ||
           flags[jloc] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag)) {
        PGE[jloc] = PGEMDI;  // Do not apply buddy check in this case.
        flags[jloc] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
      }
    }
  }
}  // namespace ufo
