/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGEPRESSURE_H_
#define UFO_PROFILE_PROFILEAVERAGEPRESSURE_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "ufo/profile/ProfileCheckBase.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
  class ProfileDataHandler;
  class ProfileDataHolder;
}

namespace ufo {

  /// \brief Profile QC: apply various transformations to observed and model pressures.
  /// The transformed pressures are used in subsequent profile averaging routines.
  class ProfileAveragePressure : public ProfileCheckBase {
   public:
    explicit ProfileAveragePressure(const ConventionalProfileProcessingParameters &options);

    /// Run check on all profiles.
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Fill variables in validator.
    void fillValidationData(ProfileDataHolder &profileDataHolder,
                            bool extended_obs_space);

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({oops::Variable{ufo::VariableNames::geovals_pressure},
            oops::Variable{ufo::VariableNames::geovals_pressure_rho_minus_one}});}

    /// List of names of GeoVaLs used in check validation.
    oops::Variables getValidationGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_testreference_logP,
            ufo::VariableNames::geovals_testreference_ExnerP,
            ufo::VariableNames::geovals_testreference_logP_rho,
            ufo::VariableNames::geovals_testreference_ExnerP_rho});}

   private:  // functions
    /// Calculate log(pressure).
    void logPressure(const std::vector <float> &pressures,
                     std::vector <float> &logP);

    /// Calculate Exner pressure.
    void ExnerPressure(const std::vector <float> &pressures,
                       std::vector <float> &ExnerP);

    /// Calculate big gap for each pressure.
    void bigPressureGaps(const std::vector <float> &pressures,
                         const int ObsType,
                         std::vector <float> &bigPgaps);

    /// Run check on a profile in the original ObsSpace and
    /// put the averaged data into the corresponding profile in the extended ObsSpace.
    void runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                            ProfileDataHolder &profileExtended);
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGEPRESSURE_H_
