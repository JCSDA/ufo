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
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Profile QC: apply various transformations to observed and model pressures.
  /// The transformed pressures are used in subsequent profile averaging routines.
  class ProfileAveragePressure : public ProfileCheckBase {
   public:
    explicit ProfileAveragePressure(const ProfileConsistencyCheckParameters &options);

    /// Run check.
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Fill variables in validator.
    void fillValidationData(ProfileDataHandler &profileDataHandler) override;

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_pressure,
            ufo::VariableNames::geovals_pressure_rho});}

    /// List of names of GeoVaLs used in check validation.
    oops::Variables getValidationGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_logP,
            ufo::VariableNames::geovals_ExnerP,
            ufo::VariableNames::geovals_logP_rho,
            ufo::VariableNames::geovals_ExnerP_rho});}

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
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGEPRESSURE_H_
