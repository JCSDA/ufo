/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGETEMPERATURE_H_
#define UFO_PROFILE_PROFILEAVERAGETEMPERATURE_H_

#include <algorithm>
#include <cmath>
#include <vector>

#include "oops/base/Variables.h"

#include "ufo/profile/ProfileAverageUtils.h"
#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class ConventionalProfileProcessingParameters;
}

namespace ufo {

  /// \brief Profile QC: average temperature observations onto model levels.
  ///
  /// Vectors produced by the AveragePressure routine must be present
  /// otherwise the exception eckit::BadValue will be thrown.
  ///
  /// The vertical processing of temperature is based on calculating the thickness
  /// of the model layers (rather than just averaging the temperatures).
  /// The potential temperature in each layer is converted to temperature
  /// by multiplying by the Exner pressure.
  ///
  /// When the model layer is not completely covered by observations,
  /// a potential temperature observation-minus-background increment is computed
  /// using linear interpolation of temperature between the layer boundaries.
  /// This increment is added to the background value to produce
  /// the averaged observation value.
  class ProfileAverageTemperature : public ProfileCheckBase {
   public:
    explicit ProfileAverageTemperature(const ConventionalProfileProcessingParameters &options);

    /// Average temperature observations onto model levels and store the results.
    /// \throws eckit::BadValue if vectors produced by the AveragePressure routine
    /// are not present.
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_potential_temperature});}

    /// List of names of GeoVaLs used in check validation.
    oops::Variables getValidationGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_air_temperature,
            ufo::VariableNames::geovals_testreference_air_temperature,
            ufo::VariableNames::geovals_testreference_air_temperature_qcflags
            });}

   private:
    /// Run check on a profile in the original ObsSpace and
    /// put the averaged data into the corresponding profile in the extended ObsSpace.
    void runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                            ProfileDataHolder &profileExtended);
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGETEMPERATURE_H_
