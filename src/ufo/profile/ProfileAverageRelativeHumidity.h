/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEAVERAGERELATIVEHUMIDITY_H_
#define UFO_PROFILE_PROFILEAVERAGERELATIVEHUMIDITY_H_

#include <algorithm>
#include <cmath>
#include <vector>

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

  /// \brief Profile QC: average relative humidity observations onto model levels.
  ///
  /// Vectors produced by the AveragePressure routine must be present
  /// otherwise the exception eckit::BadValue will be thrown.
  ///
  /// By default, relative humidities are interpolated onto model layer boundaries rather than
  /// averaged across layers in order to avoid unwanted smoothing.
  /// This behaviour can be controlled with the \p AvgRH_Interp option.
  ///
  /// The interpolated/averaged relative humidity values are rejected at any layer where
  /// the averaged temperature value is less than or equal to the threshold \p AvgRH_AvgTThreshold.
  /// This threshold can be modified to an instrument-dependent value with the
  /// parameter \p AvgRH_InstrTThresholds, which is a map between WMO sonde instrument codes
  /// and the associated temperature thresholds.
  class ProfileAverageRelativeHumidity : public ProfileCheckBase {
   public:
    explicit ProfileAverageRelativeHumidity(const ConventionalProfileProcessingParameters &options);

    /// Average relative humidity observations onto model levels and store the results.
    /// \throws eckit::BadValue if vectors produced by the AveragePressure routine
    /// are not present.
    void runCheck(ProfileDataHandler &profileDataHandler) override;

    /// Run this check on the entire sample?
    bool runOnEntireSample() override {return true;}

    /// List of names of required GeoVaLs.
    oops::Variables getGeoVaLNames() override {
      return oops::Variables({ufo::VariableNames::geovals_relative_humidity});}

    /// List of names of GeoVaLs used in check validation.
    oops::Variables getValidationGeoVaLNames() override {
      return oops::Variables({oops::Variable
                              {ufo::VariableNames::geovals_testreference_relative_humidity},
            oops::Variable{ufo::VariableNames::geovals_testreference_relative_humidity_qcflags}
            });}

   private:
    /// Run check on a profile in the original ObsSpace and
    /// put the averaged data into the corresponding profile in the extended ObsSpace.
    void runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                            ProfileDataHolder &profileExtended);
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEAVERAGERELATIVEHUMIDITY_H_
