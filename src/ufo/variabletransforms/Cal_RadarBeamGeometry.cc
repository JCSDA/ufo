/*
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_RadarBeamGeometry.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"


namespace ufo {

/************************************************************************************/
//  Cal_RadarBeamGeometry
/************************************************************************************/

static TransformMaker<Cal_RadarBeamGeometry> makerCal_RadarBeamGeometry_("RadarBeamGeometry");

Cal_RadarBeamGeometry::Cal_RadarBeamGeometry(
        const Parameters_ & options,
        const ObsFilterData & data,
        const std::shared_ptr<ioda::ObsDataVector<int>> & flags,
        const std::shared_ptr<ioda::ObsDataVector<float>> & obserr)
  : TransformBase(options, data, flags, obserr),
    params_(options)
{}

/************************************************************************************/
/*
 Calculate radar beam geometry
*/
void Cal_RadarBeamGeometry::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " Calculate radar beam geometry" << std::endl;

  SetUseValidDataOnly(true);

  const size_t nlocs = obsdb_.nlocs();
  const float missing = util::missingValue<float>();
  const float deg2rad = Constants::deg2rad;

  // Compute effective Earth radius. If OPS compatibility mode is used then
  // there are two slightly different values used at different points.
  const double effectiveEarthRadius = params_.opsCompatibilityMode ?
    ufo::Constants::mean_earth_rad * 1000.0 * 1.33 :
    ufo::Constants::mean_earth_rad * 1000.0 * 4/3;
  const double effectiveEarthRadiusBeamTilt = params_.opsCompatibilityMode ?
    ufo::Constants::mean_earth_rad * 1000.0 * 1.3 :
    effectiveEarthRadius;

  // Input data
  std::vector<float> beamTilt(nlocs);  // Expect this angle in degrees
  std::vector<float> beamAzimuth(nlocs);  // Expect this angle in degrees
  std::vector<float> gateRange(nlocs);
  std::vector<float> radarAltitude(nlocs);

  getObservation("MetaData", "beamTiltAngle", beamTilt, true);
  getObservation("MetaData", "beamAzimuthAngle", beamAzimuth, true);
  getObservation("MetaData", "gateRange", gateRange, true);
  getObservation("MetaData", "stationElevation", radarAltitude, true);

  // Output data
  std::vector<float> sinTilt(nlocs);
  std::vector<float> cosAzimuthCosTilt(nlocs);
  std::vector<float> sinAzimuthCosTilt(nlocs);
  std::vector<float> gateHeight(nlocs);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (beamTilt[jloc] == missing ||
        beamAzimuth[jloc] == missing ||
        gateRange[jloc] == missing ||
        radarAltitude[jloc] == missing) {
      sinTilt[jloc] = missing;
      cosAzimuthCosTilt[jloc] = missing;
      sinAzimuthCosTilt[jloc] = missing;
      gateHeight[jloc] = missing;
      continue;
    }

    const double initialBeamTilt = beamTilt[jloc] * deg2rad;
    const double initialBeamAzimuth = beamAzimuth[jloc] * deg2rad;
    const double sinInitialBeamTilt = std::sin(initialBeamTilt);
    const double cosInitialBeamTilt = std::cos(initialBeamTilt);

    // Height of observation accounting for Earth's curvature
    const double observationHeightSquared = std::pow(gateRange[jloc], 2.0) +
      2.0 * effectiveEarthRadius * gateRange[jloc] * sinInitialBeamTilt +
      std::pow(effectiveEarthRadius, 2.0);

    gateHeight[jloc] = std::sqrt(observationHeightSquared) -
      effectiveEarthRadius + radarAltitude[jloc];

    // Effective angles accounting for earth's curvature
    const double beamTiltCorrection =
      std::atan(gateRange[jloc] * cosInitialBeamTilt /
                (gateRange[jloc] * sinInitialBeamTilt +
                 effectiveEarthRadiusBeamTilt + radarAltitude[jloc]));
    const double finalBeamTilt = initialBeamTilt + beamTiltCorrection;
    sinTilt[jloc] = std::sin(finalBeamTilt);
    const float cosFinalBeamTilt = std::cos(finalBeamTilt);
    cosAzimuthCosTilt[jloc] = std::cos(initialBeamAzimuth) * cosFinalBeamTilt;
    sinAzimuthCosTilt[jloc] = std::sin(initialBeamAzimuth) * cosFinalBeamTilt;
  }

  obsdb_.put_db("MetaData", "sinTilt", sinTilt);
  obsdb_.put_db("MetaData", "cosAzimuthCosTilt", cosAzimuthCosTilt);
  obsdb_.put_db("MetaData", "sinAzimuthCosTilt", sinAzimuthCosTilt);
  obsdb_.put_db("MetaData", "height", gateHeight);
}

}  // namespace ufo
