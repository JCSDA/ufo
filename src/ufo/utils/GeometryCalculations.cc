/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "eckit/exception/Exceptions.h"

#include "ufo/utils/Constants.h"
#include "ufo/utils/GeometryCalculations.h"

namespace ufo {

double haversine(const double lat1,
                 const double lon1,
                 const double lat2,
                 const double lon2) {
  const double earthRadius = ufo::Constants::mean_earth_rad * 1000.0;
  const double deg2rad = ufo::Constants::deg2rad;
  const double deltaLat = lat1 - lat2;
  const double deltaLon = lon1 - lon2;
  const double a =
    std::pow(std::sin(0.5 * deltaLat * deg2rad), 2.0) +
    std::pow(std::sin(0.5 * deltaLon * deg2rad), 2.0) *
    std::cos(lat1 * deg2rad) * std::cos(lat2 * deg2rad);
  const double c = std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
  return 2.0 * earthRadius * c;
}

void convertRangeAzimToLatLon(const double range,
                              const double azim,
                              const double stationLatitude,
                              const double stationLongitude,
                              const double stationElevation,
                              const double beamTiltAngle,
                              double & obLat,
                              double & obLon) {
  if (std::abs(stationLatitude) > 89.9) {
    throw eckit::UserError("Absolute value of station latitude cannot be "
                           "greater than 89.9 degrees", Here());
  }
  const double deg2rad = ufo::Constants::deg2rad;
  const double earthRadius = ufo::Constants::mean_earth_rad * 1000.0;
  const double earthCircumference = 2.0 * M_PI * earthRadius;
  // Earth radius accounting for curvature.
  const double effectiveEarthRadius = earthRadius * 4.0 / 3.0;
  const double height =
    std::sqrt(std::pow(range, 2.0) +
              std::pow(effectiveEarthRadius, 2.0) +
              2.0 * range * effectiveEarthRadius * std::sin(beamTiltAngle * deg2rad)) -
    effectiveEarthRadius + stationElevation;
  const double arcLength =
    effectiveEarthRadius * std::asin((range * std::cos(beamTiltAngle * deg2rad)) /
                                     (effectiveEarthRadius + height));

  obLat = stationLatitude + arcLength * std::cos(azim * deg2rad) * 360.0 / earthCircumference;
  obLon = stationLongitude + arcLength * std::sin(azim * deg2rad) * 360.0 /
    (earthCircumference * std::cos(stationLatitude * deg2rad));

  // Ensure longitude lies between 0 and 360 degrees.
  if (obLon > 360.0) {
    obLon -= 360.0;
  } else if (obLon < 0.0) {
    obLon += 360.0;
  }
}

std::pair<double, double> convertRangeAzimToLatLon(const double range,
                                                   const double azim,
                                                   const double stationLatitude,
                                                   const double stationLongitude,
                                                   const double stationElevation,
                                                   const double beamTiltAngle) {
  double obLat, obLon;
  convertRangeAzimToLatLon(range, azim, stationLatitude, stationLongitude, stationElevation,
                           beamTiltAngle, obLat, obLon);
  return {obLat, obLon};
}

}  // namespace ufo
