/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_GEOMETRYCALCULATIONS_H_
#define UFO_UTILS_GEOMETRYCALCULATIONS_H_

#include <utility>

namespace ufo {

/// The haversine function computes the great-circle distance between two points on the globe.
///
/// Function parameters:
/// * Latitude of the first point [deg].
/// * Longitude of the first point [deg].
/// * Latitude of the second point [deg].
/// * Longitude of the second point [deg].
/// Return value:
/// * Great-circle distance [m].
double haversine(const double,
                 const double,
                 const double,
                 const double);

/// Compute the latitude and longitude of an observation in a beam (typically a radar gate)
/// given the 3D location of the originating station, the range and azimuth of the observation
/// relative to the station, and the tilt angle of the beam relative to the horizontal.
///
/// Function parameters:
/// * Observation distance from the station [m].
/// * Observation azimuth relative to true North [deg].
/// * Station latitude [deg].
/// * Station longitude [deg].
/// * Station elevation above sea level [m].
/// * Beam tilt angle relative to the horizontal [deg].
/// Output values modified in place:
/// * Calculated observation latitude [deg].
/// * Calculated observation longitude [deg].
void convertRangeAzimToLatLon(const double,
                              const double,
                              const double,
                              const double,
                              const double,
                              const double,
                              double &,
                              double &);

/// Compute the latitude and longitude of an observation in a beam (typically a radar gate)
/// given the 3D location of the originating station, the range and azimuth of the observation
/// relative to the station, and the tilt angle of the beam relative to the horizontal.
/// The computation in this function is identical to that used in the function above,
/// but for convenience the values are returned in a pair rather than modified in place.
///
/// Function parameters:
/// * Observation distance from the station [m].
/// * Observation azimuth relative to true North [deg].
/// * Station latitude [deg].
/// * Station longitude [deg].
/// * Station elevation above sea level [m].
/// * Beam tilt angle relative to the horizontal [deg].
/// Return values (in a pair):
/// * Calculated observation latitude [deg].
/// * Calculated observation longitude [deg].
std::pair<double, double> convertRangeAzimToLatLon(const double,
                                                   const double,
                                                   const double,
                                                   const double,
                                                   const double,
                                                   const double);

}  // namespace ufo

#endif  // UFO_UTILS_GEOMETRYCALCULATIONS_H_
