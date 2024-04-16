/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "eckit/exception/Exceptions.h"

#include "ufo/superob/SuperObRadarTemplate.h"

#include "ufo/utils/Constants.h"
#include "ufo/utils/GeometryCalculations.h"

namespace ufo {

SuperObRadarTemplate::SuperObRadarTemplate(const int numBeamsInSuperObRegion,
                                           const double superObRegionRadialExtent,
                                           const SuperObRadarTemplateData & scanData) {
  if (numBeamsInSuperObRegion < 1 || numBeamsInSuperObRegion > scanData.numBeams) {
    throw eckit::UserError("Invalid number of beams in the superob region", Here());
  }
  if (superObRegionRadialExtent < 1.0 ||
      superObRegionRadialExtent > scanData.numGates * scanData.gateWidth) {
    throw eckit::UserError("Invalid superob region radial extent", Here());
  }

  // Number of gates in each superob template region.
  const int numGatesInSuperObRegion =
    static_cast<int>(scanData.numGates * scanData.gateWidth / superObRegionRadialExtent);

  // Azimuthal extent of superob template region.
  const double superObRegionAzimuthalExtent = numBeamsInSuperObRegion * scanData.beamWidth;

  // Fill arrays of latitude and longitude for each cell in the superob template region.
  Eigen::ArrayXXd templateLat(scanData.numGates, numBeamsInSuperObRegion);
  Eigen::ArrayXXd templateLon(scanData.numGates, numBeamsInSuperObRegion);
  for (int i = 0; i < scanData.numGates; ++i) {
    for (int j = 0; j < numBeamsInSuperObRegion; ++j) {
      const double range = (i + 0.5) * scanData.gateWidth + scanData.minGateRange;
      const double azim = (j + 0.5) * scanData.beamWidth;
      ufo::convertRangeAzimToLatLon(range,
                                    azim,
                                    scanData.stationLatitude,
                                    scanData.stationLongitude,
                                    scanData.stationElevation,
                                    scanData.beamTiltAngle,
                                    templateLat(i, j),
                                    templateLon(i, j));
    }
  }

  // Arrays holding superob template information.
  // These arrays are only computed for one template (i.e. just one set of beams).
  // Once filled, these arrays are copied into all other templates around the scan.
  Eigen::ArrayXXi mask = Eigen::ArrayXXi::Zero(scanData.numGates, numBeamsInSuperObRegion);
  Eigen::ArrayXXd distance =
    Eigen::ArrayXXd::Zero(scanData.numGates, numBeamsInSuperObRegion);

  // Loop over each gate in each superob region.
  // Determine the radius of the circle with which the mask is computed,
  // then fill the mask.
  double range = 0.0;  // Range of superob region [m].
  for (int irang = 0; irang < numGatesInSuperObRegion; ++irang) {
    // Center of superob region.
    const auto[latSuperObRegionCentre, lonSuperObRegionCentre] =
      ufo::convertRangeAzimToLatLon(range + 0.5 * superObRegionRadialExtent,
                                    0.5 * superObRegionAzimuthalExtent,
                                    scanData.stationLatitude,
                                    scanData.stationLongitude,
                                    scanData.stationElevation,
                                    scanData.beamTiltAngle);
    // Upper edge of superob region in the radial direction.
    const auto[latSuperObRegionUpperRange, lonSuperObRegionUpperRange] =
      ufo::convertRangeAzimToLatLon(range + superObRegionRadialExtent,
                                    0.5 * superObRegionAzimuthalExtent,
                                    scanData.stationLatitude,
                                    scanData.stationLongitude,
                                    scanData.stationElevation,
                                    scanData.beamTiltAngle);
    // Upper edge of superob region in the azimuthal direction.
    const auto[latSuperObRegionUpperAzim, lonSuperObRegionUpperAzim] =
      ufo::convertRangeAzimToLatLon(range + 0.5 * superObRegionRadialExtent,
                                    superObRegionAzimuthalExtent,
                                    scanData.stationLatitude,
                                    scanData.stationLongitude,
                                    scanData.stationElevation,
                                    scanData.beamTiltAngle);

    // Determine reference distance, equal to whichever is the smaller
    // of the two distances between the upper edges and the centre
    // of the superob region.
    const double acceptanceRegionRadiusCandidate1 =
      ufo::haversine(latSuperObRegionCentre,
                     lonSuperObRegionCentre,
                     latSuperObRegionUpperRange,
                     lonSuperObRegionUpperRange);
    const double acceptanceRegionRadiusCandidate2 =
      ufo::haversine(latSuperObRegionCentre,
                     lonSuperObRegionCentre,
                     latSuperObRegionUpperAzim,
                     lonSuperObRegionUpperAzim);
    const double acceptanceRegionRadius =
      acceptanceRegionRadiusCandidate2 > acceptanceRegionRadiusCandidate1 ?
      acceptanceRegionRadiusCandidate1 : acceptanceRegionRadiusCandidate2;

    // Fill in template arrays based on whether each cell lies within a circle,
    // centred on the cell centre, with radius as calculated above.
    for (size_t i = 0; i < scanData.numGates; ++i) {
      for (size_t j = 0; j < numBeamsInSuperObRegion; ++j) {
        const double superObCellDistance = ufo::haversine(latSuperObRegionCentre,
                                                          lonSuperObRegionCentre,
                                                          templateLat(i, j),
                                                          templateLon(i, j));
        if (superObCellDistance < acceptanceRegionRadius) {
          mask(i, j) = 1;
          distance(i, j) = superObCellDistance;
        }
      }
    }

    // Move to the next superob region in the beam.
    range += superObRegionRadialExtent;
  }

  // Copy the template to all of the other regions in the scan.
  // Add buffer zone around the edges. This matches the Met Office OPS code.
  // todo(ctgh): Check whether the choice of extension could be modified.
  const int imax = scanData.numGates + 2;
  const int jmax = scanData.numBeams + 3;
  mask_ = Eigen::ArrayXXi::Zero(imax, jmax);
  distance_ = Eigen::ArrayXXd::Zero(imax, jmax);
  for (size_t i = 0; i < scanData.numBeams; i += numBeamsInSuperObRegion) {
    for (size_t j = 0; j < scanData.numGates; ++j) {
      for (size_t k = i; k < i + numBeamsInSuperObRegion; ++k) {
        mask_(j + 1, k + 1) = mask(j, k - i);
        distance_(j + 1, k + 1) = distance(j, k - i);
      }
    }
  }
}

}  // namespace ufo
