/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBRADARTEMPLATE_H_
#define UFO_SUPEROB_SUPEROBRADARTEMPLATE_H_

#include <Eigen/Dense>

#include <utility>

namespace ufo {

/// Contains information on a radar scan.
struct SuperObRadarTemplateData {
  /// Radar station latitude [deg].
  double stationLatitude;
  /// Radar station longitude [deg].
  double stationLongitude;
  /// Radar station elevation above mean sea level [m].
  double stationElevation;
  /// Radar beam tilt angle relative to horizontal [deg].
  double beamTiltAngle;
  /// Minimum distance of radar gates from the instrument [m].
  double minGateRange;
  /// Number of beams in this scan.
  int numBeams;
  /// Number of gates in a beam (constant for the entire scan).
  int numGates;
  /// Width of radar beam [deg].
  double beamWidth;
  /// Width of radar gate [m].
  double gateWidth;
};

/// \brief Utility class that produces a template for use in radar superobbing.
///
/// \details This class can be used in conjunction with a radar superobbing algorithm.
/// Given a radar scan comprised of several beams, each of which contains a certain number of gates,
/// this utility will produce a template that can be used by the superobbing algorithm to compute
/// the value of the superob at each location.
///
/// The template spans a chosen number of beams in the azimuthal direction and a chosen distance
/// in the radial direction. The template is subdivided into adjacent cells which span the beam
/// width in the azimuthal direction and the gate width in the radial direction.
/// The superob is computed by averaging O-B values in all of the cells in the template region.
/// However, not all cells inside the region are used in the superob calculation;
/// cells closer to the edge of the region are rejected.
/// The spatial centre of the region is determined and a circle is defined with radius equal to the
/// smallest distance from the centre to the edge of the region.
/// Only points inside that circle are used in the superob calculation.
///
/// This class has two functions, both of which return a two-dimensional Eigen array.
/// The first function (`getMask`) returns a mask of zeros and ones, and the second (`getDistance`)
/// returns an array of distances from the central point of the region.
/// Both of these arrays can then be used in the calling superob routine.
/// The template produced in this class is identical for each group of beams, so only has to be
/// computed once. It can then be mapped around the entire scan.
///
/// The parameters of the constructor are:
/// * numBeamsInSuperObRegion: the number of beams in each superob region.
/// * superObRegionRadialExtent: the extent [m] of each superob region in the radial direction.
/// * scanData: metadata associated with the radar scan.
class SuperObRadarTemplate {
 public:
  SuperObRadarTemplate(const int,
                       const double,
                       const SuperObRadarTemplateData &);

  /// Return a mask which indicates whether or not to use values inside each cell to compute
  /// a superob.
  Eigen::ArrayXXi getMask() const {return mask_;}
  /// Return the distance of each cell to the centre of the superob template region.
  Eigen::ArrayXXd getDistance() const {return distance_;}

 private:
  Eigen::ArrayXXi mask_;
  Eigen::ArrayXXd distance_;
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBRADARTEMPLATE_H_
