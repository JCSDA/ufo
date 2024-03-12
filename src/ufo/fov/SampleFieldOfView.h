/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"

#include "ufo/Fortran.h"

namespace util {
class DateTime;
template <typename T> class Range;
}


namespace ufo {
namespace fov {

// -----------------------------------------------------------------------------

/*!
 * \brief Computes the sampling locations and weights for a series of observations, assuming a
 * space-based sensor with a conical field of view.
 *
 * A FieldOfViewSampler object is associated with a particular sensor and satellite platform. It
 * computes a series of sampling locations associated with each observation from this sensor, these
 * samples serving to represent the extended field of view in the OOPS GeoVaLs interface. It also
 * computes the sampling weight associated with each sampling location.
 *
 * In principle, any number of methods could be used to compute the sampling locations and weights.
 *
 * The current implementation is an interface to the GSI `calc_fov` modules, acting as follows:
 *
 * - each observation's field of view forms an ellipse on the earth surface, centered on the
 *   observation's "location". A GSI subroutine is used to compute the ellipse.
 * - within this ellipse, a lat-lon grid is used for the sampling locations. The resolution of this
 *   lat-lon grid is set by requesting a particular number of points per semi-axis of the
 *   field-of-view ellipse.
 * - at each sample (= each point of the lat-lon grid that falls within the ellipse), a GSI
 *   subroutine is used to compute the sampling weight (= the GSI antenna power).
 *
 * WARNINGS -- Some known shortcomings of this design are:
 *
 * - close to the poles (within 2 degrees), where the lat-lon grid has a singularity, a SINGLE
 *   sample is used at the obs location, so no field-of-view effects are accounted for.
 * - the GSI code used by the FieldOfViewSampler class uses module-level variables in its
 *   computations. This means the UFO code written on top of this GSI code should not be considered
 *   thread-safe, because all FieldOfViewSampler objects point to the SAME memory in the GSI code.
 *   For this reason, access to the FieldOfViewSampler class is currently channeled through a
 *   helper method which adds mutex locks; this restores thread-safety when called correctly
 *   through the helper method, at the expense of runtime performance.
 *
 * NOTES -- Some considerations for future development are:
 *
 * - the FieldOfViewSampler class currently delegates to GSI code for computing the ellipse and
 *   weights of each observation, but this is not a fundamental restriction. It's mostly for the
 *   convenience of not porting over the large lookup tables of geometric factors and antenna power
 *   coefficients for each sensor.
 * - this class could be generalized to allow for different algorithms when computing the sampling
 *   locations and weights. For example a much simpler model of the underlying observation could be
 *   acheived by sampling the field of view ellipse with a few points located close to the obs
 *   location (where the antenna power is typically largest) and using a Gaussian falloff for the
 *   antenna power. This may or may not lead to worse H(x) than the current more complicated code.
 *   The class could also be extended to support switching between different sampling algorithms.
 */
class FieldOfViewSampler : private util::ObjectCounter<FieldOfViewSampler> {
 public:
  static const std::string classname() { return "ufo::fov::FieldOfViewSampler"; }

  FieldOfViewSampler(const std::string& sensor, const std::string& platform);
  ~FieldOfViewSampler();

 private:
  // Marking `sampleFieldOfView` and private and then marking `getSampleLocationsAndWeights` is a
  // workaround to ensure the sampling locations and weights are accessed only via the
  // helper function getSampleLocationsAndWeights which uses a threadsafe calling pattern.
  // If in the future the use of GSI code is eliminated, all this can be simplified.
  friend void getSampleLocationsAndWeights(
      std::vector<float> & sample_lons,
      std::vector<float> & sample_lats,
      std::vector<util::DateTime> & sample_times,
      std::vector<util::Range<size_t>> & sample_ranges,
      std::vector<double> & sample_weights,
      const std::string & sensor,
      const std::string & platform,
      int sample_points_per_semi_axis,
      const std::vector<float> & lons,
      const std::vector<float> & lats,
      const std::vector<util::DateTime> & times,
      const std::vector<int> & scan_positions,
      const std::vector<float> & sensor_azimuth_angles);
  void sampleFieldOfView(size_t& sample_nlocs,
      std::vector<float>& sample_lons,
      std::vector<float>& sample_lats,
      std::vector<double>& sample_weights,
      size_t sample_points_per_semi_axis,
      int scan_position, float sensor_azimuth_angle,
      float longitude, float latitude) const;

 private:
  F90fov keyFov_;
  const std::string sensor_;
  const std::string platform_;
  int gsi_npoly_;
};

// -----------------------------------------------------------------------------

/// Compute FOV sample locations and weights for each obs. Returns locations and weights for all
/// obs contiguously in two vectors, with sampling range of each obs indicated by sample_ranges.
///
/// Expect following vector sizes:
/// - sample_ranges : resized to size nlocs; one per obs
/// - sample_* : resized to max(sample_ranges.end) on output; one per FOV sample over all obs
/// - lons : size nlocs; one per obs
/// - lats : size nlocs
/// - times : size nlocs
/// - scan_positions : size nlocs
/// - sensor_azimuth_angles : size nlocs
void getSampleLocationsAndWeights(
    std::vector<float> & sample_lons,
    std::vector<float> & sample_lats,
    std::vector<util::DateTime> & sample_times,
    std::vector<util::Range<size_t>> & sample_ranges,
    std::vector<double> & sample_weights,
    const std::string & sensor,
    const std::string & platform,
    int sample_points_per_semi_axis,
    const std::vector<float> & lons,
    const std::vector<float> & lats,
    const std::vector<util::DateTime> & times,
    const std::vector<int> & scan_positions,
    const std::vector<float> & sensor_azimuth_angles);

// -----------------------------------------------------------------------------

}  // namespace fov
}  // namespace ufo
