/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/fov/SampleFieldOfView.h"

#include <algorithm>
#include <mutex>  // NOLINT(build/c++11)
#include <ostream>
#include <utility>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Range.h"

#include "ufo/fov/FieldOfView.interface.h"

namespace {

// A global mutex, matching the scope of the GSI memory buffers we're trying to protect
std::mutex gsiFovGlobalsMutex;

std::pair<std::vector<float>, std::vector<float>> lonlat_grid_over_ellipse(
    const size_t sample_points_per_semi_axis,
    const float lon, const float lat, const std::vector<double>& ellipse_lons,
    const std::vector<double>& ellipse_lats) {
  ASSERT(sample_points_per_semi_axis > 1);
  ASSERT(sample_points_per_semi_axis < 10);  // avoid expensive config with diminishing returns
  ASSERT(lon >= 0.0);
  ASSERT(lon < 360.0);
  ASSERT(lat >= -90.0);
  ASSERT(lat < 90.0);

  // Check for ellipse center too close to north/south pole, where the GSI routines have
  // coordinate problems.
  // - long term, should rewrite trig to avoid these coordinate issues.
  // - for now, we just return obs coords when obs is within a few degrees of pole
  {
    const float margin = 2.0;  // this is the safety margin used in GSI
    const float max_lat = 90.0 - margin;
    if (lat < -max_lat || lat > max_lat) {
      return std::make_pair(std::vector<float>{lon}, std::vector<float>{lat});
    }
  }

  // Grid extends around the central position symmetrically, so there is a point
  // at the central position. Find half-width of grid in each direction:
  // WARNING: the algorithm below will not work if the ellipse encloses the North or South pole
  const float max_lon = *std::max_element(ellipse_lons.begin(), ellipse_lons.end());
  const float min_lon = *std::min_element(ellipse_lons.begin(), ellipse_lons.end());
  const float max_lat = static_cast<float>(
      *std::max_element(ellipse_lats.begin(), ellipse_lats.end()));
  const float min_lat = static_cast<float>(
      *std::min_element(ellipse_lats.begin(), ellipse_lats.end()));

  const float safety = 0.99;
  const float grid_half_lon = safety * std::max(std::abs(max_lon - lon), std::abs(lon - min_lon));
  const float grid_half_lat = safety * std::max(std::abs(max_lat - lat), std::abs(lat - min_lat));

  const size_t points_per_direction = 2 * sample_points_per_semi_axis - 1;
  std::vector<float> grid_lons(points_per_direction * points_per_direction);
  std::vector<float> grid_lats(points_per_direction * points_per_direction);
  for (size_t i = 0; i < points_per_direction; ++i) {
    const float this_lon = lon - grid_half_lon
                           + (2.0 * grid_half_lon) * i / (points_per_direction - 1.0);
    for (size_t j = 0; j < points_per_direction; ++j) {
      const float this_lat = lat - grid_half_lat
                             + (2.0 * grid_half_lat) * j / (points_per_direction - 1.0);
      grid_lons.push_back(this_lon);
      grid_lats.push_back(this_lat);
    }
  }
  return std::make_pair(grid_lons, grid_lats);
}
}  // namespace

namespace ufo {
namespace fov {

// -----------------------------------------------------------------------------

FieldOfViewSampler::FieldOfViewSampler(const std::string& sensor, const std::string& platform)
  : keyFov_(0), sensor_(sensor), platform_(platform), gsi_npoly_(0)
{
  oops::Log::trace() << "ufo::fov::FieldOfViewSampler constructor starting" << std::endl;

  bool gsi_valid_instr;

  const int sensor_len = sensor_.size();
  const char* sensor_cstr = sensor_.c_str();
  const int platform_len = platform_.size();
  const char* platform_cstr = platform_.c_str();
  ufo_fov_setup_f90(keyFov_, sensor_len, sensor_cstr, platform_len, platform_cstr,
      gsi_valid_instr, gsi_npoly_);

  // Make sure initialization was successful
  ASSERT(gsi_npoly_ > 0);
  ASSERT(gsi_valid_instr);

  oops::Log::trace() << "ufo::fov::FieldOfViewSampler constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

FieldOfViewSampler::~FieldOfViewSampler() {
  oops::Log::trace() << "ufo::fov::FieldOfViewSampler destructor starting" << std::endl;
  ufo_fov_delete_f90(keyFov_);
  oops::Log::trace() << "ufo::fov::FieldOfViewSampler destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldOfViewSampler::sampleFieldOfView(
    size_t& sample_nlocs,
    std::vector<float>& sample_lons,
    std::vector<float>& sample_lats,
    std::vector<double>& sample_weights,
    const size_t sample_points_per_semi_axis,
    const int scan_position, const float sensor_azimuth_angle,
    const float longitude, const float latitude) const {
  oops::Log::trace() << "ufo::fov::FieldOfViewSampler::sampleFieldOfView starting" << std::endl;

  // Prepare output variables
  sample_nlocs = 0;
  sample_lons.clear();
  sample_lats.clear();
  sample_weights.clear();

  std::vector<double> ellipse_lons(gsi_npoly_);
  std::vector<double> ellipse_lats(gsi_npoly_);

  const int sensor_len = sensor_.size();
  const char* sensor_cstr = sensor_.c_str();
  ufo_fov_ellipse_f90(keyFov_, sensor_len, sensor_cstr, scan_position, sensor_azimuth_angle,
                      longitude, latitude, gsi_npoly_, ellipse_lons[0], ellipse_lats[0]);

  const auto grid = lonlat_grid_over_ellipse(sample_points_per_semi_axis,
                                             longitude, latitude, ellipse_lons, ellipse_lats);

  for (size_t i = 0; i < grid.first.size(); ++i) {
    const float lon = grid.first[i];  // lons vector
    const float lat = grid.second[i];  // lats vector
    double weight;

    ufo_antenna_power_within_fov_f90(keyFov_, sensor_len, sensor_cstr, scan_position,
                                     sensor_azimuth_angle, longitude, latitude, lon, lat, weight);

    // Antenna power > 0 means point is inside FOV ellipse
    //
    // Possible improvement: eliminate samples that have a positive but tiny power, because this
    // would eliminate CRTM calculations for sample locations that contribute little to the result.
    // Note this would require normalizing the weights and identifying a reasonable threshold.
    if (weight > 0.0) {
      sample_lons.push_back(lon);
      sample_lats.push_back(lat);
      sample_weights.push_back(weight);
      sample_nlocs++;
    }
  }

  oops::Log::trace() << "ufo::fov::FieldOfViewSampler::sampleFieldOfView done" << std::endl;
}

// -----------------------------------------------------------------------------

void getSampleLocationsAndWeights(
    std::vector<float> & sample_lons,
    std::vector<float> & sample_lats,
    std::vector<util::DateTime> & sample_times,
    std::vector<util::Range<size_t>> & sample_ranges,
    std::vector<double> & sample_weights,
    const std::string & sensor,
    const std::string & platform,
    const int sample_points_per_semi_axis,
    const std::vector<float> & lons,
    const std::vector<float> & lats,
    const std::vector<util::DateTime> & times,
    const std::vector<int> & scan_positions,
    const std::vector<float> & sensor_azimuth_angles) {
  if (!sample_lons.empty()) sample_lons.clear();
  if (!sample_lats.empty()) sample_lats.clear();
  if (!sample_times.empty()) sample_times.clear();
  if (!sample_ranges.empty()) sample_ranges.clear();
  if (!sample_weights.empty()) sample_weights.clear();

  ASSERT(scan_positions.size() == sensor_azimuth_angles.size());
  const size_t nlocs = scan_positions.size();

  // Lock the mutex, ensuring multiple threads can NOT simultaneously work with the global variables
  // in the underlying GSI FOV code.
  std::lock_guard<std::mutex> guard(gsiFovGlobalsMutex);

  const FieldOfViewSampler sampler(sensor, platform);

  size_t counter = 0;
  for (size_t i = 0; i < nlocs; ++i) {
    size_t num_samples = 0;
    std::vector<float> tmp_lons{};
    std::vector<float> tmp_lats{};
    std::vector<double> tmp_weights{};
    sampler.sampleFieldOfView(num_samples, tmp_lons, tmp_lats, tmp_weights,
                              sample_points_per_semi_axis,
                              scan_positions[i], sensor_azimuth_angles[i], lons[i], lats[i]);
    ASSERT(num_samples >= 1);  // the obs disappears if not sampled at all
    sample_lons.insert(sample_lons.end(), tmp_lons.begin(), tmp_lons.end());
    sample_lats.insert(sample_lats.end(), tmp_lats.begin(), tmp_lats.end());
    sample_times.insert(sample_times.end(), num_samples, times[i]);  // append num_samples copies
    sample_ranges.push_back({counter, counter + num_samples});
    sample_weights.insert(sample_weights.end(), tmp_weights.begin(), tmp_weights.end());
    counter += num_samples;
  }
}

// -----------------------------------------------------------------------------

}  // namespace fov
}  // namespace ufo
