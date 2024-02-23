/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Range.h"
#include "test/TestEnvironment.h"
#include "ufo/fov/FieldOfView.interface.h"
#include "ufo/fov/ReduceOverFieldOfView.h"
#include "ufo/fov/SampleFieldOfView.h"

namespace ufo {
namespace test {

/// Parameters defining the FOV and giving the reference to test against
class FovSamplerTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FovSamplerTestParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> sensor{"sensor", this};
  oops::RequiredParameter<std::string> satellite{"satellite", this};

  // scan position is for cross-track scanning instruments only
  oops::Parameter<std::vector<int>> scan_positions{"scan positions", std::vector<int>{}, this};
  oops::RequiredParameter<std::vector<float>> sensor_azimuth_angles{"sensor azimuth angles", this};
  oops::RequiredParameter<std::vector<float>> longitudes{"longitudes", this};
  oops::RequiredParameter<std::vector<float>> latitudes{"latitudes", this};
  oops::RequiredParameter<std::vector<util::DateTime>> times{"times", this};

  // For each obs, reference values of...
  // min/max sample index
  oops::RequiredParameter<std::vector<int>> reference_range_begins{"reference range begins", this};
  oops::RequiredParameter<std::vector<int>> reference_range_ends{"reference range ends", this};

  // min/avg/max of lon, lat, and weight
  oops::RequiredParameter<std::vector<double>> reference_lon_mins{"reference longitude mins", this};
  oops::RequiredParameter<std::vector<double>> reference_lon_avgs{"reference longitude avgs", this};
  oops::RequiredParameter<std::vector<double>> reference_lon_maxs{"reference longitude maxs", this};
  oops::RequiredParameter<std::vector<double>> reference_lat_mins{"reference latitude mins", this};
  oops::RequiredParameter<std::vector<double>> reference_lat_avgs{"reference latitude avgs", this};
  oops::RequiredParameter<std::vector<double>> reference_lat_maxs{"reference latitude maxs", this};
  oops::RequiredParameter<std::vector<double>> reference_weight_mins{"reference weight mins", this};
  oops::RequiredParameter<std::vector<double>> reference_weight_avgs{"reference weight avgs", this};
  oops::RequiredParameter<std::vector<double>> reference_weight_maxs{"reference weight maxs", this};

  // absolute tolerance for FOV tests
  oops::RequiredParameter<double> abs_tol{"abs tol", this};
};

// -----------------------------------------------------------------------------

void fovSamplerTestHelper(const std::string sensor, const std::string satellite,
                          const std::vector<int>& scan_positions,
                          const std::vector<float>& sensor_azimuth_angles,
                          const std::vector<float>& longitudes,
                          const std::vector<float>& latitudes,
                          const std::vector<util::DateTime>& times,
                          const std::vector<util::Range<size_t>>& reference_ranges,
                          const std::vector<double>& reference_lon_mins,
                          const std::vector<double>& reference_lon_avgs,
                          const std::vector<double>& reference_lon_maxs,
                          const std::vector<double>& reference_lat_mins,
                          const std::vector<double>& reference_lat_avgs,
                          const std::vector<double>& reference_lat_maxs,
                          const std::vector<double>& reference_weight_mins,
                          const std::vector<double>& reference_weight_avgs,
                          const std::vector<double>& reference_weight_maxs,
                          const double abs_tol) {
  const int sample_points_per_semi_axis = 4;

  std::vector<float> sample_lons_floats;
  std::vector<float> sample_lats_floats;
  std::vector<util::DateTime> sample_times;
  std::vector<util::Range<size_t>> sample_ranges;
  std::vector<double> sample_weights;
  fov::getSampleLocationsAndWeights(sample_lons_floats, sample_lats_floats, sample_times,
                                    sample_ranges, sample_weights,
                                    sensor, satellite,
                                    sample_points_per_semi_axis,
                                    longitudes, latitudes, times,
                                    scan_positions, sensor_azimuth_angles);

  const int nlocs = longitudes.size();
  const int nsamples = sample_lons_floats.size();

  // Cast floats to doubles for easier comparisons below (until templated lambdas!)
  std::vector<double> sample_lons(nsamples);
  for (size_t i = 0; i < nsamples; ++i) {
    sample_lons[i] = static_cast<double>(sample_lons_floats[i]);
  }
  std::vector<double> sample_lats(nsamples);
  for (size_t i = 0; i < nsamples; ++i) {
    sample_lats[i] = static_cast<double>(sample_lats_floats[i]);
  }

  // Sanity checks on returned vectors
  EXPECT(sample_ranges.size() == nlocs);  // one range per obs
  EXPECT(sample_lats.size() == nsamples);
  EXPECT(sample_times.size() == nsamples);
  EXPECT(sample_weights.size() == nsamples);

  // Check sample ranges vs expected
  for (size_t i = 0; i < nlocs; ++i) {
    EXPECT(reference_ranges[i].begin == sample_ranges[i].begin);
    EXPECT(reference_ranges[i].end == sample_ranges[i].end);
  }

  // Check times, all samples should share obs time
  for (size_t i = 0; i < nlocs; ++i) {
    for (size_t isample = sample_ranges[i].begin; isample < sample_ranges[i].end; ++isample) {
      EXPECT(times[i] == sample_times[isample]);
    }
  }

  const auto min_per_fov = [&](std::vector<double> & data) -> std::vector<double> {
    std::vector<double> min(nlocs);
    for (size_t i = 0; i < nlocs; ++i) {
      min[i] = *std::min_element(data.begin() + sample_ranges[i].begin,
                                 data.begin() + sample_ranges[i].end);
    }
    return min;
  };
  const auto avg_per_fov = [&](std::vector<double> & data) -> std::vector<double> {
    std::vector<double> avg(nlocs);
    for (size_t i = 0; i < nlocs; ++i) {
      avg[i] = std::accumulate(data.begin() + sample_ranges[i].begin,
                               data.begin() + sample_ranges[i].end,
                               0.0) / (sample_ranges[i].end - sample_ranges[i].begin);
    }
    return avg;
  };
  const auto max_per_fov = [&](std::vector<double> & data) -> std::vector<double> {
    std::vector<double> max(nlocs);
    for (size_t i = 0; i < nlocs; ++i) {
      max[i] = *std::max_element(data.begin() + sample_ranges[i].begin,
                                 data.begin() + sample_ranges[i].end);
    }
    return max;
  };

  // Sample locations and weights are too numerous to check individually, so check statistics
  EXPECT(oops::are_all_close_absolute(reference_lon_mins, min_per_fov(sample_lons), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_lon_avgs, avg_per_fov(sample_lons), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_lon_maxs, max_per_fov(sample_lons), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_lat_mins, min_per_fov(sample_lats), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_lat_avgs, avg_per_fov(sample_lats), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_lat_maxs, max_per_fov(sample_lats), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_weight_mins, min_per_fov(sample_weights), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_weight_avgs, avg_per_fov(sample_weights), abs_tol));
  EXPECT(oops::are_all_close_absolute(reference_weight_maxs, max_per_fov(sample_weights), abs_tol));
}

// -----------------------------------------------------------------------------

// Check that field of view samples and weights match expected references
// This function mostly is in charge of reading config, setting up params, then
// delegates to the test helper function for the test work and comparisons.
void testSampleFieldOfView() {
  const eckit::Configuration& conf = ::test::TestEnvironment::config();
  const eckit::LocalConfiguration fovConf = conf.getSubConfiguration("field of view sampling");
  FovSamplerTestParameters params{};
  params.validateAndDeserialize(fovConf);

  // Sanity check on inputs
  const int nlocs = params.longitudes.value().size();
  ASSERT(params.scan_positions.value().empty() || params.scan_positions.value().size() == nlocs);
  ASSERT(params.sensor_azimuth_angles.value().size() == nlocs);
  ASSERT(params.latitudes.value().size() == nlocs);
  ASSERT(params.times.value().size() == nlocs);

  ASSERT(params.reference_lon_mins.value().size() == nlocs);
  ASSERT(params.reference_lon_avgs.value().size() == nlocs);
  ASSERT(params.reference_lon_maxs.value().size() == nlocs);
  ASSERT(params.reference_lat_mins.value().size() == nlocs);
  ASSERT(params.reference_lat_avgs.value().size() == nlocs);
  ASSERT(params.reference_lat_maxs.value().size() == nlocs);
  ASSERT(params.reference_weight_mins.value().size() == nlocs);
  ASSERT(params.reference_weight_avgs.value().size() == nlocs);
  ASSERT(params.reference_weight_maxs.value().size() == nlocs);

  std::vector<util::Range<size_t>> reference_ranges(nlocs);
  for (size_t i = 0; i < nlocs; ++i) {
    reference_ranges[i].begin = params.reference_range_begins.value()[i];
    reference_ranges[i].end = params.reference_range_ends.value()[i];
  }

  fovSamplerTestHelper(params.sensor, params.satellite,
                       params.scan_positions, params.sensor_azimuth_angles,
                       params.longitudes, params.latitudes, params.times,
                       reference_ranges,
                       params.reference_lon_mins,
                       params.reference_lon_avgs,
                       params.reference_lon_maxs,
                       params.reference_lat_mins,
                       params.reference_lat_avgs,
                       params.reference_lat_maxs,
                       params.reference_weight_mins,
                       params.reference_weight_avgs,
                       params.reference_weight_maxs,
                       params.abs_tol);
}

// -----------------------------------------------------------------------------

// Reduction helpers are simple algebra and can be tested vs analytic computations
void testReduceOverFieldOfView() {
  const double tol = 1e-15;

  // Make up fake samples for 3 obs with {4,5,6} samples each
  const int nlocs = 3;
  const int nsamples = 15;
  const std::vector<util::Range<size_t>> ranges({{{0, 4}, {4, 9}, {9, 15}}});
  std::vector<double> samples({{0.2, 0.3, 0.6, 0.9,
                                1.7, 1.2, 1.9, 4.2, 2.9,
                                1020.2, 1982.3, 1375.9, 4082.3, 2159.8, 2901.1}});
  std::vector<double> weights({{1.0, 1.0, 1.0, 1.0,
                                0.2, 0.3, 0.4, 0.5, 0.6,
                                100.0, 1.0, 1.0, 1.0, 1.0, 1.0}});

  std::vector<double> result(nlocs);
  fov::average(result, ranges, samples, weights);

  EXPECT(oops::is_close_absolute(result[0], 0.5, tol));
  EXPECT(oops::is_close_absolute(result[1], 2.65, tol));
  EXPECT(oops::is_close_absolute(result[2], 1090.68, tol));

  // Fake masks for these 3 obs: obs0 is unmasked, obs1 is partly masked, obs2 is fully masked away
  std::vector<double> mask_per_obs({{1.0, 0.2, 0.0}});
  std::vector<double> mask_per_sample({{1.0, 1.0, 1.0, 1.0,
                                        0.8, 0.2, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0}});
  const double fallback = 1.23;
  fov::averageWithMask(result, ranges, samples, weights, mask_per_obs, mask_per_sample, fallback);

  EXPECT(oops::is_close_absolute(result[0], 0.5, tol));
  EXPECT(oops::is_close_absolute(result[1], 1.5636363636363636, tol));
  EXPECT(oops::is_close_absolute(result[2], 1.23, tol));

  // Similar for identifying the dominant (most heavily weighted) integer (but stored as doubles..)
  std::vector<double> int_samples({{1.0, 2.0, 1.0, 3.0,
                                    1.0, 2.0, 1.0, 3.0, 1.0,
                                    1.0, 2.0, 1.0, 3.0, 1.0, 4.0}});
  const int int_fallback = 6;
  fov::dominantIntegerValue(result, ranges, int_samples, weights, mask_per_obs, mask_per_sample,
                            int_fallback);

  EXPECT(result[0] == 1.0);
  EXPECT(result[1] == 1.0);
  EXPECT(result[2] == 6.0);
}

// -----------------------------------------------------------------------------

class SampleAndReduceOverFieldOfView : public oops::Test {
 public:
  SampleAndReduceOverFieldOfView() = default;
  ~SampleAndReduceOverFieldOfView() = default;

 private:
  std::string testid() const override { return "ufo::test::SampleAndReduceOverFieldOfView"; }

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/fov/testSampleFieldOfView") { testSampleFieldOfView(); });
    ts.emplace_back(CASE("ufo/fov/testReduceOverFieldOfView") { testReduceOverFieldOfView(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo
