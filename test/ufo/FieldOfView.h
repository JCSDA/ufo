/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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
#include "test/TestEnvironment.h"
#include "ufo/fov/FieldOfView.interface.h"

namespace ufo {
namespace test {

/// Parameters defining the FOV and giving the reference to test against
class FovWrapperTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FovWrapperTestParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> sensor{"sensor", this};
  oops::RequiredParameter<std::string> satellite{"satellite", this};
  // scan position is for cross-track scanning instruments only
  oops::Parameter<int> scan_position{"scan position", -999, this};
  oops::RequiredParameter<float> sensor_azimuth_angle{"sensor azimuth angle", this};
  oops::RequiredParameter<float> longitude{"longitude", this};
  oops::RequiredParameter<float> latitude{"latitude", this};
  // reference longitudes of FOV ellipse polygon, from GSI call
  oops::RequiredParameter<std::vector<double>> reference_ellipse_longitudes{
      "reference ellipse longitudes", this};
  // reference latitudes of FOV ellipse polygon, from GSI call
  oops::RequiredParameter<std::vector<double>> reference_ellipse_latitudes{
      "reference ellipse latitudes", this};
  // longitudes at which to sample and test the antenna power
  oops::RequiredParameter<std::vector<float>> sample_longitudes{"sample longitudes", this};
  // latitudes at which to sample and test the antenna power
  oops::RequiredParameter<std::vector<float>> sample_latitudes{"sample latitudes", this};
  // reference antenna power values, from GSI call
  oops::RequiredParameter<std::vector<double>> reference_sample_powers{
      "reference sample powers", this};
  // absolute tolerance for FOV tests
  oops::RequiredParameter<double> abs_tol{"abs tol", this};
};

// -----------------------------------------------------------------------------

/// Check that field of view ellipse and antenna power match the references within abs_tol
void fovWrapperTestHelper(const std::string sensor, const std::string satellite,
                          const int scan_position,
                          const float sensor_azimuth_angle,
                          const float longitude, const float latitude,
                          const std::vector<double>& reference_ellipse_lons,
                          const std::vector<double>& reference_ellipse_lats,
                          const std::vector<float>& sample_lons,
                          const std::vector<float>& sample_lats,
                          const std::vector<double>& reference_sample_powers,
                          const double abs_tol) {
  const int sensor_len = sensor.size();
  const char* sensor_cstr = sensor.c_str();
  const int satellite_len = satellite.size();
  const char* satellite_cstr = satellite.c_str();
  F90fov key;
  bool gsi_valid_instr;
  int gsi_npoly;

  fov::ufo_fov_setup_f90(key, sensor_len, sensor_cstr, satellite_len, satellite_cstr,
                         gsi_valid_instr, gsi_npoly);

  EXPECT(gsi_valid_instr);
  EXPECT(gsi_npoly == 30);

  std::vector<double> ellipse_lons(gsi_npoly);
  std::vector<double> ellipse_lats(gsi_npoly);
  fov::ufo_fov_ellipse_f90(key, sensor_len, sensor_cstr, scan_position, sensor_azimuth_angle,
                           longitude, latitude, gsi_npoly, ellipse_lons[0], ellipse_lats[0]);

  EXPECT(oops::are_all_close_absolute(ellipse_lons, reference_ellipse_lons, abs_tol));
  EXPECT(oops::are_all_close_absolute(ellipse_lats, reference_ellipse_lats, abs_tol));

  for (size_t i = 0; i < reference_sample_powers.size(); ++i) {
    double sample_power;
    fov::ufo_antenna_power_within_fov_f90(key, sensor_len, sensor_cstr, scan_position,
                                          sensor_azimuth_angle, longitude, latitude, sample_lons[i],
                                          sample_lats[i], sample_power);
    EXPECT(oops::is_close_absolute(sample_power, reference_sample_powers[i], abs_tol));
  }

  fov::ufo_fov_delete_f90(key);
}

// -----------------------------------------------------------------------------

// Read config, set up params, call test helper function
void testFieldOfViewFromConf() {
  const eckit::Configuration& conf = ::test::TestEnvironment::config();
  const eckit::LocalConfiguration fovConf = conf.getSubConfiguration("field of view");
  FovWrapperTestParameters params{};
  params.validateAndDeserialize(fovConf);

  fovWrapperTestHelper(params.sensor, params.satellite, params.scan_position,
                       params.sensor_azimuth_angle, params.longitude, params.latitude,
                       params.reference_ellipse_longitudes, params.reference_ellipse_latitudes,
                       params.sample_longitudes, params.sample_latitudes,
                       params.reference_sample_powers, params.abs_tol);
}

// -----------------------------------------------------------------------------

class FieldOfView : public oops::Test {
 public:
  FieldOfView() = default;
  ~FieldOfView() = default;

 private:
  std::string testid() const override { return "ufo::test::FieldOfView"; }

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/FieldOfView/testFieldOfView") { testFieldOfViewFromConf(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo
