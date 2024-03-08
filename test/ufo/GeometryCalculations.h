/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_GEOMETRYCALCULATIONS_H_
#define TEST_UFO_GEOMETRYCALCULATIONS_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "test/TestEnvironment.h"

#include "ufo/utils/GeometryCalculations.h"

namespace ufo {
namespace test {

CASE("ufo/GeometryCalculations/haversine") {
  // Tests of the haversine function.
  // Reference values have been obtained from an independent python script.

  // Distance between two points located at the South pole.
  const double result1 = ufo::haversine(-90.0, 0.0, -90.0, 0.0);
  EXPECT(oops::is_close_relative(result1, 0.0, 1e-14));

  // Distance between the South pole and the North pole.
  const double result2 = ufo::haversine(-90.0, 0.0, 90.0, 90.0);
  EXPECT(oops::is_close_relative(result2, 20015086.79602057, 1e-14));

  // Distance around the equator.
  const double result3 = ufo::haversine(0.0, 0.0, 0.0, 180.0);
  EXPECT(oops::is_close_relative(result3, 20015086.79602057, 1e-14));

  // Distances between some arbitrarily-chosen points.
  const double result4 = ufo::haversine(-60.0, 70.0, 30.0, -80.0);
  EXPECT(oops::is_close_relative(result4, 16001196.682441674, 1e-14));
  const double result5 = ufo::haversine(45.0, -34.0, 23.0, 12.0);
  EXPECT(oops::is_close_relative(result5, 4808553.133012492, 1e-14));
  const double result6 = ufo::haversine(100.0, 80.0, 100.001, 79.999);
  EXPECT(oops::is_close_relative(result6, 112.85910788562497, 1e-14));
}

CASE("ufo/GeometryCalculations/convertRangeAzimToLatLon") {
  // Tests of the range/azimuth to lat/lon conversion routine which returns a pair of values.
  const auto[lat1, lon1] = ufo::convertRangeAzimToLatLon(0.0, 90.0, 60.0, 30.0, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat1, 60.0, 1e-14) &&
         oops::is_close_relative(lon1, 30.0, 1e-14));
  const auto[lat2, lon2] = ufo::convertRangeAzimToLatLon(10000.0, 90.0, 60.0, 30.0, 1000.0, 0.0);
  EXPECT(oops::is_close_relative(lat2, 60.0, 1e-14) &&
         oops::is_close_relative(lon2, 30.179843066810278, 1e-14));
  const auto[lat3, lon3] = ufo::convertRangeAzimToLatLon(10000.0, 90.0, 60.0, 30.0, 1000.0, 45.0);
  EXPECT(oops::is_close_relative(lat3, 60.0, 1e-14) &&
         oops::is_close_relative(lon3, 30.127062525508897, 1e-14));
  const auto[lat4, lon4] = ufo::convertRangeAzimToLatLon(10000.0, 90.0, 60.0, 30.0, 1000.0, 90.0);
  EXPECT(oops::is_close_relative(lat4, 60.0, 1e-14) &&
         oops::is_close_relative(lon4, 30.0, 1e-14));
  const auto[lat5, lon5] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, 0.0, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat5, 0.06359161122573972, 1e-14) &&
         oops::is_close_relative(lon5, 359.9364083887743, 1e-14));
  const auto[lat6, lon6] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 89.9, 0.0, 10000.0, -45.0);
  EXPECT(oops::is_close_relative(lat6, 89.94495057109484, 1e-14) &&
         oops::is_close_relative(lon6, 334.2452068200381, 1e-14));

  // Tests of the range/azimuth to lat/lon conversion routine which modifies values in place.
  double lat1b, lon1b;
  ufo::convertRangeAzimToLatLon(0.0, 90.0, 60.0, 30.0, 0.0, 0.0, lat1b, lon1b);
  EXPECT(oops::is_close_relative(lat1, lat1b, 1e-14) &&
         oops::is_close_relative(lon1, lon1b, 1e-14));

  // Expect an exception to be thrown if the latitude is too close to either of the poles.
  EXPECT_THROWS_MSG(ufo::convertRangeAzimToLatLon(0.0, 0.0, 89.999, 0.0, 0.0, 0.0),
                    "Absolute value of station latitude cannot be greater than 89.9 degrees");
  EXPECT_THROWS_MSG(ufo::convertRangeAzimToLatLon(0.0, 0.0, -89.999, 0.0, 0.0, 0.0),
                    "Absolute value of station latitude cannot be greater than 89.9 degrees");
}

class GeometryCalculations : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::GeometryCalculations";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_GEOMETRYCALCULATIONS_H_
