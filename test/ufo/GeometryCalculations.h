/*
 * (C) Crown copyright 2024-2026, Met Office
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
  EXPECT(oops::is_close_relative(result2, 20015114.352186374, 1e-14));

  // Distance around the equator.
  const double result3 = ufo::haversine(0.0, 0.0, 0.0, 180.0);
  EXPECT(oops::is_close_relative(result3, 20015114.352186374, 1e-14));

  // Distances between some arbitrarily-chosen points.
  const double result4 = ufo::haversine(-60.0, 70.0, 30.0, -80.0);
  EXPECT(oops::is_close_relative(result4, 16001218.712405037, 1e-14));
  const double result5 = ufo::haversine(45.0, -34.0, 23.0, 12.0);
  EXPECT(oops::is_close_relative(result5, 4808559.7532829279, 1e-14));
  const double result6 = ufo::haversine(100.0, 80.0, 100.001, 79.999);
  EXPECT(oops::is_close_relative(result6, 112.85926326662934, 1e-14));
}

CASE("ufo/GeometryCalculations/convertRangeAzimToLatLon") {
  // Tests of the range/azimuth to lat/lon conversion routine which returns a pair of values.
  const auto[lat1, lon1] = ufo::convertRangeAzimToLatLon(0.0, 90.0, 60.0, 30.0, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat1, 60.0, 1e-14) &&
         oops::is_close_relative(lon1, 30.0, 1e-14));
  const auto[lat2, lon2] = ufo::convertRangeAzimToLatLon(10000.0, 90.0, 60.0, 30.0, 1000.0, 0.0);
  EXPECT(oops::is_close_relative(lat2, 60.0, 1e-14) &&
         oops::is_close_relative(lon2, 30.1798428192375, 1e-14));
  const auto[lat3, lon3] = ufo::convertRangeAzimToLatLon(10000.0, 90.0, 60.0, 30.0, 1000.0, 45.0);
  EXPECT(oops::is_close_relative(lat3, 60.0, 1e-14) &&
         oops::is_close_relative(lon3, 30.127062350739433, 1e-14));
  const auto[lat4, lon4] = ufo::convertRangeAzimToLatLon(10000.0, 90.0, 60.0, 30.0, 1000.0, 90.0);
  EXPECT(oops::is_close_relative(lat4, 60.0, 1e-14) &&
         oops::is_close_relative(lon4, 30.0, 1e-14));
  const auto[lat5, lon5] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, 0.0, 0.0, 0.0);
  // The station longitude (0) is valid in both conventions, so the longitude closest to the
  // station is returned. The longitude is wrapped via fmod, hence the looser tolerance.
  EXPECT(oops::is_close_relative(lat5, 0.06359152367493523, 1e-14) &&
         oops::is_close_relative(lon5, -0.063591523674915607, 1e-12));
  const auto[lat6, lon6] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 89.9, 0.0, 10000.0, -45.0);
  EXPECT(oops::is_close_relative(lat6, 89.944950509229699, 1e-14) &&
         oops::is_close_relative(lon6, 334.24524226617109 - 360.0, 1e-12));

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

CASE("ufo/GeometryCalculations/convertRangeAzimToLatLonLongitudeConvention") {
  // The same station expressed as 0 and 360 degrees. The [0, 360) convention is preserved
  // for the latter rather than the result being wrapped into the 350s or below zero.
  const auto[lat1, lon1] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, 0.0, 0.0, 0.0);
  const auto[lat1b, lon1b] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, 360.0, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat1b, lat1, 1e-14) &&
         oops::is_close_relative(lon1b, 359.93640847632508, 1e-14));

  // The same station expressed as -10 and 350 degrees: each result stays in the convention
  // of its station longitude and the two differ by exactly 360 degrees.
  const auto[lat2, lon2] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, -10.0, 0.0, 0.0);
  const auto[lat2b, lon2b] = ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, 350.0, 0.0, 0.0);
  EXPECT(lon2 < 0.0 && lon2b > 180.0);
  EXPECT(oops::is_close_relative(lat2, lat2b, 1e-14) &&
         oops::is_close_relative(lon2, lon2b - 360.0, 1e-12));

  // A station just west of the +180/-180 seam, offset mostly eastwards so that the result
  // crosses the seam. The station longitude is valid in both conventions, so the result is
  // the value closest to the station, i.e. just above 180.
  const auto[lat3, lon3] = ufo::convertRangeAzimToLatLon(10000.0, 89.5, 0.0, 179.95, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat3, 0.0007847947488409176, 1e-14) &&
         oops::is_close_relative(lon3, 180.03992857089824, 1e-14));

  // The same station expressed in the [-180, 180) convention: the result crosses the seam
  // and is wrapped to just above -180.
  const auto[lat4, lon4] = ufo::convertRangeAzimToLatLon(10000.0, 89.5, 0.0, -180.05, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat4, lat3, 1e-14) &&
         oops::is_close_relative(lon4, lon3 - 360.0, 1e-14));

  // The mirror image: a station just east of the seam in the [-180, 180) convention,
  // offset mostly westwards, is wrapped to just below +180.
  const auto[lat5, lon5] = ufo::convertRangeAzimToLatLon(10000.0, -89.5, 0.0, -179.95, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat5, lat3, 1e-14) &&
         oops::is_close_relative(lon5, 360.0 - lon3, 1e-14));

  // A station just west of 360 in the [0, 360) convention, offset eastwards, is wrapped
  // to just above 0.
  const auto[lat6, lon6] = ufo::convertRangeAzimToLatLon(10000.0, 89.5, 0.0, 359.99, 0.0, 0.0);
  EXPECT(oops::is_close_relative(lat6, lat3, 1e-14) &&
         oops::is_close_relative(lon6, lon3 - 179.95 + 359.99 - 360.0, 1e-10));

  // The in-place overload should follow the same convention.
  double lat2c, lon2c;
  ufo::convertRangeAzimToLatLon(10000.0, -45.0, 0.0, -10.0, 0.0, 0.0, lat2c, lon2c);
  EXPECT(oops::is_close_relative(lat2, lat2c, 1e-14) &&
         oops::is_close_relative(lon2, lon2c, 1e-14));
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
