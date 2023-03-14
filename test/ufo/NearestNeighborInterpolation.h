/*
 * (C) Copyright 2023 NOAA, UCAR, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_NEARESTNEIGHBORINTERPOLATION_H_
#define TEST_UFO_NEARESTNEIGHBORINTERPOLATION_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/utils/NearestNeighborInterpolation.h"

namespace ufo {
namespace test {

CASE("ufo/NearestNeighborInterpolation/atInterpolationPoints") {
  ufo::NearestNeighborInterpolation interp({1.0, 2.5, 5.0}, {2.0, 4.0, 1.0});

  EXPECT_EQUAL(interp(1.0), 2.0);
  EXPECT_EQUAL(interp(2.5), 4.0);
  EXPECT_EQUAL(interp(5.0), 1.0);
}

CASE("ufo/NearestNeighborInterpolation/betweenInterpolationPoints") {
  ufo::NearestNeighborInterpolation interp({1.0, 2.5, 5.0}, {2.0, 4.0, 1.0});

  EXPECT_EQUAL(interp(1.5), 2.0);
  EXPECT_EQUAL(interp(2.0), 4.0);
  EXPECT_EQUAL(interp(4.0), 1.0);
}

CASE("ufo/NearestNeighborInterpolation/outsideInterpolationPoints") {
  ufo::NearestNeighborInterpolation interp({1.0, 2.5, 5.0}, {2.0, 4.0, 1.0});

  EXPECT_EQUAL(interp(0.0), 2.0);
  EXPECT_EQUAL(interp(10.0), 1.0);
}

CASE("ufo/NearestNeighborInterpolation/singleInterpolationPoint") {
  ufo::NearestNeighborInterpolation interp({-1.0}, {2.0});

  EXPECT_EQUAL(interp(-10.0), 2.0);
  EXPECT_EQUAL(interp(-1.0), 2.0);
  EXPECT_EQUAL(interp(10.0), 2.0);
}

CASE("ufo/NearestNeighborInterpolation/noInterpolationPoints") {
  EXPECT_THROWS(ufo::NearestNeighborInterpolation({}, {}));
}

CASE("ufo/NearestNeighborInterpolation/differentNumberOfAbscissasAndOrdinates") {
  EXPECT_THROWS(ufo::NearestNeighborInterpolation({1.0, 2.0}, {1.0, 2.0, 3.0}));
  EXPECT_THROWS(ufo::NearestNeighborInterpolation({1.0, 2.0, 3.0}, {1.0, 2.0}));
}

class NearestNeighborInterpolation : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::NearestNeighborInterpolation";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_NEARESTNEIGHBORINTERPOLATION_H_
