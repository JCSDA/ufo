/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PIECEWISELINEARINTERPOLATION_H_
#define TEST_UFO_PIECEWISELINEARINTERPOLATION_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {
namespace test {

CASE("ufo/PiecewiseLinearInterpolation/atInterpolationPoints") {
  ufo::PiecewiseLinearInterpolation interp({-1.0, 1.0, 5.0}, {2.0, 4.0, 0.0});

  EXPECT_EQUAL(interp(-1.0), 2.0);
  EXPECT_EQUAL(interp(1.0), 4.0);
  EXPECT_EQUAL(interp(5.0), 0.0);
}

CASE("ufo/PiecewiseLinearInterpolation/betweenInterpolationPoints") {
  ufo::PiecewiseLinearInterpolation interp({-1.0, 1.0, 5.0}, {2.0, 4.0, 0.0});

  EXPECT_EQUAL(interp(0.0), 3.0);
  EXPECT_EQUAL(interp(2.0), 3.0);
}

CASE("ufo/PiecewiseLinearInterpolation/outsideInterpolationPoints") {
  ufo::PiecewiseLinearInterpolation interp({-1.0, 1.0, 5.0}, {2.0, 4.0, 0.0});

  EXPECT_EQUAL(interp(-10.0), 2.0);
  EXPECT_EQUAL(interp(10.0), 0.0);
}

CASE("ufo/PiecewiseLinearInterpolation/singleInterpolationPoint") {
  ufo::PiecewiseLinearInterpolation interp({-1.0}, {2.0});

  EXPECT_EQUAL(interp(-10.0), 2.0);
  EXPECT_EQUAL(interp(-1.0), 2.0);
  EXPECT_EQUAL(interp(10.0), 2.0);
}

CASE("ufo/PiecewiseLinearInterpolation/noInterpolationPoints") {
  EXPECT_THROWS(ufo::PiecewiseLinearInterpolation({}, {}));
}

CASE("ufo/PiecewiseLinearInterpolation/differentNumberOfAbscissasAndOrdinates") {
  EXPECT_THROWS(ufo::PiecewiseLinearInterpolation({1.0, 2.0}, {1.0, 2.0, 3.0}));
  EXPECT_THROWS(ufo::PiecewiseLinearInterpolation({1.0, 2.0, 3.0}, {1.0, 2.0}));
}

class PiecewiseLinearInterpolation : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::PiecewiseLinearInterpolation";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PIECEWISELINEARINTERPOLATION_H_
