/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_SUPEROBRADARTEMPLATE_H_
#define TEST_UFO_SUPEROBRADARTEMPLATE_H_

#include <Eigen/Dense>

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "test/TestEnvironment.h"

#include "ufo/superob/SuperObRadarTemplate.h"

namespace ufo {
namespace test {

CASE("ufo/SuperObRadarTemplate/typicalExample") {
  // A typical example of how this template may be set up.

  ufo::SuperObRadarTemplateData SOtemplateData;
  SOtemplateData.stationLatitude = 0.0;
  SOtemplateData.stationLongitude = 0.0;
  SOtemplateData.stationElevation = 0.0;
  SOtemplateData.beamTiltAngle = 6.0;
  SOtemplateData.minGateRange = 0.0;
  SOtemplateData.numBeams = 10;
  SOtemplateData.numGates = 10;
  SOtemplateData.beamWidth = 36.0;
  SOtemplateData.gateWidth = 100.0;

  const int numBeamsInSuperObRegion = 5;
  const double superObRegionRadialExtent = 500.0;

  const ufo::SuperObRadarTemplate SOtemplate(numBeamsInSuperObRegion,
                                             superObRegionRadialExtent,
                                             SOtemplateData);

  const auto mask = SOtemplate.getMask();
  const auto distance = SOtemplate.getDistance();

  // Reference values of the mask and distance arrays produced
  // by an independent python script.
  Eigen::ArrayXXi refMask(12, 13);
  refMask << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  Eigen::ArrayXXd refDistance(12, 13);
  refDistance << 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 238.00929918, 210.44013314, 198.90364474, 210.44013314, 238.00929918,
    238.00929918, 210.44013314, 198.90364474, 210.44013314, 238.00929918, 0., 0.,
    0., 247.28096846, 155.10574595, 99.45169997, 155.10574595, 247.28096846,
    247.28096846, 155.10574595, 99.45169997, 155.10574595, 247.28096846, 0., 0.,
    0., 0., 153.66161072, 0., 153.66161072, 0., 0., 153.66161072, 0.,
    153.66161072, 0., 0., 0.,
    0., 0., 207.23705773, 99.45145515, 207.23705773, 0., 0., 207.23705773,
    99.45145515, 207.23705773, 0., 0., 0.,
    0., 0., 0., 198.90266545, 0., 0., 0., 0., 198.90266545, 0., 0., 0., 0.,
    0., 0., 0., 198.90119613, 0., 0., 0., 0., 198.90119613, 0., 0., 0., 0.,
    0., 0., 0., 99.4504756, 0., 0., 0., 0., 99.4504756, 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., 0., 99.45023065, 0., 0., 0., 0., 99.45023065, 0., 0., 0., 0.,
    0., 0., 0., 198.90021632, 0., 0., 0., 0., 198.90021632, 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;

  for (int i = 0; i < mask.rows(); ++i) {
    for (int j = 0; j < mask.cols(); ++j) {
      EXPECT(mask(i, j) == refMask(i, j));
      EXPECT(oops::is_close_relative(distance(i, j), refDistance(i, j), 1e-10));
    }
  }
}

CASE("ufo/SuperObRadarTemplate/exception") {
  // Trigger exceptions by misconfiguring the template.

  ufo::SuperObRadarTemplateData SOtemplateData;
  SOtemplateData.stationLatitude = 0.0;
  SOtemplateData.stationLongitude = 0.0;
  SOtemplateData.stationElevation = 0.0;
  SOtemplateData.beamTiltAngle = 6.0;
  SOtemplateData.minGateRange = 0.0;
  SOtemplateData.numBeams = 10;
  SOtemplateData.numGates = 10;
  SOtemplateData.beamWidth = 36.0;
  SOtemplateData.gateWidth = 100.0;

  EXPECT_THROWS_MSG(ufo::SuperObRadarTemplate SOtemplateExcept1(0, 500.0, SOtemplateData),
                    "Invalid number of beams in the superob region");

  EXPECT_THROWS_MSG(ufo::SuperObRadarTemplate SOtemplateExcept2(1000, 500.0, SOtemplateData),
                    "Invalid number of beams in the superob region");

  EXPECT_THROWS_MSG(ufo::SuperObRadarTemplate SOtemplateExcept3(5, 0.0, SOtemplateData),
                    "Invalid superob region radial extent");

  EXPECT_THROWS_MSG(ufo::SuperObRadarTemplate SOtemplateExcept4(5, 1.0e9, SOtemplateData),
                    "Invalid superob region radial extent");
};

class SuperObRadarTemplate : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::SuperObRadarTemplate";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_SUPEROBRADARTEMPLATE_H_
