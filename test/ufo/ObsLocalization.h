/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSLOCALIZATION_H_
#define TEST_UFO_OBSLOCALIZATION_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/geometry/Point3.h"
#include "eckit/testing/Test.h"

#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"

#include "test/TestEnvironment.h"

#include "ufo/instantiateObsLocFactory.h"
#include "ufo/obslocalization/ObsLocalization.h"

namespace ufo {
namespace test {

/// \brief Tests the point3/point3 computeLocalization method of ufo obs localizations.
///
/// For each obs localization configured under "obs localizations" in the YAML, the test:
/// 1. Constructs the localization object.
/// 2. Loops over "point pairs" specified in the config.
/// 3. For each pair, calls computeLocalization(p1, p2) and compares to a
/// reference value.
/// NOTE this does not test the geometry iterator / obsvector localization....
/// that wasn't Travis' job that day.
void testObsLocalizationPoint3() {
  typedef ioda::ObsIterator Iterator_;

  instantiateObsLocFactory<Iterator_>();

  const eckit::Configuration & conf = ::test::TestEnvironment::config();
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  ioda::ObsSpace obsspace(obsconf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  const std::vector<eckit::LocalConfiguration> locConfs =
      conf.getSubConfigurations("obs localizations");

  for (const auto & locConf : locConfs) {
    const eckit::LocalConfiguration locSubConf(locConf, "localization");
    const std::string name = locSubConf.getString("localization method");
    oops::Log::info() << "Testing obs localization point3/point3: " << name << std::endl;

    ObsLocalization<Iterator_> obsloc(locSubConf, obsspace);

    const std::vector<eckit::LocalConfiguration> pairConfs =
        locConf.getSubConfigurations("point pair tests");

    for (const auto & pairConf : pairConfs) {
      const std::vector<double> coords1 = pairConf.getDoubleVector("p1");
      const std::vector<double> coords2 = pairConf.getDoubleVector("p2");
      ASSERT(coords1.size() == 3);
      ASSERT(coords2.size() == 3);
      const eckit::geometry::Point3 p1(coords1[0], coords1[1], coords1[2]);
      const eckit::geometry::Point3 p2(coords2[0], coords2[1], coords2[2]);

      const double tol = pairConf.getDouble("tolerance");
      const double ref = pairConf.getDouble("reference value");
      const double result = obsloc.computeLocalization(p1, p2);

      oops::Log::info() << "  p1=" << p1 << " p2=" << p2
                        << " result=" << result << " ref=" << ref << std::endl;
      EXPECT(std::abs(result - ref) <= tol);
    }
  }
}

// -----------------------------------------------------------------------------

class ObsLocalization : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~ObsLocalization() = default;

 private:
  std::string testid() const override { return "ufo::test::ObsLocalization"; }

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsLocalization/testObsLocalizationPoint3")
      { testObsLocalizationPoint3(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSLOCALIZATION_H_
