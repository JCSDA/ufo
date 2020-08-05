/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_GEOVALS_H_
#define TEST_UFO_GEOVALS_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsOperator.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testGeoVaLs() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  std::vector<eckit::LocalConfiguration> confs;
  conf.get("geovals test", confs);
  for (size_t jconf = 0; jconf < confs.size(); ++jconf) {
/// Setup ObsSpace
    const eckit::LocalConfiguration obsconf(confs[jconf], "obs space");
    ioda::ObsSpace ospace(obsconf, oops::mpi::comm(), bgn, end);

/// Setup GeoVaLs
    const eckit::LocalConfiguration gconf(confs[jconf], "geovals");
    const oops::Variables ingeovars(gconf, "state variables");
    const GeoVaLs gval(gconf, ospace, ingeovars);

    const double tol = gconf.getDouble("tolerance");

/// Check that GeoVaLs default constructor works
     oops::Log::trace() <<
      "GeoVaLs default constructor - does not allocate fields" << std::endl;
    GeoVaLs gv_temp(ospace.comm());

/// Check that GeoVaLs constructor to create a GeoVaLs with one location works
    if (gconf.has("one location check")) {
      oops::Log::trace() << "Check that GeoVaLs constructor for one location works" << std::endl;
      const eckit::LocalConfiguration gconfone(gconf, "one location check");
      const std::string var = gconfone.getString("variable");
      const std::vector<int> ind = gconfone.getIntVector("indices");
      const std::vector<float> values = gconfone.getFloatVector("values");
      const float oneloctol = gconfone.getFloat("tolerance");

      // Loop over each location and test just the lowest level
      oops::TestVerbosity verbosity = oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE;
      for (std::size_t i = 0; i < ind.size(); ++i) {
        GeoVaLs gv_one(gval, ind[i]);
        std::vector<float> gv_val(1);
        gv_one.get(gv_val, var, 1);
        EXPECT(oops::is_close_absolute(gv_val[0], values[i], oneloctol, verbosity));
      }
    } else {
      oops::Log::trace() << "Test just the constructor for a one location GeoVaLs" << std::endl;
      int index = 0;
      GeoVaLs gv_one(gval, index);
    }

/// Check that GeoVaLs merge followed by a split gives back the original geovals
    oops::Log::trace() <<
      "GeoVaLs merge followed by a split gives back the original geovals" << std::endl;

    GeoVaLs gv(gval);

    double dp_gval = gval.dot_product_with(gval);
    oops::Log::debug()<< "initial gval dot product with itself " << dp_gval << std::endl;

    gv.zero();
    gv.merge(gval, gval);

    double dp_gv_merged = gv.dot_product_with(gv);
    oops::Log::debug()<< "gv dot product with itself after merge " << dp_gv_merged << std::endl;
    EXPECT(abs(dp_gv_merged - 2.0 * dp_gval)/dp_gv_merged < tol);

    GeoVaLs gv1(ospace.comm()), gv2(ospace.comm());
    gv.split(gv1, gv2);

    double dp_gv1_split = gv1.dot_product_with(gv1);
    double dp_gv2_split = gv2.dot_product_with(gv2);

    oops::Log::debug()<< "gv1 gv2 dot products with itself after split "
                      << dp_gv1_split << " " << dp_gv2_split << std::endl;

    EXPECT(gv1.rms() == gv2.rms());
    EXPECT(gv1.rms() == gval.rms());
    EXPECT(abs(dp_gv1_split - dp_gv2_split)/dp_gv1_split < tol);
    EXPECT(abs(dp_gv1_split - dp_gval)/dp_gv1_split < tol);

    oops::Log::trace() <<
      "GeoVaLs merge followed by a split test succeeded" << std::endl;

///  Check that  GeoVaLs & operator *= (const std::vector<float>);
    oops::Log::trace() <<
      "Check that GeoVaLs & operator *= (const std::vector<float>);" << std::endl;

    std::size_t nlocs = ospace.nlocs();
    {
      std::vector<float> tw(nlocs, 1.0f);
      gv1 *= tw;
      EXPECT(gv1.rms() == gval.rms());
    }
    {
      std::vector<float> tw(nlocs, 2.0f);
      gv1 *= tw;
      double rms1, rms2;
      rms1 = gv1.rms();
      rms2 = 2.0 * gval.rms();
      oops::Log::debug()<< "rms1, rms2 = " <<  rms1  << "  " << rms2 << std::endl;
      EXPECT(rms1 == rms2);
    }
    oops::Log::trace() <<
      "GeoVaLs & operator *= (const std::vector<float>); test succeeded" << std::endl;
  }
}
// -----------------------------------------------------------------------------

class GeoVaLs : public oops::Test {
 public:
  GeoVaLs() {}
  virtual ~GeoVaLs() {}
 private:
  std::string testid() const {return "ufo::test::GeoVaLs";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/GeoVaLs/testGeoVaLs")
      { testGeoVaLs(); });
  }
};

// ----------;

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_GEOVALS_H_
