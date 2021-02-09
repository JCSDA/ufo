/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_GEOVALS_H_
#define TEST_UFO_GEOVALS_H_

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
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
    ioda::ObsSpace ospace(obsconf, oops::mpi::world(), bgn, end, oops::mpi::myself());

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

    GeoVaLs gv(gval);
    if (gconf.has("reorderzdir check")) {
      const eckit::LocalConfiguration gconfchk(gconf, "reorderzdir check");
      const std::string flipto = gconfchk.getString("direction");
      std::string flipback = (flipto == "bottom2top") ? "top2bottom" : "bottom2top";
      std::size_t nobs = ospace.nlocs();
      gv.reorderzdir("air_pressure_levels", flipto);
      std::vector<float> gvar(nobs);
      std::vector<float> gvarref(nobs);
      float sum;
      for (size_t i = 0; i < ingeovars.size(); ++i) {
        size_t nlevs = gv.nlevs(ingeovars[i]);
        sum = 0;
        for (size_t k = 0; k < nlevs; ++k) {
          size_t kk = nlevs - k;
          gv.get(gvar, ingeovars[i], k+1);
          gval.get(gvarref, ingeovars[i], kk);
          for (size_t  j = 0; j < nobs; ++j) {
            gvar[j] = gvar[j] - gvarref[j];
            sum += sum + gvar[j];
          }
        }
      }
      gv.reorderzdir("air_pressure_levels", flipback);
      const double tol = gconfchk.getDouble("tolerance");
      EXPECT(abs(sum) < tol);
    }

/// Check that GeoVaLs merge followed by a split gives back the original geovals
    oops::Log::trace() <<
      "GeoVaLs merge followed by a split gives back the original geovals" << std::endl;

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

/// \brief Tests GeoVaLs::allocate, GeoVaLs::put and GeoVaLs::get method (for 2D variables)
void testGeoVaLsAllocatePutGet() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  const eckit::LocalConfiguration testconf(conf, "geovals get test");

  /// Test 2D variables
  const std::string var1 = "variable1";
  const std::string var2 = "variable2";
  const oops::Variables testvars({var1, var2});

  /// Setup GeoVaLs that are not filled in or allocated; test that they are not allocated
  const Locations locs(testconf, oops::mpi::world());
  GeoVaLs gval(locs, testvars);
  oops::Log::test() << "GeoVals allocate test: created empty GeoVaLs with " << testvars <<
                       " variables and nlocs=" << gval.nlocs() << std::endl;
  EXPECT_EQUAL(gval.nlevs(var1), 0);
  EXPECT_EQUAL(gval.nlevs(var2), 0);

  /// Allocate only the first variable, and test that it's allocated correctly
  const oops::Variables testvar1({var1});
  const int nlevs1 = 10;
  gval.allocate(nlevs1, testvar1);
  oops::Log::test() << "Allocated " << var1 << "; nlevs(var1) = " << gval.nlevs(var1) <<
                                               "; nlevs(var2) = " << gval.nlevs(var2) << std::endl;
  EXPECT_EQUAL(gval.nlevs(var1), nlevs1);
  EXPECT_EQUAL(gval.nlevs(var2), 0);

  /// Allocate only the second variable (2D variable), test that it's allocated correctly
  const oops::Variables testvar2({var2});
  const int nlevs2 = 1;
  gval.allocate(nlevs2, testvar2);
  oops::Log::test() << "Allocated " << var2 << "; nlevs(var1) = " << gval.nlevs(var1) <<
                                               "; nlevs(var2) = " << gval.nlevs(var2) << std::endl;
  EXPECT_EQUAL(gval.nlevs(var1), nlevs1);
  EXPECT_EQUAL(gval.nlevs(var2), nlevs2);

  /// set all values for the first level of 2D variable to an arbitrary number
  const float fillvalue = 3.0;
  oops::Log::test() << "Put fill value: " << fillvalue << std::endl;
  const std::vector<double> refvalues(gval.nlocs(), fillvalue);
  const std::vector<float>  refvalues_float(refvalues.begin(), refvalues.end());
  const std::vector<int>    refvalues_int(refvalues.begin(), refvalues.end());
  gval.put(refvalues, var2, 1);

  /// check that get method returns what we put in geovals
  std::vector<double> testvalues(gval.nlocs(), 0);
  gval.get(testvalues, var2);
  oops::Log::test() << "Get result:        " << testvalues << std::endl;
  EXPECT_EQUAL(testvalues, refvalues);
  std::vector<float> testvalues_float(gval.nlocs(), 0);
  gval.get(testvalues_float, var2);
  oops::Log::test() << "Get(float) result: " << testvalues_float << std::endl;
  EXPECT_EQUAL(testvalues_float, refvalues_float);
  std::vector<int> testvalues_int(gval.nlocs(), 0);
  gval.get(testvalues_int, var2);
  oops::Log::test() << "Get(int) result:   " << testvalues_int << std::endl;
  EXPECT_EQUAL(testvalues_int, refvalues_int);

  /// test 3D put and get
  for (size_t jlev = 0; jlev < nlevs1; ++jlev) {
    const float fillvalue = 3.0*(jlev+1);
    const std::vector<double> refvalues(gval.nlocs(), fillvalue);
    gval.put(refvalues, var1, jlev+1);
    oops::Log::test() << jlev << " level: put fill value: " << fillvalue << std::endl;
    std::vector<double> testvalues(gval.nlocs(), 0);
    gval.get(testvalues, var1, jlev+1);
    oops::Log::test() << jlev << " level: get result: " << testvalues << std::endl;
    EXPECT_EQUAL(testvalues, refvalues);
  }
}

// -----------------------------------------------------------------------------

class GeoVaLs : public oops::Test {
 public:
  GeoVaLs() = default;
  virtual ~GeoVaLs() = default;
 private:
  std::string testid() const override {return "ufo::test::GeoVaLs";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/GeoVaLs/testGeoVaLs")
      { testGeoVaLs(); });
    ts.emplace_back(CASE("ufo/GeoVaLs/testGeoVaLsAllocatePutGet")
      { testGeoVaLsAllocatePutGet(); });
  }

  void clear() const override {}
};

// ----------;

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_GEOVALS_H_
