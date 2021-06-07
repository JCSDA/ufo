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
    GeoVaLs gv_temp(ospace.distribution(), ingeovars);

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
        gv_one.getAtLevel(gv_val, var, 1);
        EXPECT(oops::is_close_absolute(gv_val[0], values[i], oneloctol, 0, verbosity));
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
          gv.getAtLevel(gvar, ingeovars[i], k+1);
          gval.getAtLevel(gvarref, ingeovars[i], kk);
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

    GeoVaLs gv1(ospace.distribution(), gv.getVars());
    GeoVaLs gv2(ospace.distribution(), gv.getVars());
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

/// \brief Tests GeoVaLs::allocate, GeoVals::put, GeoVaLs::get,
/// GeoVaLs::putAtLevel, GeoVaLs::getAtLevel,
/// GeoVaLs::putAtLocation and GeoVaLs::getAtLocation.
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

  /// Set all values for the first level of a 2D variable to an arbitrary number with the
  /// put method. Then check that the get method returns the expected values in each case.
  /// Do this for several data types.
  /// (1) doubles
  const double fillvalue_double = 3.01234567890123;
  oops::Log::test() << "Put(double) fill value: " << fillvalue_double << std::endl;
  const std::vector<double> refvalues_double(gval.nlocs(), fillvalue_double);
  gval.putAtLevel(refvalues_double, var2, 1);
  std::vector<double> testvalues_double(gval.nlocs(), 0);
  gval.get(testvalues_double, var2);
  oops::Log::test() << "Get(double) result: " << testvalues_double << std::endl;
  EXPECT_EQUAL(testvalues_double, refvalues_double);
  /// (2) floats
  const float fillvalue_float = 4.1f;
  oops::Log::test() << "Put(float) fill value: " << fillvalue_float << std::endl;
  const std::vector<float> refvalues_float(gval.nlocs(), fillvalue_float);
  gval.putAtLevel(refvalues_float, var2, 1);
  std::vector<float> testvalues_float(gval.nlocs(), 0);
  gval.get(testvalues_float, var2);
  oops::Log::test() << "Get(float) result: " << testvalues_float << std::endl;
  EXPECT_EQUAL(testvalues_float, refvalues_float);
  /// (3) ints
  const int fillvalue_int = 5;
  oops::Log::test() << "Put(int) fill value: " << fillvalue_int << std::endl;
  const std::vector<int> refvalues_int(gval.nlocs(), fillvalue_int);
  gval.putAtLevel(refvalues_int, var2, 1);
  std::vector<int> testvalues_int(gval.nlocs(), 0);
  gval.get(testvalues_int, var2);
  oops::Log::test() << "Get(int) result: " << testvalues_int << std::endl;
  EXPECT_EQUAL(testvalues_int, refvalues_int);

  /// Check that the getAtLocation method returns what we put in the GeoVaLs.
  /// The reference GeoVaLs at indices (jlev, jloc) are equal to jlev + jloc.
  for (size_t jlev = 0; jlev < nlevs1; ++jlev) {
    std::vector<double> refvalues_loc_double(gval.nlocs());
    std::iota(refvalues_loc_double.begin(), refvalues_loc_double.end(), jlev);
    gval.putAtLevel(refvalues_loc_double, var1, jlev+1);
  }
  for (size_t jloc = 0; jloc < gval.nlocs(); ++jloc) {
    // Get the test vector at this location.
    std::vector<double> testvalues_loc_double(gval.nlevs(var1));
    gval.getAtLocation(testvalues_loc_double, var1, jloc);
    // Recreate reference vector for this location.
    std::vector<double> refvalues_loc_double(gval.nlevs(var1));
    std::iota(refvalues_loc_double.begin(), refvalues_loc_double.end(), jloc);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_loc_double, refvalues_loc_double);
    // Repeat the test for floats.
    std::vector<float> testvalues_loc_float(gval.nlevs(var1));
    gval.getAtLocation(testvalues_loc_float, var1, jloc);
    std::vector<float> refvalues_loc_float(refvalues_loc_double.begin(),
                                           refvalues_loc_double.end());
    EXPECT_EQUAL(testvalues_loc_float, refvalues_loc_float);
    // Repeat the test for ints.
    std::vector<int> testvalues_loc_int(gval.nlevs(var1));
    gval.getAtLocation(testvalues_loc_int, var1, jloc);
    std::vector<int> refvalues_loc_int(refvalues_loc_double.begin(),
                                       refvalues_loc_double.end());
    EXPECT_EQUAL(testvalues_loc_int, refvalues_loc_int);
  }

  /// Check that the putAtLocation method correctly puts values in the GeoVaLs.
  /// This is similar to the previous test but the putting and getting routines
  /// are transposed.
  /// This is done separately for each data type, using different fill values each time.
  /// (1) doubles
  /// The reference GeoVaLs at indices (jlev, jloc) are equal to jlev + jloc.
  oops::Log::test() << "putAtLoction with doubles" << std::endl;
  for (size_t jloc = 0; jloc < gval.nlocs(); ++jloc) {
    std::vector<double> refvalues_double(gval.nlevs(var1));
    std::iota(refvalues_double.begin(), refvalues_double.end(), jloc);
    gval.putAtLocation(refvalues_double, var1, jloc);
  }
  oops::Log::test() << "testing putAtLoction with doubles" << std::endl;
  for (size_t jlev = 0; jlev < gval.nlevs(var1); ++jlev) {
    // Get the test vector on this level.
    std::vector <double> testvalues_double(gval.nlocs());
    gval.getAtLevel(testvalues_double, var1, jlev + 1);
    // Recreate reference vector for this level.
    std::vector <double> refvalues_double(gval.nlocs());
    std::iota(refvalues_double.begin(), refvalues_double.end(), jlev);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_double, refvalues_double);
  }
  /// (2) floats
  /// The reference GeoVaLs at indices (jlev, jloc) are equal to jlev + jloc + 1.
  oops::Log::test() << "putAtLoction with floats" << std::endl;
  for (size_t jloc = 0; jloc < gval.nlocs(); ++jloc) {
    std::vector<float> refvalues_float(gval.nlevs(var1));
    std::iota(refvalues_float.begin(), refvalues_float.end(), jloc + 1);
    gval.putAtLocation(refvalues_float, var1, jloc);
  }
  oops::Log::test() << "testing putAtLoction with floats" << std::endl;
  for (size_t jlev = 0; jlev < gval.nlevs(var1); ++jlev) {
    // Get the test vector on this level.
    std::vector <float> testvalues_float(gval.nlocs());
    gval.getAtLevel(testvalues_float, var1, jlev + 1);
    // Recreate reference vector for this level.
    std::vector <float> refvalues_float(gval.nlocs());
    std::iota(refvalues_float.begin(), refvalues_float.end(), jlev + 1);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_float, refvalues_float);
  }
  /// (3) ints
  /// The reference GeoVaLs at indices (jlev, jloc) are equal to jlev + jloc + 2.
  oops::Log::test() << "putAtLoction with ints" << std::endl;
  for (size_t jloc = 0; jloc < gval.nlocs(); ++jloc) {
    std::vector<int> refvalues_int(gval.nlevs(var1));
    std::iota(refvalues_int.begin(), refvalues_int.end(), jloc + 2);
    gval.putAtLocation(refvalues_int, var1, jloc);
  }
  oops::Log::test() << "testing putAtLoction with ints" << std::endl;
  for (size_t jlev = 0; jlev < gval.nlevs(var1); ++jlev) {
    // Get the test vector on this level.
    std::vector <int> testvalues_int(gval.nlocs());
    gval.getAtLevel(testvalues_int, var1, jlev + 1);
    // Recreate reference vector for this level.
    std::vector <int> refvalues_int(gval.nlocs());
    std::iota(refvalues_int.begin(), refvalues_int.end(), jlev + 2);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_int, refvalues_int);
  }

  /// Check code paths that throw exceptions for the getAtLocation method.
  std::vector<double> testvalues_loc_wrongsize(gval.nlevs(var1) + 1, 0.0);
  EXPECT_THROWS(gval.getAtLocation(testvalues_loc_wrongsize, var1, 1));
  std::vector<double> testvalues_loc(gval.nlevs(var1), 0.0);
  EXPECT_THROWS(gval.getAtLocation(testvalues_loc, var1, -1));
  EXPECT_THROWS(gval.getAtLocation(testvalues_loc, var1, gval.nlocs()));

  /// Check code paths that throw exceptions for the putAtLocation method.
  EXPECT_THROWS(gval.putAtLocation(testvalues_loc_wrongsize, var1, 1));
  EXPECT_THROWS(gval.putAtLocation(testvalues_loc, var1, -1));
  EXPECT_THROWS(gval.putAtLocation(testvalues_loc, var1, gval.nlocs()));

  /// test 3D put and get
  for (size_t jlev = 0; jlev < nlevs1; ++jlev) {
    const float fillvalue = 3.0*(jlev+1);
    const std::vector<double> refvalues(gval.nlocs(), fillvalue);
    gval.putAtLevel(refvalues, var1, jlev + 1);
    oops::Log::test() << jlev << " level: put fill value: " << fillvalue << std::endl;
    std::vector<double> testvalues(gval.nlocs(), 0);
    gval.getAtLevel(testvalues, var1, jlev + 1);
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
