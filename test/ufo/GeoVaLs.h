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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/distribution/RoundRobin.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Locations.h"
#include "oops/base/SamplingMethodSelector.h"
#include "oops/interface/SampledLocations.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsOperator.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace eckit
{
  // Don't use the contracted output for this type: the current implementation works only
  // with integer types.
  template <typename T>
  struct VectorPrintSelector<util::Range<T>> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testGeoVaLs() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  std::vector<eckit::LocalConfiguration> confs;
  conf.get("geovals test", confs);
  for (size_t jconf = 0; jconf < confs.size(); ++jconf) {
/// Setup ObsSpace
    const eckit::LocalConfiguration obsconf(confs[jconf], "obs space");
    ioda::ObsSpace ospace(obsconf, oops::mpi::world(), timeWindow, oops::mpi::myself());

/// Setup GeoVaLs
    const eckit::LocalConfiguration geovalsconf(confs[jconf], "geovals");
    const eckit::LocalConfiguration gconf(confs[jconf], "geovals test");
    const oops::Variables ingeovars(gconf, "state variables");
    const GeoVaLs gval(geovalsconf, ospace, ingeovars);

    const double tol = gconf.getDouble("tolerance");

/// Check that GeoVaLs default constructor works
     oops::Log::trace() <<
      "GeoVaLs default constructor - does not allocate fields" << std::endl;
    GeoVaLs gv_temp(ospace.distribution(), ingeovars);

/// Check that GeoVaLs constructor to create a GeoVaLs with one location works
    if (gconf.has("one location check")) {
      oops::Log::trace() << "Check that GeoVaLs constructor for one location works" << std::endl;
      const eckit::LocalConfiguration gconfone(gconf, "one location check");
      const oops::Variable var{gconfone.getString("variable")};
      const std::vector<int> ind = gconfone.getIntVector("indices");
      const std::vector<float> values = gconfone.getFloatVector("values");
      const float oneloctol = gconfone.getFloat("tolerance");

      // Loop over each location and test just the lowest level
      oops::TestVerbosity verbosity = oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE;
      for (std::size_t i = 0; i < ind.size(); ++i) {
        GeoVaLs gv_one(gval, ind[i]);
        std::vector<float> gv_val(1);
        gv_one.getAtLevel(gv_val, var, 0);
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
          size_t kk = nlevs - k - 1;
          gv.getAtLevel(gvar, ingeovars[i], k);
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

///  Check that  GeoVaLs & operator *= (const std::vector<float>);
    oops::Log::trace() <<
      "Check that GeoVaLs & operator *= (const std::vector<float>);" << std::endl;

    GeoVaLs gv1(gval);
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

enum class GetMethod {
  GET_PROFILE,
  GET_AT_LOCATION
};

/// Helper function checking that the getProfile or getAtLocation method returns what we put in the
/// GeoVaLs.
/// The reference GeoVaLs at indices (jlev, jprofile) are equal to seed + jlev + jprofile.
void testPutAtLevelAndGetProfileOrGetAtLocation(
    GeoVaLs &gval, const oops::Variable &var, double seed, GetMethod getMethodToTest) {
  const size_t nprofiles = gval.nprofiles(var);
  const size_t nlevs = gval.nlevs(var);
  for (size_t jlev = 0; jlev < nlevs; ++jlev) {
    std::vector<double> refvalues_path_double(nprofiles);
    std::iota(refvalues_path_double.begin(), refvalues_path_double.end(), seed + jlev);
    gval.putAtLevel(refvalues_path_double, var, jlev);
  }
  for (size_t jprofile = 0; jprofile < nprofiles; ++jprofile) {
    // Get the test vector at this location.
    std::vector<double> testvalues_path_double(nlevs);
    if (getMethodToTest == GetMethod::GET_PROFILE)
      gval.getProfile(testvalues_path_double, var, jprofile);
    else
      gval.getAtLocation(testvalues_path_double, var, jprofile);
    // Recreate reference vector for this location.
    std::vector<double> refvalues_path_double(nlevs);
    std::iota(refvalues_path_double.begin(), refvalues_path_double.end(), seed + jprofile);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_path_double, refvalues_path_double);

    // Repeat the test for floats.
    std::vector<float> testvalues_path_float(nlevs);
    if (getMethodToTest == GetMethod::GET_PROFILE)
      gval.getProfile(testvalues_path_float, var, jprofile);
    else
      gval.getAtLocation(testvalues_path_float, var, jprofile);
    std::vector<float> refvalues_path_float(refvalues_path_double.begin(),
                                           refvalues_path_double.end());
    EXPECT_EQUAL(testvalues_path_float, refvalues_path_float);

    // Repeat the test for ints.
    std::vector<int> testvalues_path_int(nlevs);
    if (getMethodToTest == GetMethod::GET_PROFILE)
      gval.getProfile(testvalues_path_int, var, jprofile);
    else
      gval.getAtLocation(testvalues_path_int, var, jprofile);
    std::vector<int> refvalues_path_int(refvalues_path_double.begin(),
                                       refvalues_path_double.end());
    EXPECT_EQUAL(testvalues_path_int, refvalues_path_int);
  }
}

enum class PutMethod {
  PUT_PROFILE,
  PUT_AT_LOCATION
};

/// Helper function checking that the putProfile or putAtLocation method puts correct values in
/// the GeoVaLs. This is similar to the function above, but the putting and getting routines
/// are transposed.
/// This is done separately for each data type, using different fill values each time.
void testPutProfileOrPutAtLocationAndGetAtLevel(
    GeoVaLs &gval, const oops::Variable &var, double seed, PutMethod putMethodToTest) {
  const size_t nprofiles = gval.nprofiles(var);
  const size_t nlevs = gval.nlevs(var);

  /// (1) doubles
  /// The reference GeoVaLs at indices (jlev, jprofile) are equal to seed + jlev + jprofile.
  oops::Log::test() << "putProfile/putAtLocation with doubles" << std::endl;
  for (size_t jprofile = 0; jprofile < nprofiles; ++jprofile) {
    std::vector<double> refvalues_double(nlevs);
    std::iota(refvalues_double.begin(), refvalues_double.end(), seed + jprofile);
    if (putMethodToTest == PutMethod::PUT_PROFILE)
      gval.putProfile(refvalues_double, var, jprofile);
    else
      gval.putAtLocation(refvalues_double, var, jprofile);
  }
  oops::Log::test() << "testing putProfile/putAtLocation with doubles" << std::endl;
  for (size_t jlev = 0; jlev < nlevs; ++jlev) {
    // Get the test vector on this level.
    std::vector <double> testvalues_double(nprofiles);
    gval.getAtLevel(testvalues_double, var, jlev);
    // Recreate reference vector for this level.
    std::vector <double> refvalues_double(nprofiles);
    std::iota(refvalues_double.begin(), refvalues_double.end(), seed + jlev);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_double, refvalues_double);
  }
  /// (2) floats
  /// The reference GeoVaLs at indices (jlev, jprofile) are equal to seed + jlev + jprofile + 1.
  oops::Log::test() << "putProfile/putAtLocation with floats" << std::endl;
  for (size_t jprofile = 0; jprofile < nprofiles; ++jprofile) {
    std::vector<float> refvalues_float(nlevs);
    std::iota(refvalues_float.begin(), refvalues_float.end(), seed + jprofile + 1);
    if (putMethodToTest == PutMethod::PUT_PROFILE)
      gval.putProfile(refvalues_float, var, jprofile);
    else
      gval.putAtLocation(refvalues_float, var, jprofile);
  }
  oops::Log::test() << "testing putProfile/putAtLocation with floats" << std::endl;
  for (size_t jlev = 0; jlev < nlevs; ++jlev) {
    // Get the test vector on this level.
    std::vector <float> testvalues_float(nprofiles);
    gval.getAtLevel(testvalues_float, var, jlev);
    // Recreate reference vector for this level.
    std::vector <float> refvalues_float(nprofiles);
    std::iota(refvalues_float.begin(), refvalues_float.end(), seed + jlev + 1);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_float, refvalues_float);
  }
  /// (3) ints
  /// The reference GeoVaLs at indices (jlev, jprofile) are equal to seed + jlev + jprofile + 2.
  oops::Log::test() << "putProfile/putAtLocation with ints" << std::endl;
  for (size_t jprofile = 0; jprofile < nprofiles; ++jprofile) {
    std::vector<int> refvalues_int(nlevs);
    std::iota(refvalues_int.begin(), refvalues_int.end(), seed + jprofile + 2);
    if (putMethodToTest == PutMethod::PUT_PROFILE)
      gval.putProfile(refvalues_int, var, jprofile);
    else
      gval.putAtLocation(refvalues_int, var, jprofile);
  }
  oops::Log::test() << "testing putProfile/putAtLocation with ints" << std::endl;
  for (size_t jlev = 0; jlev < nlevs; ++jlev) {
    // Get the test vector on this level.
    std::vector <int> testvalues_int(nprofiles);
    gval.getAtLevel(testvalues_int, var, jlev);
    // Recreate reference vector for this level.
    std::vector <int> refvalues_int(nprofiles);
    std::iota(refvalues_int.begin(), refvalues_int.end(), seed + jlev + 2);
    // Compare the two vectors.
    EXPECT_EQUAL(testvalues_int, refvalues_int);
  }
}

/// \brief Tests GeoVaLs::allocate, GeoVals::put, GeoVaLs::get,
/// GeoVaLs::putAtLevel, GeoVaLs::getAtLevel,
/// GeoVaLs::putProfile and GeoVaLs::getProfile.
void testGeoVaLsAllocatePutGet() {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  const eckit::LocalConfiguration testconf(conf, "geovals get test");

  /// Test 2D variables
  const oops::Variable var1{"variable1"};
  const oops::Variable var2{"variable2"};
  const oops::Variables testvars({var1, var2});

  /// Setup GeoVaLs that are not filled in or allocated; test that they are not allocated
  SampledLocations_ sampledLocations(testconf, oops::mpi::world());
  const size_t nprofiles = sampledLocations.sampledLocations().size();
  GeoVaLs gval(std::move(sampledLocations), testvars);
  oops::Log::test() << "GeoVals allocate test: created empty GeoVaLs with " << testvars <<
                       " variables " << std::endl;
  EXPECT_EQUAL(gval.nlevs(var1), 0);
  EXPECT_EQUAL(gval.nlevs(var2), 0);
  std::cout << "gval.nprofiles(var1): " << gval.nprofiles(var1) << std::endl;

  EXPECT_EQUAL(gval.nprofiles(var1), nprofiles);
  EXPECT_EQUAL(gval.nprofiles(var2), nprofiles);

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
  const std::vector<double> refvalues_double(nprofiles, fillvalue_double);
  gval.putAtLevel(refvalues_double, var2, 0);
  std::vector<double> testvalues_double(nprofiles, 0);
  gval.get(testvalues_double, var2);
  oops::Log::test() << "Get(double) result: " << testvalues_double << std::endl;
  EXPECT_EQUAL(testvalues_double, refvalues_double);
  /// (2) floats
  const float fillvalue_float = 4.1f;
  oops::Log::test() << "Put(float) fill value: " << fillvalue_float << std::endl;
  const std::vector<float> refvalues_float(nprofiles, fillvalue_float);
  gval.putAtLevel(refvalues_float, var2, 0);
  std::vector<float> testvalues_float(nprofiles, 0);
  gval.get(testvalues_float, var2);
  oops::Log::test() << "Get(float) result: " << testvalues_float << std::endl;
  EXPECT_EQUAL(testvalues_float, refvalues_float);
  /// (3) ints
  const int fillvalue_int = 5;
  oops::Log::test() << "Put(int) fill value: " << fillvalue_int << std::endl;
  const std::vector<int> refvalues_int(nprofiles, fillvalue_int);
  gval.putAtLevel(refvalues_int, var2, 0);
  std::vector<int> testvalues_int(nprofiles, 0);
  gval.get(testvalues_int, var2);
  oops::Log::test() << "Get(int) result: " << testvalues_int << std::endl;
  EXPECT_EQUAL(testvalues_int, refvalues_int);

  // Test the putAtLevel and getProfile methods.
  testPutAtLevelAndGetProfileOrGetAtLocation(gval, var1, 1.3 /*seed*/, GetMethod::GET_PROFILE);
  // When each location is sampled only by a single interpolation path (as in this test case),
  // getProfile and getAtLocation are equivalent, so let's test the latter too.
  testPutAtLevelAndGetProfileOrGetAtLocation(gval, var1, 2.7 /*seed*/, GetMethod::GET_AT_LOCATION);

  // Test the putProfile and getAtLevel methods.
  testPutProfileOrPutAtLocationAndGetAtLevel(gval, var1, 1.3 /*seed*/, PutMethod::PUT_PROFILE);
  // Again, in this test case putProfile and putAtLocation are equivalent, so let's test the latter
  // too.
  testPutProfileOrPutAtLocationAndGetAtLevel(gval, var1, 2.7 /*seed*/, PutMethod::PUT_AT_LOCATION);

  /// Check code paths that throw exceptions for the getProfile method.
  std::vector<double> testvalues_path_wrongsize(gval.nlevs(var1) + 1, 0.0);
  EXPECT_THROWS(gval.getProfile(testvalues_path_wrongsize, var1, 1));
  std::vector<double> testvalues_path(gval.nlevs(var1), 0.0);
  EXPECT_THROWS(gval.getProfile(testvalues_path, var1, -1));
  EXPECT_THROWS(gval.getProfile(testvalues_path, var1, nprofiles));

  /// Check code paths that throw exceptions for the putProfile method.
  EXPECT_THROWS(gval.putProfile(testvalues_path_wrongsize, var1, 1));
  EXPECT_THROWS(gval.putProfile(testvalues_path, var1, -1));
  EXPECT_THROWS(gval.putProfile(testvalues_path, var1, nprofiles));

  /// test 3D put and get
  for (size_t jlev = 0; jlev < nlevs1; ++jlev) {
    const float fillvalue = 3.0*(jlev+1);
    const std::vector<double> refvalues(nprofiles, fillvalue);
    gval.putAtLevel(refvalues, var1, jlev);
    oops::Log::test() << jlev << " level: put fill value: " << fillvalue << std::endl;
    std::vector<double> testvalues(nprofiles, 0);
    gval.getAtLevel(testvalues, var1, jlev);
    oops::Log::test() << jlev << " level: get result: " << testvalues << std::endl;
    EXPECT_EQUAL(testvalues, refvalues);
  }
}

/// \brief Tests GeoVaLs(const Locations &, const Variables &, const std::vector<size_t> &)
/// constructor. Tests that levels get correctly allocated, and that the GeoVaLs are zeroed out.
void testGeoVaLsConstructor() {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  const eckit::LocalConfiguration testconf(conf, "geovals get test");

  const oops::Variable var1{"variable1"};
  const oops::Variable var2{"variable2"};
  const oops::Variables testvars({var1, var2});

  /// Setup GeoVaLs; test that they are allocated
  SampledLocations_ sampledLocations(testconf, oops::mpi::world());
  const size_t npaths = sampledLocations.sampledLocations().size();
  const std::vector<size_t> testnlevs({10, 1});
  GeoVaLs gval(std::move(sampledLocations), testvars, testnlevs);
  oops::Log::test() << "Created and allocated GeoVaLs: " << var1 << "; nlevs(var1) = " <<
                       gval.nlevs(var1) << ", " << var2 << "; nlevs(var2) = " <<
                       gval.nlevs(var2) << std::endl;
  EXPECT_EQUAL(gval.nlevs(var1), 10);
  EXPECT_EQUAL(gval.nlevs(var2), 1);
  EXPECT_EQUAL(gval.nprofiles(var1), npaths);
  EXPECT_EQUAL(gval.nprofiles(var2), npaths);

  const std::vector<double> refvalues(npaths, 0.0);
  const std::vector<float>  refvalues_float(npaths, 0.0);
  const std::vector<int>    refvalues_int(npaths, 0);

  /// check that get method returns zeroes (initial values in GeoVaLs)
  std::vector<double> testvalues(npaths, 1.0);
  gval.get(testvalues, var2);
  oops::Log::test() << "Get result:        " << testvalues << std::endl;
  EXPECT_EQUAL(testvalues, refvalues);
  std::vector<float> testvalues_float(npaths, 1.0);
  gval.get(testvalues_float, var2);
  oops::Log::test() << "Get(float) result: " << testvalues_float << std::endl;
  EXPECT_EQUAL(testvalues_float, refvalues_float);
  std::vector<int> testvalues_int(npaths, 1);
  gval.get(testvalues_int, var2);
  oops::Log::test() << "Get(int) result:   " << testvalues_int << std::endl;
  EXPECT_EQUAL(testvalues_int, refvalues_int);
}

// Classes and functions needed by the following test

class GeoVaLParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeoVaLParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> name{"variable", this};
  oops::RequiredParameter<size_t> nlevels{"nlevels", this};
  oops::RequiredParameter<size_t> samplingMethod{"sampling method", this};
};

class SamplingMethodParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SamplingMethodParameters, Parameters)
 public:
  oops::RequiredParameter<std::vector<float>> longitudes{"longitudes", this};
  oops::RequiredParameter<std::vector<float>> latitudes{"latitudes", this};
  oops::RequiredParameter<std::vector<util::DateTime>> datetimes{"datetimes", this};
  oops::RequiredParameter<std::vector<std::pair<size_t, size_t>>> pathsGroupedByLocation{
    "paths grouped by location", this};
};

class MultipleSamplingMethodsTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MultipleSamplingMethodsTestParameters, Parameters)
 public:
  oops::RequiredParameter<std::vector<SamplingMethodParameters>> samplingMethods{
      "sampling methods", this};
  oops::RequiredParameter<std::vector<GeoVaLParameters>> geovals{"geovals", this};
  oops::RequiredParameter<size_t> gnlocs{"gnlocs", this};
};

class Selector : public oops::SamplingMethodSelector {
 public:
  explicit Selector(std::map<oops::Variable, size_t> samplingMethodIndexByVar) :
    samplingMethodIndexByVar_(std::move(samplingMethodIndexByVar))
  {}

  size_t methodIndex(const oops::Variable &var) const override {
    return samplingMethodIndexByVar_.at(var);
  }

 private:
  std::map<oops::Variable, size_t> samplingMethodIndexByVar_;
};

std::shared_ptr<ioda::Distribution> makeDistribution(size_t gnlocs) {
  auto dist = std::make_shared<ioda::RoundRobin>(
        oops::mpi::world(), ioda::EmptyDistributionParameters());
  const eckit::geometry::Point2 point(0, 0);
  for (size_t loc = 0; loc < gnlocs; ++loc) {
    dist->assignRecord(loc, loc, point);
  }
  dist->computePatchLocs();
  return dist;
}

/// Create a GeoVaLs object holding variables interpolated along paths produced by multiple
/// location sampling methods.
/// Verify that the GeoVaLs have been allocated to the correct size and they are linked to the
/// correct sampling methods.
/// Verify that the put... and get... methods work correctly.
void testGeoVaLsWithMultipleSamplingMethods() {
  typedef oops::Locations<ObsTraits> Locations_;
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  const eckit::LocalConfiguration testconf(conf, "multiple sampling methods test");

  MultipleSamplingMethodsTestParameters params;
  params.validateAndDeserialize(testconf);

  std::shared_ptr<ioda::Distribution> dist = makeDistribution(params.gnlocs);

  oops::Variables variables;
  std::vector<size_t> nlevels;
  std::map<oops::Variable, size_t> samplingMethodIndexByVar;
  for (const GeoVaLParameters & geovalParams : params.geovals.value()) {
    variables.push_back(geovalParams.name);
    nlevels.push_back(geovalParams.nlevels);
    samplingMethodIndexByVar.emplace(oops::Variable{geovalParams.name},
                                     geovalParams.samplingMethod);
  }

  std::vector<SampledLocations_> sampledLocationsVec;
  for (const SamplingMethodParameters & samplingMethodParams : params.samplingMethods.value()) {
    std::vector<util::Range<size_t>> pathsGroupedByLocation;
    for (const std::pair<size_t, size_t> &range :
         samplingMethodParams.pathsGroupedByLocation.value())
      pathsGroupedByLocation.emplace_back(util::Range<size_t>{range.first, range.second});

    sampledLocationsVec.emplace_back(
          std::make_unique<SampledLocations>(samplingMethodParams.longitudes,
                                             samplingMethodParams.latitudes,
                                             samplingMethodParams.datetimes,
                                             dist,
                                             pathsGroupedByLocation));
  }

  // Reference values to compare against later
  std::map<oops::Variable, size_t> expectedNumProfilesByVar;
  std::map<oops::Variable, std::vector<util::Range<size_t>>> expectedGroupedProfileIndicesByVar;
  for (size_t ivar = 0; ivar < variables.size(); ++ivar) {
    const SampledLocations &sampledLocations =
        sampledLocationsVec[samplingMethodIndexByVar.at(variables[ivar])].sampledLocations();
    expectedNumProfilesByVar[variables[ivar]] = sampledLocations.size();
    expectedGroupedProfileIndicesByVar[variables[ivar]] = sampledLocations.pathsGroupedByLocation();
  }

  // Create an Locations object and pass it to a GeoVaLs constructor.
  Locations_ locations(std::move(sampledLocationsVec),
                       std::make_unique<Selector>(samplingMethodIndexByVar));
  GeoVaLs geovals(std::move(locations), variables, nlevels);

  // Test that GeoVaLs have the correct size and are linked to the correct sampling methods.
  for (size_t ivar = 0; ivar < variables.size(); ++ivar) {
    EXPECT_EQUAL(geovals.nlevs(variables[ivar]), nlevels[ivar]);
    EXPECT_EQUAL(geovals.nprofiles(variables[ivar]), expectedNumProfilesByVar.at(variables[ivar]));

    std::vector<util::Range<size_t>> groupedProfileIndices;
    geovals.getProfileIndicesGroupedByLocation(variables[ivar], groupedProfileIndices);
    const std::vector<util::Range<size_t>> &expectedGroupedProfileIndices =
        expectedGroupedProfileIndicesByVar.at(variables[ivar]);

    EXPECT_EQUAL(groupedProfileIndices, expectedGroupedProfileIndices);
  }

  for (size_t ivar = 0; ivar < variables.size(); ++ivar) {
    testPutAtLevelAndGetProfileOrGetAtLocation(geovals, variables[ivar], 1.3 /*seed*/,
                                               GetMethod::GET_PROFILE);
    testPutProfileOrPutAtLocationAndGetAtLevel(geovals, variables[ivar], 2.7 /*seed*/,
                                               PutMethod::PUT_PROFILE);
  }
}

oops::Locations<ObsTraits> createLocationsForFormatsTest() {
  typedef oops::Locations<ObsTraits> Locations_;
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  const size_t gnlocs = 4;
  const size_t nlocs = 4;
  std::shared_ptr<ioda::Distribution> dist = makeDistribution(gnlocs);

  std::vector<SampledLocations_> sampledLocationsVec;
  {  // sampling method # 0
    const size_t npaths = nlocs;
    std::vector<util::Range<size_t>> pathsGroupedByLocation{{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    // The actual longitudes, latitudes and times don't matter
    std::vector<float> longitudes(npaths, 0.0);
    std::vector<float> latitudes(npaths, 0.0);
    std::vector<util::DateTime> times(npaths, util::DateTime(2001, 1, 1, 0, 0, 0));
    sampledLocationsVec.emplace_back(
          std::make_unique<SampledLocations>(longitudes, latitudes, times, dist,
                                             pathsGroupedByLocation));
  }
  {  // sampling method # 1
    const size_t npaths = 2 * nlocs;
    std::vector<util::Range<size_t>> pathsGroupedByLocation{{0, 2}, {2, 4}, {4, 6}, {6, 8}};
    // The actual longitudes, latitudes and times don't matter
    std::vector<float> longitudes(npaths, 0.0);
    std::vector<float> latitudes(npaths, 0.0);
    std::vector<util::DateTime> times(npaths, util::DateTime(2001, 1, 1, 0, 0, 0));
    sampledLocationsVec.emplace_back(
          std::make_unique<SampledLocations>(longitudes, latitudes, times, dist,
                                             pathsGroupedByLocation));
  }

  std::map<oops::Variable, size_t> samplingMethodIndexByVar{
    {oops::Variable{"var1_sampled"}, 1},
    {oops::Variable{"var2_reduced"}, 0},
    {oops::Variable{"var3_sampled_aliased_with_reduced"}, 0},
    {oops::Variable{"var4_sampled_distinct_from_reduced"}, 1},
  };

  return Locations_(std::move(sampledLocationsVec),
                    std::make_unique<Selector>(samplingMethodIndexByVar));
}

void testSampledAndReducedFormats(const GeoVaLs &geovals) {
  EXPECT(geovals.has(oops::Variable{"var1_sampled"}, GeoVaLFormat::SAMPLED));
  EXPECT_NOT(geovals.has(oops::Variable{"var1_sampled"}, GeoVaLFormat::REDUCED));
  EXPECT_NOT(geovals.has(oops::Variable{"var2_reduced"}, GeoVaLFormat::SAMPLED));
  EXPECT(geovals.has(oops::Variable{"var2_reduced"}, GeoVaLFormat::REDUCED));
  EXPECT(geovals.has(oops::Variable{"var3_sampled_aliased_with_reduced"}, GeoVaLFormat::SAMPLED));
  EXPECT(geovals.has(oops::Variable{"var3_sampled_aliased_with_reduced"}, GeoVaLFormat::REDUCED));
  EXPECT(geovals.has(oops::Variable{"var4_sampled_distinct_from_reduced"}, GeoVaLFormat::SAMPLED));
  EXPECT(geovals.has(oops::Variable{"var4_sampled_distinct_from_reduced"}, GeoVaLFormat::REDUCED));

  EXPECT_NOT(geovals.areReducedAndSampledFormatsAliased(oops::Variable{"var1_sampled"}));
  EXPECT_NOT(geovals.areReducedAndSampledFormatsAliased(oops::Variable{"var2_reduced"}));
  EXPECT(geovals.areReducedAndSampledFormatsAliased(oops::Variable
                                                    {"var3_sampled_aliased_with_reduced"}));
  EXPECT_NOT(geovals.areReducedAndSampledFormatsAliased(oops::Variable
                                                        {"var4_sampled_distinct_from_reduced"}));

  {
    size_t nprofiles = geovals.nprofiles(oops::Variable{"var1_sampled"}, GeoVaLFormat::SAMPLED);
    std::vector<double> expectedValues(nprofiles, 1.0);
    geovals.putAtLevel(expectedValues, oops::Variable{"var1_sampled"}, 0, GeoVaLFormat::SAMPLED);
    std::vector<double> values(nprofiles, 0.0);
    geovals.getAtLevel(values, oops::Variable{"var1_sampled"}, 0, GeoVaLFormat::SAMPLED);
    EXPECT(values == expectedValues);
  }

  {
    size_t nprofiles = geovals.nprofiles(oops::Variable{"var2_reduced"}, GeoVaLFormat::REDUCED);
    std::vector<double> expectedValues(nprofiles, 2.0);
    geovals.putAtLevel(expectedValues, oops::Variable{"var2_reduced"}, 0, GeoVaLFormat::REDUCED);
    std::vector<double> values(nprofiles, 0.0);
    geovals.getAtLevel(values, oops::Variable{"var2_reduced"}, 0, GeoVaLFormat::REDUCED);
    EXPECT(values == expectedValues);
  }

  {
    size_t nSampledProfiles =
        geovals.nprofiles(oops::Variable{"var3_sampled_aliased_with_reduced"},
          GeoVaLFormat::SAMPLED);
    size_t nReducedProfiles =
        geovals.nprofiles(oops::Variable{"var3_sampled_aliased_with_reduced"},
          GeoVaLFormat::SAMPLED);
    EXPECT_EQUAL(nSampledProfiles, nReducedProfiles);

    {
      std::vector<double> expectedValues(nSampledProfiles, 3.0);
      geovals.putAtLevel(expectedValues, oops::Variable{"var3_sampled_aliased_with_reduced"}, 0,
                         GeoVaLFormat::SAMPLED);
      {
        std::vector<double> values(nSampledProfiles, 0.0);
        geovals.getAtLevel(values, oops::Variable{"var3_sampled_aliased_with_reduced"}, 0,
          GeoVaLFormat::SAMPLED);
        EXPECT(values == expectedValues);
      }
      {
        std::vector<double> values(nReducedProfiles, 0.0);
        geovals.getAtLevel(values, oops::Variable{"var3_sampled_aliased_with_reduced"}, 0,
          GeoVaLFormat::REDUCED);
        EXPECT(values == expectedValues);
      }
    }

    {
      std::vector<double> expectedValues(nReducedProfiles, 3.5);
      geovals.putAtLevel(expectedValues, oops::Variable{"var3_sampled_aliased_with_reduced"}, 0,
                         GeoVaLFormat::REDUCED);
      {
        std::vector<double> values(nSampledProfiles, 0.0);
        geovals.getAtLevel(values, oops::Variable{"var3_sampled_aliased_with_reduced"}, 0,
          GeoVaLFormat::SAMPLED);
        EXPECT(values == expectedValues);
      }
      {
        std::vector<double> values(nReducedProfiles, 0.0);
        geovals.getAtLevel(values, oops::Variable{"var3_sampled_aliased_with_reduced"}, 0,
          GeoVaLFormat::REDUCED);
        EXPECT(values == expectedValues);
      }
    }
  }

  {
    size_t nSampledProfiles =
        geovals.nprofiles(oops::Variable{"var4_sampled_distinct_from_reduced"},
          GeoVaLFormat::SAMPLED);
    size_t nReducedProfiles =
        geovals.nprofiles(oops::Variable{"var4_sampled_distinct_from_reduced"},
          GeoVaLFormat::REDUCED);
    EXPECT(nSampledProfiles != nReducedProfiles);

    std::vector<double> expectedSampledValues(nSampledProfiles, 4.0);
    geovals.putAtLevel(expectedSampledValues, oops::Variable{"var4_sampled_distinct_from_reduced"},
                       0, GeoVaLFormat::SAMPLED);
    std::vector<double> expectedReducedValues(nReducedProfiles, 4.5);
    geovals.putAtLevel(expectedReducedValues, oops::Variable{"var4_sampled_distinct_from_reduced"},
                        0, GeoVaLFormat::REDUCED);
    {
      std::vector<double> values(nSampledProfiles, 0.0);
      geovals.getAtLevel(values, oops::Variable{"var4_sampled_distinct_from_reduced"}, 0,
        GeoVaLFormat::SAMPLED);
      EXPECT(values == expectedSampledValues);
    }
    {
      std::vector<double> values(nReducedProfiles, 0.0);
      geovals.getAtLevel(values, oops::Variable{"var4_sampled_distinct_from_reduced"}, 0,
        GeoVaLFormat::REDUCED);
      EXPECT(values == expectedReducedValues);
    }
  }
}

void testSampledAndReducedFormatsAfterOneStageConstruction() {
  oops::Variables variables({oops::Variable{"var1_sampled"},
                             oops::Variable{"var3_sampled_aliased_with_reduced"},
                             oops::Variable{"var4_sampled_distinct_from_reduced"}});
  std::vector<size_t> nlevels{1, 1, 1};
  GeoVaLs geovals(createLocationsForFormatsTest(), variables, nlevels);
  oops::Variables reducedVariables({oops::Variable{"var2_reduced"},
                                    oops::Variable{"var3_sampled_aliased_with_reduced"},
                                    oops::Variable{"var4_sampled_distinct_from_reduced"}});
  std::vector<size_t> nReducedLevels{1, 1, 1};
  geovals.addReducedVars(reducedVariables, nReducedLevels);

  testSampledAndReducedFormats(geovals);
}

void testSampledAndReducedFormatsAfterTwoStageConstruction() {
  oops::Variables variables({oops::Variable{"var1_sampled"},
                             oops::Variable{"var3_sampled_aliased_with_reduced"},
                             oops::Variable{"var4_sampled_distinct_from_reduced"}});
  GeoVaLs geovals(createLocationsForFormatsTest(), variables);
  geovals.allocate(1, variables);
  oops::Variables reducedVariables({oops::Variable{"var2_reduced"},
                                    oops::Variable{"var3_sampled_aliased_with_reduced"},
                                    oops::Variable{"var4_sampled_distinct_from_reduced"}});
  std::vector<size_t> nReducedLevels{1, 1, 1};
  geovals.addReducedVars(reducedVariables, nReducedLevels);

  testSampledAndReducedFormats(geovals);
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
    ts.emplace_back(CASE("ufo/GeoVaLs/testGeoVaLsConstructor")
      { testGeoVaLsConstructor(); });
    ts.emplace_back(CASE("ufo/GeoVaLs/testGeoVaLsWithMultipleSamplingMethods")
      { testGeoVaLsWithMultipleSamplingMethods(); });
    ts.emplace_back(CASE("ufo/GeoVaLs/testSampledAndReducedFormatsAfterOneStageConstruction")
      { testSampledAndReducedFormatsAfterOneStageConstruction(); });
    ts.emplace_back(CASE("ufo/GeoVaLs/testSampledAndReducedFormatsAfterTwoStageConstruction")
      { testSampledAndReducedFormatsAfterTwoStageConstruction(); });
  }

  void clear() const override {}
};

// ----------;

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_GEOVALS_H_
