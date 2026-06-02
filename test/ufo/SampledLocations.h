/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_SAMPLEDLOCATIONS_H_
#define TEST_UFO_SAMPLEDLOCATIONS_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"
#include "ufo/SampledLocations.h"

namespace ufo {
namespace test {

/// Parameters describing SampledLocations/TimeMask test
class SampledLocationsTimeTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SampledLocationsTimeTestParameters, Parameters)
 public:
  /// t1 & t2 for subsetting interpolation paths
  oops::RequiredParameter<util::DateTime> t1{"t1", this};
  oops::RequiredParameter<util::DateTime> t2{"t2", this};
  /// reference mask
  oops::RequiredParameter<std::vector<bool>> refTimeMask{"reference mask", this};
};

/// Parameters describing SampledLocations test
class SampledLocationsTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SampledLocationsTestParameters, Parameters)
 public:
  /// reference longitudes
  oops::RequiredParameter<std::vector<float>> refLons{"reference lons", this};
  /// reference latitudes
  oops::RequiredParameter<std::vector<float>> refLats{"reference lats", this};
  /// set of tests for time masks
  oops::RequiredParameter<std::vector<SampledLocationsTimeTestParameters>>
        timeMaskTests{"time mask tests", this};
};

/// Code shared by all SampledLocations tests
class SampledLocationsTestFixture : private boost::noncopyable {
 public:
  static const ioda::ObsSpace & obsspace() {return *getInstance().obsspace_;}
  static const eckit::LocalConfiguration & testconfig() {return getInstance().testconfig_;}
  static const SampledLocationsTestParameters & testparams() {
    return getInstance().testparams_;
  }

 private:
  static SampledLocationsTestFixture & getInstance() {
    static SampledLocationsTestFixture theSampledLocationsTestFixture;
    return theSampledLocationsTestFixture;
  }

  SampledLocationsTestFixture() {
    const eckit::Configuration & conf = ::test::TestEnvironment::config();
    const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));
    const eckit::LocalConfiguration obsconf(conf, "obs space");
    obsspace_.reset(new ioda::ObsSpace(obsconf, oops::mpi::world(), timeWindow,
                                       oops::mpi::myself()));
    testconfig_ = conf.getSubConfiguration("locations test");
    testparams_.validateAndDeserialize(testconfig_);
  }

  std::shared_ptr<ioda::ObsSpace> obsspace_;
  eckit::LocalConfiguration       testconfig_;
  SampledLocationsTestParameters         testparams_;
};

// -----------------------------------------------------------------------------
/// Tests SampledLocations constructors and method size()
void testSampledLocations() {
  typedef SampledLocationsTestFixture Test_;

  const size_t nlocs = Test_::obsspace().nlocs();

  // testConstructor:: SampledLocations(const eckit::Configuration &)
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  SampledLocations locs1(conf, oops::mpi::world());
  EXPECT(locs1.size() == nlocs);
  oops::Log::test() << "SampledLocations(configuration, eckit mpi communicator): "
                    << locs1 << std::endl;

  // testConstructor::SampledLocations(const ioda::ObsSpace &);
/*  SampledLocations locs2(Test_::obsspace());
  EXPECT(locs2.size() == nlocs);
  oops::Log::test() << "SampledLocations(obsspace) constructor: "
                    << locs2 << std::endl;*/
}

// -----------------------------------------------------------------------------
/// Tests SampledLocations accessors lons() and lats()
void testLonsLats() {
  typedef SampledLocationsTestFixture Test_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  SampledLocations locs(conf, oops::mpi::world());

  const SampledLocationsTestParameters & params = Test_::testparams();
  const float abstol = 1.0e-8;
  EXPECT(oops::are_all_close_absolute(params.refLons.value(), locs.lons(), abstol));
  EXPECT(oops::are_all_close_absolute(params.refLats.value(), locs.lats(), abstol));
}

// -----------------------------------------------------------------------------
extern "C" {
  /// similar to LonsLats test, on Fortran level
  /// Returns 1 if the test passes, 0 if the test fails
  int test_sampled_locations_lonslats_f90(const SampledLocations &,
                                                      const eckit::Configuration &);
}

/// Tests SampledLocations accessors lons() and lats() from Fortran
void testFortranLonsLats() {
  typedef SampledLocationsTestFixture Test_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  SampledLocations locs(conf, oops::mpi::world());

  EXPECT(test_sampled_locations_lonslats_f90(locs, Test_::testconfig()));
}

// -----------------------------------------------------------------------------
/// Tests SampledLocations::isInTimeWindow method
void testTimeMask() {
  typedef SampledLocationsTestFixture Test_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  SampledLocations locs(conf, oops::mpi::world());

  const SampledLocationsTestParameters & params = Test_::testparams();
  for (const auto & test : params.timeMaskTests.value()) {
    EXPECT_EQUAL(test.refTimeMask.value(), locs.isInTimeWindow(test.t1, test.t2));
  }
}

// -----------------------------------------------------------------------------
extern "C" {
  /// similar to TimeMask test, on Fortran level
  /// Returns 1 if the test passes, 0 if the test fails
  int test_sampled_locations_timemask_f90(const SampledLocations &,
                                                      const eckit::Configuration &);
}

/// Tests SampledLocations::isInTimeWindow method from Fortran
void testFortranTimeMask() {
  typedef SampledLocationsTestFixture Test_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  SampledLocations locs(conf, oops::mpi::world());

  EXPECT(test_sampled_locations_timemask_f90(locs, Test_::testconfig()));
}

// -----------------------------------------------------------------------------
/// Tests operator+= (concatenation of two SampledLocations)
void testConcatenate() {
  typedef SampledLocationsTestFixture Test_;

  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  SampledLocations locs1(conf, oops::mpi::world());
  SampledLocations locs2(conf, oops::mpi::world());
  const size_t nlocs = locs1.size();

  locs1 += locs2;

  EXPECT_EQUAL(locs1.size(), 2*nlocs);

  const float abstol = 1.0e-8;

  /// test that concatenated longitudes have two copies of original lons
  const std::vector<float> & lons = locs1.lons();
  std::vector<float> lons1(lons.begin(), lons.begin() + nlocs);
  std::vector<float> lons2(lons.begin() + nlocs, lons.end());
  EXPECT(oops::are_all_close_absolute(lons1, lons2, abstol));
  EXPECT(oops::are_all_close_absolute(lons1, locs2.lons(), abstol));

  /// test that concatenated latitudes have two copies of original lats
  const std::vector<float> & lats = locs1.lats();
  std::vector<float> lats1(lats.begin(), lats.begin() + nlocs);
  std::vector<float> lats2(lats.begin() + nlocs, lats.end());
  EXPECT(oops::are_all_close_absolute(lats1, lats2, abstol));
  EXPECT(oops::are_all_close_absolute(lats1, locs2.lats(), abstol));
}

// -----------------------------------------------------------------------------

class SampledLocations : public oops::Test {
 public:
  SampledLocations() {}
  virtual ~SampledLocations() = default;

 private:
  std::string testid() const override {return "ufo::test::SampledLocations";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/SampledLocations/testSampledLocations")
      { testSampledLocations(); });
    ts.emplace_back(CASE("ufo/SampledLocations/testLonsLats")
      { testLonsLats(); });
    ts.emplace_back(CASE("ufo/SampledLocations/testFortranLonsLats")
      { testFortranLonsLats(); });
    ts.emplace_back(CASE("ufo/SampledLocations/testTimeMask")
      { testTimeMask(); });
    ts.emplace_back(CASE("ufo/SampledLocations/testFortranTimeMask")
      { testFortranTimeMask(); });
    ts.emplace_back(CASE("ufo/SampledLocations/testConcatenation")
      { testConcatenate(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_SAMPLEDLOCATIONS_H_
