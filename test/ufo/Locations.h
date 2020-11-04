/*
 * (C) Copyright 2019 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_LOCATIONS_H_
#define TEST_UFO_LOCATIONS_H_

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
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/Locations.h"

namespace ufo {
namespace test {
// -----------------------------------------------------------------------------

void testLocations() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  //  Setup ObsSpace
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  ioda::ObsSpace odb(obsconf, oops::mpi::world(), bgn, end, oops::mpi::myself());
  const size_t nlocs = odb.nlocs();

  // testConstructor:: Locations():
  Locations locs(oops::mpi::world());
  EXPECT(locs.nobs() == 0);
  oops::Log::test() << "Locs(eckit mpi communicator): " << locs << std::endl;

  // testConstructor:: Locations(const eckit::Configuration &)
  Locations locs1(conf, oops::mpi::world());
  EXPECT(locs1.nobs() == nlocs);
  oops::Log::test() << "Locs(eckit constructor, eckit mpi communicator): " << locs1 << std::endl;

  // testConstructor::Locations(const ioda::ObsSpace &, const util::DateTime &,
  //                            const util::DateTime &);
  Locations locs_t(odb);
  EXPECT(locs_t.nobs() == nlocs);
  oops::Log::test() << "Locs(odb,t1,t2) constructor): " << locs_t << std::endl;
}

// -----------------------------------------------------------------------------

class Locations : public oops::Test {
 public:
  Locations() {}
  virtual ~Locations() {}
 private:
  std::string testid() const override {return "ufo::test::Locations";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Locations/testLocations")
      { testLocations(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_LOCATIONS_H_
