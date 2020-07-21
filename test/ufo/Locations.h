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
#include "oops/parallel/mpi/mpi.h"
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
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));
  const eckit::LocalConfiguration obsconf(conf, "ObsSpace");
  ioda::ObsSpace odb(obsconf, oops::mpi::comm(), bgn, end);
  const size_t nlocs = odb.nlocs();

  // testConstructor:: Locations():
  Locations locs(oops::mpi::comm());
  EXPECT(locs.nobs() == 0);
  oops::Log::test() << "Locs(eckit mpi communicator): " << locs << std::endl;

  // testConstructor:: Locations(const eckit::Configuration &)
  Locations locs1(conf, oops::mpi::comm());
  EXPECT(locs1.nobs() == nlocs);
  oops::Log::test() << "Locs(eckit constructor, eckit mpi communicator): " << locs1 << std::endl;

  // testConstructor::Locations(const ioda::ObsSpace &, const util::DateTime &,
  //                            const util::DateTime &);
  Locations locs_t(odb, bgn, end);
  EXPECT(locs_t.nobs() == nlocs);
  oops::Log::test() << "Locs(odb,t1,t2) constructor): " << locs_t << std::endl;

  // test operator+=(const Locations & other)
  util::Duration twin = end - bgn;
  util::DateTime stateTime = bgn + twin/2;

  Locations locs_b(odb, stateTime, end);
  const size_t nlocs_b = locs_b.nobs();
  oops::Log::test() << "Locs(odb,t1+(t2-t1)/2,t2) constructor): " << locs_b << std::endl;

  {
    Locations locs_a(odb, bgn, stateTime);
    const size_t nlocs_a = locs_a.nobs();
    oops::Log::test() << "Locs(odb,t1,t1+(t2-t1)/2) constructor): " << locs_a << std::endl;

    EXPECT(locs_t.nobs() == nlocs_a + nlocs_b);

    locs_a += locs_b;
    EXPECT(locs_t.nobs() == locs_a.nobs());
    oops::Log::test() << "Locs(odb,t1,t1+(t2-t1)/2) + "
                      << "Locs(odb,t1+(t2-t1)/2,t2) concatenated: " << locs_a << std::endl;
  }
  Locations locs_a(odb, bgn, stateTime);
  locs_b += locs_a;
  EXPECT(locs_t.nobs() == locs_b.nobs());
  oops::Log::test() << "Locs(odb,t1+(t2-t1)/2,t2) + "
                    << "Locs(odb,t1,t1+(t2-t1)/2) concatenated: " << locs_b << std::endl;
}

// -----------------------------------------------------------------------------

class Locations : public oops::Test {
 public:
  Locations() {}
  virtual ~Locations() {}
 private:
  std::string testid() const {return "ufo::test::Locations";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Locations/testLocations")
      { testLocations(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_LOCATIONS_H_
