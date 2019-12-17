/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_TEMPORALTHINNING_H_
#define TEST_UFO_TEMPORALTHINNING_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_shared.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "ufo/filters/TemporalThinning.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testTemporalThinning(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "ObsSpace");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::comm(), bgn, end);

  if (conf.has("category")) {
    const std::vector<int> categories = conf.getIntVector("category");
    obsspace.put_db("MetaData", "category", categories);
  }

  if (conf.has("priority")) {
    const std::vector<int> priorities = conf.getIntVector("priority");
    obsspace.put_db("MetaData", "priority", priorities);
  }

  auto obserr = boost::make_shared<ioda::ObsDataVector<float>>(
      obsspace, obsspace.obsvariables(), "ObsError");
  auto qcflags = boost::make_shared<ioda::ObsDataVector<int>>(
      obsspace, obsspace.obsvariables());

  const eckit::LocalConfiguration filterConf(conf, "TemporalThinning");
  ufo::TemporalThinning filter(obsspace, filterConf, qcflags, obserr);
  filter.preProcess();

  const std::vector<size_t> expectedThinnedObsIndices =
      conf.getUnsignedVector("expected_thinned_obs_indices");
  std::vector<size_t> thinnedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::thinned)
      thinnedObsIndices.push_back(i);
  EXPECT_EQUAL(thinnedObsIndices, expectedThinnedObsIndices);
}

CASE("ufo/TemporalThinning/Min_spacing below observation spacing") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Min_spacing below observation spacing"));
}

CASE("ufo/TemporalThinning/Min spacing equal to observation spacing") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Min spacing equal to observation spacing"));
}

CASE("ufo/TemporalThinning/Min spacing above observation spacing") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Min spacing above observation spacing"));
}

CASE("ufo/TemporalThinning/Categories") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Categories"));
}

CASE("ufo/TemporalThinning/Categories, observations sorted in descending order") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Categories, observations sorted "
                                                 "in descending order"));
}

CASE("ufo/TemporalThinning/Categories, where clause") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Categories, where clause"));
}

CASE("ufo/TemporalThinning/Tolerance and priorities, "
     "first observation in each group to be retained") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Tolerance and priorities, first observation in "
                                                 "each group to be retained"));
}

CASE("ufo/TemporalThinning/Tolerance and priorities, "
     "second observation in each group to be retained") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Tolerance and priorities, second observation in "
                                                 "each group to be retained"));
}

CASE("ufo/TemporalThinning/Tolerance and priorities, "
     "third observation in each group to be retained") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Tolerance and priorities, third observation in "
                                                 "each group to be retained"));
}

CASE("ufo/TemporalThinning/Tolerance but no priorities") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Tolerance but no priorities"));
}

CASE("ufo/TemporalThinning/Seed time inside observation time range "
     "(should pick preceding observation)") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time inside observation time range (should "
                                                 "pick preceding observation)"));
}

CASE("ufo/TemporalThinning/Seed time inside observation time range "
     "(should pick following observation)") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time inside observation time range "
                                                 "(should pick following observation)"));
}

CASE("ufo/TemporalThinning/Seed time midway between two observations "
     "(should pick following observation)") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time midway between two observations "
                                                 "(should pick following observation)"));
}

CASE("ufo/TemporalThinning/Seed time at earliest observation") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time at earliest observation"));
}

CASE("ufo/TemporalThinning/Seed time before earliest observation") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time before earliest observation"));
}

CASE("ufo/TemporalThinning/Seed time at latest observation") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time at latest observation"));
}

CASE("ufo/TemporalThinning/Seed time after latest observation") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Seed time after latest observation"));
}

CASE("ufo/TemporalThinning/Tolerance, priorities and seed time at a low-priority observation "
     "followed by a high-priority one") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Tolerance, priorities and seed time at "
                                                 "a low-priority observation followed by "
                                                 "a high-priority one"));
}

CASE("ufo/TemporalThinning/Tolerance, priorities and seed time at a low-priority observation "
     "preceded by a high-priority one") {
  testTemporalThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Tolerance, priorities and seed time at "
                                                 "a low-priority observation preceded by "
                                                 "a high-priority one"));
}

class TemporalThinning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::TemporalThinning";}

  void register_tests() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_TEMPORALTHINNING_H_
