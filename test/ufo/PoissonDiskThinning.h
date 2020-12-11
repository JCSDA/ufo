/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_POISSONDISKTHINNING_H_
#define TEST_UFO_POISSONDISKTHINNING_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/PoissonDiskThinning.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testPoissonDiskThinning(const eckit::LocalConfiguration &conf,
                             bool expectValidationError = false) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), bgn, end, oops::mpi::myself());

  if (conf.has("air_pressures")) {
    const std::vector<float> air_pressures = conf.getFloatVector("air_pressures");
    obsspace.put_db("MetaData", "air_pressure", air_pressures);
    const std::vector<float> air_pressure_obserrors(air_pressures.size(), 1.0f);
    obsspace.put_db("ObsError", "air_pressure", air_pressure_obserrors);
  }

  if (conf.has("min_vertical_spacing")) {
    const std::vector<float> min_vertical_spacing = conf.getFloatVector("min_vertical_spacing");
    obsspace.put_db("MetaData", "min_vertical_spacing", min_vertical_spacing);
  }

  if (conf.has("category")) {
    const std::vector<int> categories = conf.getIntVector("category");
    obsspace.put_db("MetaData", "category", categories);
  }

  if (conf.has("string_category")) {
    const std::vector<std::string> categories = conf.getStringVector("string_category");
    obsspace.put_db("MetaData", "string_category", categories);
  }

  if (conf.has("priority")) {
    const std::vector<int> priorities = conf.getIntVector("priority");
    obsspace.put_db("MetaData", "priority", priorities);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "Poisson Disk Thinning");
  ufo::PoissonDiskThinning filter(obsspace, filterConf, qcflags, obserr);
  if (expectValidationError) {
    EXPECT_THROWS(filter.preProcess());
    return;
  } else {
    filter.preProcess();
  }

  const std::vector<size_t> expectedThinnedObsIndices =
      conf.getUnsignedVector("expected_thinned_obs_indices");
  std::vector<size_t> thinnedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::thinned)
      thinnedObsIndices.push_back(i);
  EXPECT_EQUAL(thinnedObsIndices, expectedThinnedObsIndices);
}

CASE("ufo/PoissonDiskThinning/No thinning") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "No thinning"));
}

CASE("ufo/PoissonDiskThinning/"
     "Horizontal thinning, min spacing smaller than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                    "Horizontal thinning, min spacing "
                                                    "smaller than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Horizontal thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal thinning, min spacing "
                                                 "larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Vertical thinning, min spacing smaller than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical thinning, min spacing "
                                                 "smaller than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Vertical thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical thinning, min spacing "
                                                 "larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/Vertical thinning, where clause") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical thinning, where clause"));
}

CASE("ufo/PoissonDiskThinning/"
     "Time thinning, min spacing equal to nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Time thinning, min spacing "
                                                 "equal to nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Time thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Time thinning, min spacing "
                                                 "larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Horizontal and vertical thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal and vertical thinning, min spacing "
                                                 "larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Horizontal and time thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal and time thinning, min spacing "
                                                 "larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Vertical and time thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical and time thinning, min spacing "
                                                 "larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/"
     "Horizontal, vertical and time thinning, min spacing larger than nearest neighbor spacing") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal, vertical and time thinning, min "
                                                 "spacing larger than nearest neighbor spacing"));
}

CASE("ufo/PoissonDiskThinning/Priorities") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Priorities"));
}

CASE("ufo/PoissonDiskThinning/Int-valued categories") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Int-valued categories"));
}

CASE("ufo/PoissonDiskThinning/String-valued categories") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "String-valued categories"));
}

CASE("ufo/PoissonDiskThinning/Variable min spacings") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Variable min spacings"));
}

CASE("ufo/PoissonDiskThinning/Variable min spacings, shuffling") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Variable min spacings, shuffling"));
}

CASE("ufo/PoissonDiskThinning/Cylindrical exclusion volumes") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Cylindrical exclusion volumes"));
}

CASE("ufo/PoissonDiskThinning/Ellipsoidal exclusion volumes") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Ellipsoidal exclusion volumes"));
}

CASE("ufo/PoissonDiskThinning/Incorrectly ordered min horizontal spacings") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Incorrectly ordered min horizontal spacings"),
                          true /*expectValidationError?*/);
}

CASE("ufo/PoissonDiskThinning/Incorrectly ordered min vertical spacings") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Incorrectly ordered min vertical spacings"),
                          true /*expectValidationError?*/);
}

CASE("ufo/PoissonDiskThinning/Incorrectly ordered min time spacings") {
  testPoissonDiskThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Incorrectly ordered min time spacings"),
                          true /*expectValidationError?*/);
}

class PoissonDiskThinning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::PoissonDiskThinning";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_POISSONDISKTHINNING_H_
