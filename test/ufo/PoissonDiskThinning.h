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

#include <boost/make_shared.hpp>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "ufo/filters/PoissonDiskThinning.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testPoissonDiskThinning(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "ObsSpace");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::comm(), bgn, end);

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

  if (conf.has("priority")) {
    const std::vector<int> priorities = conf.getIntVector("priority");
    obsspace.put_db("MetaData", "priority", priorities);
  }

  auto obserr = boost::make_shared<ioda::ObsDataVector<float>>(
      obsspace, obsspace.obsvariables(), "ObsError");
  auto qcflags = boost::make_shared<ioda::ObsDataVector<int>>(
      obsspace, obsspace.obsvariables());

  const eckit::LocalConfiguration filterConf(conf, "Poisson Disk Thinning");
  ufo::PoissonDiskThinning filter(obsspace, filterConf, qcflags, obserr);
  filter.preProcess();

  const std::vector<size_t> expectedThinnedObsIndices =
      conf.getUnsignedVector("expected_thinned_obs_indices");
  std::vector<size_t> thinnedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::thinned)
      thinnedObsIndices.push_back(i);
  EXPECT_EQUAL(thinnedObsIndices, expectedThinnedObsIndices);
}

class PoissonDiskThinning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::PoissonDiskThinning";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/PoissonDiskThinning/" + testCaseName, testCaseConf)
                      {
                        testPoissonDiskThinning(testCaseConf);
                      });
    }
  }
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_POISSONDISKTHINNING_H_
