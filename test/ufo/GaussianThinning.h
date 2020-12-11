/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_GAUSSIANTHINNING_H_
#define TEST_UFO_GAUSSIANTHINNING_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/Gaussian_Thinning.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testGaussianThinning(const eckit::LocalConfiguration &conf) {
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

  eckit::LocalConfiguration filterConf(conf, "GaussianThinning");
  ufo::GaussianThinningParameters filterParameters;
  filterParameters.validateAndDeserialize(filterConf);
  ufo::Gaussian_Thinning filter(obsspace, filterParameters, qcflags, obserr);
  filter.preProcess();

  const std::vector<size_t> expectedThinnedObsIndices =
      conf.getUnsignedVector("expected_thinned_obs_indices");
  std::vector<size_t> thinnedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::thinned)
      thinnedObsIndices.push_back(i);
  EXPECT_EQUAL(thinnedObsIndices.size(), expectedThinnedObsIndices.size());
  const bool equal = std::equal(thinnedObsIndices.begin(), thinnedObsIndices.end(),
                                expectedThinnedObsIndices.begin());
  EXPECT(equal);
}

class GaussianThinning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::GaussianThinning";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/GaussianThinning/" + testCaseName, testCaseConf)
                      {
                        testGaussianThinning(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_GAUSSIANTHINNING_H_
