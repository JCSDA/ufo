/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_STUCKCHECK_H_
#define TEST_UFO_STUCKCHECK_H_

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
#include "ufo/filters/StuckCheck.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testStuckCheck(const eckit::LocalConfiguration &conf) {
  util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  if (conf.has("air_temperatures")) {
    const std::vector<float> airTemperatures = conf.getFloatVector("air_temperatures");
    obsspace.put_db("ObsValue", "airTemperature", airTemperatures);
  }

  if (conf.has("air_pressures")) {
    const std::vector<float> airPressures = conf.getFloatVector("air_pressures");
    obsspace.put_db("ObsValue", "pressure", airPressures);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  eckit::LocalConfiguration filterConf(conf, "Stuck Check");
  ufo::StuckCheckParameters filterParameters;
  filterParameters.validateAndDeserialize(filterConf);
  ufo::StuckCheck filter(obsspace, filterParameters, qcflags, obserr);
  filter.preProcess();

  const std::vector<size_t> expectedRejectedObsIndices =
      conf.getUnsignedVector("expected_rejected_obs_indices");
  std::vector<size_t> rejectedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::track)
      rejectedObsIndices.push_back(i);
  EXPECT_EQUAL(rejectedObsIndices, expectedRejectedObsIndices);
}

class StuckCheck : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::StuckCheck";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/StuckCheck/" + testCaseName, testCaseConf)
                      {
                        testStuckCheck(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_STUCKCHECK_H_
