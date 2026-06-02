/*
 * (C) 2021 Crown Copyright Met Office. All rights reserved.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef TEST_UFO_HISTORYCHECK_H_
#define TEST_UFO_HISTORYCHECK_H_

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
#include "ufo/filters/HistoryCheck.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {
void testHistoryCheck(const eckit::LocalConfiguration &conf) {
  util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  eckit::LocalConfiguration filterConf(conf, "History Check");
  ufo::HistoryCheckParameters filterParameters;
  filterParameters.validateAndDeserialize(filterConf);

  if (conf.has("air_temperatures")) {
    const std::vector<float> airTemperatures = conf.getFloatVector("air_temperatures");
    obsspace.put_db("ObsValue", "airTemperature", airTemperatures);
  }

  if (conf.has("station_ids")) {
    const std::vector<int> stationIds = conf.getIntVector("station_ids");
    obsspace.put_db("MetaData", "stationIdentification", stationIds);
  }

  if (conf.has("station_ids_string")) {
    const std::vector<std::string> stationIds = conf.getStringVector("station_ids_string");
    obsspace.put_db("MetaData", "stationIdentification", stationIds);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  ufo::HistoryCheck filter(obsspace, filterParameters, qcflags, obserr, conf);
  filter.preProcess();

  const std::vector<size_t> expectedRejectedObsIndices =
      conf.getUnsignedVector("expected rejected obs indices");
  std::vector<size_t> rejectedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::history)
      rejectedObsIndices.push_back(i);
  EXPECT_EQUAL(rejectedObsIndices, expectedRejectedObsIndices);
}

class HistoryCheck : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::HistoryCheck";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string &testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/HistoryCheck/" + testCaseName, testCaseConf)
                      {
                        testHistoryCheck(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_HISTORYCHECK_H_
