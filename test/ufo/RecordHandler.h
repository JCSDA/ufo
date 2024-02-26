/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_RECORDHANDLER_H_
#define TEST_UFO_RECORDHANDLER_H_

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "ioda/ObsSpace.h"

#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/missingValues.h"
#include "oops/util/TimeWindow.h"

#include "ufo/utils/RecordHandler.h"

namespace ufo {
namespace test {

void testRecordHandler(const eckit::LocalConfiguration &conf) {
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  // Obtain air_temperature and eastward_wind from configuration and save to ObsSpace.
  std::vector<float> air_temperature = conf.getFloatVector("air_temperature");
  std::vector<float> eastward_wind = conf.getFloatVector("eastward_wind");

  // Create a vector of QC flags based on the values of the observed values.
  ioda::ObsDataVector<int> qcflags(obsspace, obsspace.obsvariables());

  // Change -999 to the missing floating-point value.
  // Also set the relevant QC flag to missing.
  const float missingFloat = util::missingValue<float>();
  for (size_t jloc = 0; jloc < obsspace.nlocs(); ++jloc) {
    if (air_temperature[jloc] == -999) {
      air_temperature[jloc] = missingFloat;
      qcflags[0][jloc] = 10;
    }
    if (eastward_wind[jloc] == -999) {
      eastward_wind[jloc] = missingFloat;
      qcflags[1][jloc] = 10;
    }
  }
  obsspace.put_db("ObsValue", "air_temperature", air_temperature);
  obsspace.put_db("ObsValue", "eastward_wind", eastward_wind);

  const Variables filtervars = Variables(obsspace.obsvariables());

  // Create a record handler.
  ufo::RecordHandler recordHandler(obsspace,
                                   filtervars,
                                   qcflags,
                                   conf.getBool("retainOnlyIfAllFilterVariablesAreValid", false));

  // (1) Get input vector.
  // (There is not a getBoolVector option for eckit::Configuration, which is why the conversion from
  // int to bool is performed.)
  const std::vector<int> inputInt = conf.getIntVector("input");
  const std::vector<bool> input(inputInt.begin(), inputInt.end());

  // (2) Check treatment of apply vector

  // Get modified apply vector.
  const std::vector<bool> modified_apply = recordHandler.changeApplyIfRecordsAreSingleObs(input);

  // Get expected modified apply vector.
  const std::vector<int> expected_applyInt = conf.getIntVector("expected_apply");
  const std::vector<bool> expected_apply(expected_applyInt.begin(), expected_applyInt.end());

  EXPECT(expected_apply.size() == modified_apply.size());
  for (size_t idx = 0; idx < modified_apply.size(); ++idx)
    EXPECT(modified_apply[idx] == expected_apply[idx]);

  // (3) Check treatment of isThinned vector.

  // Get modified isThinned vector.
  const std::vector<bool> modified_isThinned =
    recordHandler.changeThinnedIfRecordsAreSingleObs(input);

  // Get expected modified isThinned vector.
  const std::vector<int> expected_isThinnedInt = conf.getIntVector("expected_isThinned");
  const std::vector<bool> expected_isThinned(expected_isThinnedInt.begin(),
                                             expected_isThinnedInt.end());

  EXPECT(expected_isThinned.size() == modified_isThinned.size());
  for (size_t idx = 0; idx < modified_isThinned.size(); ++idx)
    EXPECT(modified_isThinned[idx] == expected_isThinned[idx]);
}

class RecordHandler : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::RecordHandler";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/RecordHandler/" + testCaseName, testCaseConf)
                      {
                        testRecordHandler(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_RECORDHANDLER_H_

