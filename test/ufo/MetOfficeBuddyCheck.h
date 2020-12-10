/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_METOFFICEBUDDYCHECK_H_
#define TEST_UFO_METOFFICEBUDDYCHECK_H_

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
#include "ufo/filters/MetOfficeBuddyCheck.h"
#include "ufo/filters/Variables.h"
#include "ufo/Locations.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {
namespace test {

void testMetOfficeBuddyCheck(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsSpace(obsSpaceConf, oops::mpi::world(), bgn, end, oops::mpi::myself());

  const eckit::LocalConfiguration floatVarInitConf(conf, "FloatVariables");
  for (const std::string & varNameGroup : floatVarInitConf.keys()) {
    std::string varName, varGroup;
    ufo::splitVarGroup(varNameGroup, varName, varGroup);
    const std::vector<float> values = floatVarInitConf.getFloatVector(varNameGroup);
    obsSpace.put_db(varGroup, varName, values);
  }

  const eckit::LocalConfiguration intVarInitConf(conf, "IntVariables");
  for (const std::string & varNameGroup : intVarInitConf.keys()) {
    std::string varName, varGroup;
    ufo::splitVarGroup(varNameGroup, varName, varGroup);
    const std::vector<int> values = intVarInitConf.getIntVector(varNameGroup);
    obsSpace.put_db(varGroup, varName, values);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsSpace, obsSpace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsSpace, obsSpace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "Met Office Buddy Check");
  ufo::MetOfficeBuddyCheck filter(obsSpace, filterConf, qcflags, obserr);
  filter.preProcess();

  ioda::ObsVector hofx(obsSpace, "HofX");
  ufo::Locations locations(obsSpace);
  ufo::ObsDiagnostics obsDiags(obsSpace, locations, oops::Variables());
  filter.postFilter(hofx, obsDiags);

  const eckit::LocalConfiguration pgeConf(conf, "ExpectedGrossErrorProbabilities");
  for (const std::string & varNameGroup : pgeConf.keys()) {
    std::string varName, varGroup;
    ufo::splitVarGroup(varNameGroup, varName, varGroup);
    const std::vector<float> expectedPges = pgeConf.getFloatVector(varNameGroup);
    std::vector<float> actualPges(obsSpace.nlocs());
    obsSpace.get_db(varGroup, varName, actualPges);
    EXPECT(oops::are_all_close_absolute(actualPges, expectedPges, 1e-6f));
  }
}

class MetOfficeBuddyCheck : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::MetOfficeBuddyCheck";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/MetOfficeBuddyCheck/" + testCaseName, testCaseConf)
                      {
                        testMetOfficeBuddyCheck(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_METOFFICEBUDDYCHECK_H_
