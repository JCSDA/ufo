/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PROFILECONSISTENCYCHECKS_H_
#define TEST_UFO_PROFILECONSISTENCYCHECKS_H_

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
#include "test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "ufo/filters/ProfileConsistencyChecks.h"
#include "ufo/filters/Variables.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {
namespace test {

void testProfileConsistencyChecks(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "ObsSpace");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::comm(), bgn, end);

  ioda::ObsVector hofx(obsspace);

  const eckit::LocalConfiguration obsdiagconf(conf, "ObsDiag");
  std::vector<eckit::LocalConfiguration> varconfs;
  obsdiagconf.get("variables", varconfs);
  const Variables diagvars(varconfs);
  const ObsDiagnostics obsdiags(obsdiagconf, obsspace, diagvars.toOopsVariables());

  auto obserr = boost::make_shared<ioda::ObsDataVector<float>>(
      obsspace, obsspace.obsvariables(), "ObsError");

  auto qcflags = boost::make_shared<ioda::ObsDataVector<int>>(
      obsspace, obsspace.obsvariables());

  const eckit::LocalConfiguration filterConf(conf, "ProfileConsistencyChecks");
  ufo::ProfileConsistencyChecks filter(obsspace, filterConf, qcflags, obserr);

  filter.preProcess();
  filter.postFilter(hofx, obsdiags);

  // Check there are no mismatches between the values produced by this code and the OPS equivalents
  for (auto nMM : filter.getMismatches()) {
      EXPECT_EQUAL(nMM, 0);
    }
}

class ProfileConsistencyChecks : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::ProfileConsistencyChecks";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/ProfileConsistencyChecks/" + testCaseName, testCaseConf)
                      {
                        testProfileConsistencyChecks(testCaseConf);
                      });
    }
  }
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PROFILECONSISTENCYCHECKS_H_
