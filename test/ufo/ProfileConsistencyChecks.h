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

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ProfileConsistencyCheckParameters.h"
#include "ufo/filters/ProfileConsistencyChecks.h"
#include "ufo/filters/Variables.h"
#include "ufo/ObsDiagnostics.h"

#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ProfileCheckBase.h"
#include "ufo/profile/ProfileCheckUInterp.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {
namespace test {

void testProfileConsistencyChecks(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::comm(), bgn, end);

  ioda::ObsVector hofx(obsspace);

  const eckit::LocalConfiguration obsdiagconf(conf, "obs diagnostics");
  std::vector<eckit::LocalConfiguration> varconfs;
  obsdiagconf.get("variables", varconfs);
  const Variables diagvars(varconfs);
  const ObsDiagnostics obsdiags(obsdiagconf, obsspace, diagvars.toOopsVariables());

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));

  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "ProfileConsistencyChecks");

  // Determine whether an exception is expected to be thrown.
  // Exceptions can be thrown in two places: on instantiation of the filter,
  // and during the operation of the filter.
  bool expectThrowOnInstantiation = conf.getBool("ExpectThrowOnInstantiation", false);
  bool expectThrowDuringOperation = conf.getBool("ExpectThrowDuringOperation", false);

  if (expectThrowOnInstantiation) {
    EXPECT_THROWS(ufo::ProfileConsistencyChecks filterThrow(obsspace, filterConf, qcflags, obserr));
    // Do not proceed further in this case.
    return;
  }

  ufo::ProfileConsistencyChecks filter(obsspace, filterConf, qcflags, obserr);
  filter.preProcess();

  if (expectThrowDuringOperation)
    EXPECT_THROWS(filter.postFilter(hofx, obsdiags));
  else
    filter.postFilter(hofx, obsdiags);

  // Determine whether the mismatch check should be bypassed or not.
  // It might be necessary to disable the mismatch check in tests which are
  // designed to reach code paths that would normally result in failure.
  bool bypassMismatchComparison = conf.getBool("BypassMismatchComparison", false);

  // Check there are no mismatches between the values produced by this code and the OPS equivalents
  if (!bypassMismatchComparison) {
    for (auto nMM : filter.getMismatches())
      EXPECT_EQUAL(nMM, 0);
  }

  // === Additional tests of exceptions === //

  // Test whether adding the same check twice throws an exception.
  bool addDuplicateCheck = conf.getBool("AddDuplicateCheck", false);
  if (addDuplicateCheck) {
    static ProfileCheckMaker<ProfileCheckUInterp> makerDuplicate1_("duplicate");
    EXPECT_THROWS(static ProfileCheckMaker<ProfileCheckUInterp> makerDuplicate2_("duplicate"));
  }

  // Test whether using the get function with the wrong type throws an exception.
  bool getWrongType = conf.getBool("GetWrongType", false);
  if (getWrongType) {
    std::unique_ptr <ProfileConsistencyCheckParameters> options_;
    options_.reset(new ProfileConsistencyCheckParameters());
    options_->deserialize(conf);
    EntireSampleDataHandler entireSampleDataHandler(obsspace,
                                                    *options_);
    // Load data from obsspace
    entireSampleDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    // Attempt to access data with incorrect type
    EXPECT_THROWS(entireSampleDataHandler.get<int>(ufo::VariableNames::obs_air_pressure));
    std::vector<bool> apply(obsspace.nlocs(), true);
    ProfileIndices profileIndices(obsspace,
                                  *options_,
                                  apply);
    ProfileDataHandler profileDataHandler(obsspace,
                                          *options_,
                                          entireSampleDataHandler,
                                          profileIndices);
    // Obtain profile data
    profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    // Attempt to access data with incorrect type
    EXPECT_THROWS(profileDataHandler.get<int>(ufo::VariableNames::obs_air_pressure));
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
