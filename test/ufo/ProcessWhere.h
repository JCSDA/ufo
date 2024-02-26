/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PROCESSWHERE_H_
#define TEST_UFO_PROCESSWHERE_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ufo {
namespace test {

class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters);
 public:
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};
  oops::RequiredParameter<int> sizeWhereTrue{"size where true", this};
};

// -----------------------------------------------------------------------------

void testProcessWhere(const eckit::LocalConfiguration &conf,
                      bool is_in_usererror = false) {
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  eckit::LocalConfiguration obsconf(conf, "obs space");
  ioda::ObsSpace ospace(obsconf, oops::mpi::world(), timeWindow, oops::mpi::myself());
  ObsFilterData data(ospace);

  const int nlocs = conf.getInt("nlocs");
  EXPECT(data.nlocs() == nlocs);

  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ProcessWhere", confs);
  for (size_t jconf = 0; jconf < confs.size(); ++jconf) {
    eckit::LocalConfiguration config = confs[jconf];
    TestParameters params;
    params.validateAndDeserialize(config);
    if (is_in_usererror) {
      EXPECT_THROWS(processWhere(params.where, data, params.whereOperator));
    } else {
      std::vector<bool> result = processWhere(params.where, data, params.whereOperator);
      const int size_ref = params.sizeWhereTrue;
      const int size = std::count(result.begin(), result.end(), true);
      oops::Log::info() << "reference: " << size_ref << ", compare with " << size << std::endl;
      EXPECT(size == size_ref);
    }
  }
}

// -----------------------------------------------------------------------------

class ProcessWhere : public oops::Test {
 public:
  ProcessWhere() {}
  virtual ~ProcessWhere() {}

 private:
  std::string testid() const override {return "ufo::test::ProcessWhere";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ProcessWhere/testProcessWhere_successful") {
      testProcessWhere(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                "successful"));
    });

    ts.emplace_back(CASE("ufo/ProcessWhere/testProcessWhere_operator") {
      testProcessWhere(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                "operator"));
    });

    ts.emplace_back(CASE("ufo/ProcessWhere/testProcessWhere_isin_usererror") {
      testProcessWhere(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "user_error_type_handling_is_in"),
                       true);
    });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PROCESSWHERE_H_
