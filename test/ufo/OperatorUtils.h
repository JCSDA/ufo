/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OPERATORUTILS_H_
#define TEST_UFO_OPERATORUTILS_H_

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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/OperatorUtils.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

/// \brief Options used to configure the testing of code in OperatorUtils.h
class testOperatorUtilsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(testOperatorUtilsParameters, Parameters)

 public:
  oops::RequiredParameter<std::vector<Variable>> expectedOperatorVariables{
    "expected operator variables",
    "List of expected operator variables from test.",
    this};

  oops::RequiredParameter<std::vector<int>> expectedOperatorVariableIndices{
    "expected indices",
    "List of expected indices from test.",
    this};

  oops::OptionalParameter<std::vector<Variable>> optionalVariables{
    "variables",
    "An optional list of variables.  If this is not used the variables in the ObsSpace "
    "are used.",
    this};
};

// -----------------------------------------------------------------------------

/// Code shared by all tests
class TestFixture : private boost::noncopyable {
 public:
  static const ioda::ObsSpace & obsspace() {return *getInstance().obsspace_;}
  static const eckit::LocalConfiguration & config() {return getInstance().config_;}

 private:
  static TestFixture & getInstance() {
    static TestFixture theTestFixture;
    return theTestFixture;
  }

  TestFixture() {
    const eckit::Configuration & conf = ::test::TestEnvironment::config();
    const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));
    const eckit::LocalConfiguration obsconf(conf, "obs space");
    obsspace_.reset(new ioda::ObsSpace(obsconf, oops::mpi::world(), timeWindow,
                                       oops::mpi::myself()));
    config_ = conf.getSubConfiguration("cases");
  }

  std::shared_ptr<ioda::ObsSpace> obsspace_;
  eckit::LocalConfiguration config_;
};

void testOperatorUtils(const testOperatorUtilsParameters &parameters) {
  oops::ObsVariables expectedOperatorVariables =
    ufo::Variables(parameters.expectedOperatorVariables.value()).toOopsObsVariables();

  std::vector<int> expectedOperatorVariableIndices =
    parameters.expectedOperatorVariableIndices.value();

  oops::ObsVariables operatorVariables;
  std::vector<int> operatorVariableIndices;
  ufo::getOperatorVariables(parameters.optionalVariables.value(),
                            TestFixture::obsspace().obsvariables(),
                            operatorVariables, operatorVariableIndices);

  EXPECT_EQUAL(operatorVariables, expectedOperatorVariables);
  EXPECT_EQUAL(operatorVariableIndices, expectedOperatorVariableIndices);
}

CASE("ufo/OperatorUtils/Without 'variables' option") {
  testOperatorUtilsParameters parameters;
  parameters.validateAndDeserialize(
    TestFixture::config().getSubConfiguration("without variables"));
  testOperatorUtils(parameters);
}

CASE("ufo/OperatorUtils/With 'variables' option") {
  testOperatorUtilsParameters parameters;
  parameters.validateAndDeserialize(
    TestFixture::config().getSubConfiguration("with variables"));
  testOperatorUtils(parameters);
}

class OperatorUtils : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::OperatorUtils";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OPERATORUTILS_H_
