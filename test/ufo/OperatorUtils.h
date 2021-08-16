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
#include "test/TestEnvironment.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/OperatorUtils.h"

namespace ufo {
namespace test {

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
    util::DateTime bgn(conf.getString("window begin"));
    util::DateTime end(conf.getString("window end"));
    const eckit::LocalConfiguration obsconf(conf, "obs space");
    ioda::ObsTopLevelParameters obsparams;
    obsparams.validateAndDeserialize(obsconf);
    obsspace_.reset(new ioda::ObsSpace(obsparams, oops::mpi::world(), bgn, end,
                                       oops::mpi::myself()));
    config_ = conf.getSubConfiguration("cases");
  }

  std::shared_ptr<ioda::ObsSpace> obsspace_;
  eckit::LocalConfiguration config_;
};

void testOperatorUtils(const eckit::LocalConfiguration &conf) {
  oops::Variables expectedOperatorVariables =
      ufo::Variables(conf.getSubConfigurations("expected operator variables")).toOopsVariables();
  std::vector<int> expectedOperatorVariableIndices =
      conf.getIntVector("expected indices");

  oops::Variables operatorVariables;
  std::vector<int> operatorVariableIndices;
  ufo::getOperatorVariables(conf, TestFixture::obsspace().obsvariables(),
                            operatorVariables, operatorVariableIndices);

  EXPECT_EQUAL(operatorVariables, expectedOperatorVariables);
  EXPECT_EQUAL(operatorVariableIndices, expectedOperatorVariableIndices);
}

CASE("ufo/OperatorUtils/Without 'variables' option") {
  testOperatorUtils(TestFixture::config().getSubConfiguration("without variables"));
}

CASE("ufo/OperatorUtils/With 'variables' option") {
  testOperatorUtils(TestFixture::config().getSubConfiguration("with variables"));
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
