/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PARAMETERS_H_
#define TEST_UFO_PARAMETERS_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
namespace test {

class MyParameters : public oops::Parameters {
 public:
  oops::OptionalParameter<ufo::Variable> optVariableParameter{"opt_variable_parameter", this};
};

void testDefaultValue() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  MyParameters params;
  EXPECT(params.optVariableParameter.value() == boost::none);

  const eckit::LocalConfiguration emptyConf(conf, "empty");
  params.deserialize(emptyConf);

  EXPECT(params.optVariableParameter.value() == boost::none);
}

void testCorrectValue() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(::test::TestEnvironment::config(), "full");
  params.deserialize(fullConf);

  EXPECT(params.optVariableParameter.value() != boost::none);
  EXPECT_EQUAL(params.optVariableParameter.value().get().group(), "MetaData");
  EXPECT_EQUAL(params.optVariableParameter.value().get().variable(), "latitude");
}
void testIncorrectValue() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_opt_variable_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

class Parameters : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::Parameters";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Parameters/defaultValue") {
                      testDefaultValue();
                    });
    ts.emplace_back(CASE("ufo/Parameters/correctValue") {
                      testCorrectValue();
                    });
    ts.emplace_back(CASE("ufo/Parameters/incorrectValue") {
                      testIncorrectValue();
                    });
  }
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PARAMETERS_H_
