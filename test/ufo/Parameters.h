/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PARAMETERS_H_
#define TEST_UFO_PARAMETERS_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_shared.hpp>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "../ufo/Expect.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "ufo/utils/parameters/EnumParameter.h"
#include "ufo/utils/parameters/OptionalParameter.h"
#include "ufo/utils/parameters/Parameter.h"
#include "ufo/utils/parameters/Parameters.h"

namespace ufo {
namespace test {

enum class Fruit {
  APPLE, ORANGE
};

}  // namespace test
}  // namespace ufo

namespace ufo {

template <>
inline test::Fruit enumFromString(const std::string &s) {
  if (s == "apple")
    return test::Fruit::APPLE;
  if (s == "orange")
    return test::Fruit::ORANGE;
  throw eckit::BadParameter("Bad conversion from std::string '" + s + "' to Fruit", Here());
}

}  // namespace ufo

namespace ufo {
namespace test {

class MyParameters : public ufo::Parameters {
 public:
  Parameter<float> floatParameter{"float_parameter", 1.5f, this};
  Parameter<int> intParameter{"int_parameter", 2, this};
  Parameter<bool> boolParameter{"bool_parameter", true, this};
  OptionalParameter<float> optFloatParameter{"opt_float_parameter", this};
  OptionalParameter<ufo::Variable> optVariableParameter{"opt_variable_parameter", this};
  OptionalParameter<util::DateTime> optDateTimeParameter{"opt_date_time_parameter", this};
  OptionalParameter<util::Duration> optDurationParameter{"opt_duration_parameter", this};
  EnumParameter<Fruit> fruitParameter{"fruit_parameter", Fruit::ORANGE, this};
};

void testDefaultValues() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  MyParameters params;

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optVariableParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT(params.fruitParameter == Fruit::ORANGE);

  const eckit::LocalConfiguration emptyConf(conf, "empty");
  params.deserialize(emptyConf);

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optVariableParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT(params.fruitParameter == Fruit::ORANGE);
}

void testCorrectValues() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(::test::TestEnvironment::config(), "full");
  params.deserialize(fullConf);

  EXPECT_EQUAL(params.floatParameter, 3.5f);
  EXPECT_EQUAL(params.intParameter, 4);
  EXPECT(!params.boolParameter);
  EXPECT(params.optFloatParameter.value() != boost::none);
  EXPECT_EQUAL(params.optFloatParameter.value().get(), 5.5f);
  EXPECT(params.optVariableParameter.value() != boost::none);
  EXPECT_EQUAL(params.optVariableParameter.value().get().group(), "MetaData");
  EXPECT_EQUAL(params.optVariableParameter.value().get().variable(), "latitude");
  EXPECT(params.optDateTimeParameter.value() != boost::none);
  EXPECT_EQUAL(params.optDateTimeParameter.value().get(), util::DateTime(2010, 2, 3, 4, 5, 6));
  EXPECT(params.optDurationParameter.value() != boost::none);
  EXPECT_EQUAL(params.optDurationParameter.value().get(), util::Duration("PT01H02M03S"));
  EXPECT(params.fruitParameter == Fruit::APPLE);
}

void testIncorrectValueOfFloatParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalFloatParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_opt_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalVariableParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_opt_variable_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalDateTimeParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_opt_date_time_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalDurationParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_opt_duration_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfEnumParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(),
                                       "error_in_fruit_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

class Parameters : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::Parameters";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Parameters/defaultValues") {
                      testDefaultValues();
                    });
    ts.emplace_back(CASE("ufo/Parameters/correctValues") {
                      testCorrectValues();
                    });
    ts.emplace_back(CASE("ufo/Parameters/incorrectValueOfFloatParameter") {
                      testIncorrectValueOfFloatParameter();
                    });
    ts.emplace_back(CASE("ufo/Parameters/incorrectValueOfOptionalFloatParameter") {
                      testIncorrectValueOfOptionalFloatParameter();
                    });
    ts.emplace_back(CASE("ufo/Parameters/incorrectValueOfOptionalVariableParameter") {
                      testIncorrectValueOfOptionalVariableParameter();
                    });
    // Conversion from string to DateTime or Duration calls abort() on failure,
    // so we can't test these cases. Leaving them commented-out in case this changes in future.
    //
    // ts.emplace_back(CASE("ufo/Parameters/incorrectValueOfOptionalDateTimeParameter") {
    //                   testIncorrectValueOfOptionalDateTimeParameter();
    //                 });
    // ts.emplace_back(CASE("ufo/Parameters/incorrectValueOfOptionalDurationParameter") {
    //                   testIncorrectValueOfOptionalDurationParameter();
    //                 });
    ts.emplace_back(CASE("ufo/Parameters/incorrectValueOfEnumParameter") {
                      testIncorrectValueOfEnumParameter();
                    });
  }
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PARAMETERS_H_
