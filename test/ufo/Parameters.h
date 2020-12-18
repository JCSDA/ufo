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
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "test/TestEnvironment.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {
namespace test {

const bool validationSupported = oops::Parameters::isValidationSupported();

class MyParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MyParameters, Parameters)
 public:
  oops::OptionalParameter<ufo::Variable> optVariableParameter{"opt_variable_parameter", this};
};

template <typename ParametersType>
void doTestSerialization(const eckit::Configuration &config) {
  // We deserialize a configuration loaded from a YAML file into parameters and then serialize them
  // back into a configuration. The test verifies that the configuration objects produce the same
  // output when printed.
  //
  // For this to work, parameter names in the YAML file must be ordered alphabetically; that's
  // because the YAML parser creates configurations storing keys and values OrderedMapContent
  // objects (preserving the order in which individual options were specified in the YAML file),
  // but the LocalConfiguration::set() method stores keys and values in MapContent objects (with
  // keys ordered alphabetically).

  ParametersType params;
  params.deserialize(config);

  eckit::LocalConfiguration outputConfig;
  params.serialize(outputConfig);

  std::stringstream expectedStream;
  expectedStream << config;
  const std::string expected = expectedStream.str();

  std::stringstream receivedStream;
  receivedStream << outputConfig;
  std::string received = receivedStream.str();

  EXPECT_EQUAL(received, expected);
}

void testDefaultValue() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  MyParameters params;
  EXPECT(params.optVariableParameter.value() == boost::none);

  const eckit::LocalConfiguration emptyConf(conf, "empty");
  EXPECT_NO_THROW(params.validate(emptyConf));
  params.deserialize(emptyConf);

  EXPECT(params.optVariableParameter.value() == boost::none);
}

void testCorrectValue() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(::test::TestEnvironment::config(), "full");
  EXPECT_NO_THROW(params.validate(fullConf));
  params.deserialize(fullConf);

  EXPECT(params.optVariableParameter.value() != boost::none);
  EXPECT_EQUAL(params.optVariableParameter.value().get().group(), "MetaData");
  EXPECT_EQUAL(params.optVariableParameter.value().get().variable(), "latitude");
  EXPECT_EQUAL(params.optVariableParameter.value().get().channels(), (std::vector<int>{1, 5}));
}

void testSimpleString() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(::test::TestEnvironment::config(), "simple_string");
  EXPECT_NO_THROW(params.validate(fullConf));
  params.deserialize(fullConf);

  EXPECT(params.optVariableParameter.value() != boost::none);
  EXPECT_EQUAL(params.optVariableParameter.value().get().group(), "MetaData");
  EXPECT_EQUAL(params.optVariableParameter.value().get().variable(), "latitude");
}

void testNoChannels() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(::test::TestEnvironment::config(), "no_channels");
  EXPECT_NO_THROW(params.validate(fullConf));
  params.deserialize(fullConf);

  EXPECT(params.optVariableParameter.value() != boost::none);
  EXPECT_EQUAL(params.optVariableParameter.value().get().group(), "MetaData");
  EXPECT_EQUAL(params.optVariableParameter.value().get().variable(), "latitude");
  EXPECT_EQUAL(params.optVariableParameter.value().get().channels(), (std::vector<int>{}));
}

void testComplexChannels() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(::test::TestEnvironment::config(), "complex_channels");
  EXPECT_NO_THROW(params.validate(fullConf));
  params.deserialize(fullConf);

  EXPECT(params.optVariableParameter.value() != boost::none);
  EXPECT_EQUAL(params.optVariableParameter.value().get().group(), "MetaData");
  EXPECT_EQUAL(params.optVariableParameter.value().get().variable(), "latitude");
  EXPECT_EQUAL(params.optVariableParameter.value().get().channels(),
               (std::vector<int>{1, 5, 10, 11, 12, 13, 14, 15}));
}

void testMissingName() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(), "missing_name");
  if (validationSupported)
    EXPECT_THROWS(params.validate(conf));
  EXPECT_THROWS_MSG(params.deserialize(conf), "ConfigurationNotFound: [name]");
}

void testMisspelledProperty() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(), "misspelled_property");
  if (validationSupported)
    EXPECT_THROWS(params.validate(conf));
}

void testSerializationWithChannels() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(), "full");
  doTestSerialization<MyParameters>(conf);
}

void testSerializationWithoutChannels() {
  MyParameters params;
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config(), "no_channels");
  doTestSerialization<MyParameters>(conf);
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
    ts.emplace_back(CASE("ufo/Parameters/testSimpleString") {
                      testSimpleString();
                    });
    ts.emplace_back(CASE("ufo/Parameters/noChannels") {
                      testNoChannels();
                    });
    ts.emplace_back(CASE("ufo/Parameters/complexChannels") {
                      testComplexChannels();
                    });
    ts.emplace_back(CASE("ufo/Parameters/missingName") {
                      testMissingName();
                    });
    ts.emplace_back(CASE("ufo/Parameters/misspelledProperty") {
                      testMisspelledProperty();
                    });
    ts.emplace_back(CASE("ufo/Parameters/serializationWithChannels") {
                      testSerializationWithChannels();
                    });
    ts.emplace_back(CASE("ufo/Parameters/serializationWithoutChannels") {
                      testSerializationWithoutChannels();
                    });
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PARAMETERS_H_
