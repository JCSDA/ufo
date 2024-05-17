/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_VARIABLES_H_
#define TEST_UFO_VARIABLES_H_

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVariables.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testVariable() {
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("test variables", conf);
  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    // read variable from config
    Variable var(conf[jj]);
    // read reference vector of strings and group
    std::vector<std::string> refvars(conf[jj].getStringVector("reference names"));
    const std::string refgroup = conf[jj].getString("reference group");
    // compare the two
    EXPECT(var.size() == refvars.size());
    EXPECT(var.group() == refgroup);
    for (std::size_t jvar = 0; jvar < var.size(); ++jvar) {
      EXPECT(var.variable(jvar) == refvars[jvar]);
    }
    // test the fullName() method
    const std::string refFullName = conf[jj].getString("reference full name");
    EXPECT_EQUAL(var.fullName(), refFullName);
  }
}

// -----------------------------------------------------------------------------
void testConstructor() {
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("oops variables", conf);
  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    // read variable from config
    oops::ObsVariables oopsvars(conf[jj], "variables");
    // read reference vector of strings
    std::vector<std::string> refvars(conf[jj].getStringVector("reference names"));
    // init ufo::Variables
    Variables vars(oopsvars);
    // compare the two
    EXPECT(vars.nvars() == refvars.size());
    for (std::size_t jvar = 0; jvar < vars.nvars(); ++jvar) {
      EXPECT(vars.variable(jvar).variable() == refvars[jvar]);
    }
  }
}

// -----------------------------------------------------------------------------
// Helper function: checks if variable is in variables. Ignores functions,
// expands channels (i.e. for bt channels 1-2 bt_1, bt_2 are in vars, bt is not)
bool hasVariable(const Variables & vars, const Variable & var) {
  bool found = false;
  for (size_t jj = 0; jj < vars.nvars(); ++jj) {
    if (vars.variable(jj).variable() == var.variable() &&
        vars.variable(jj).group() == var.group()) found = true;
  }
  return found;
}


// -----------------------------------------------------------------------------
void testAllFromGroupFor(const std::string &funcName,
                         const eckit::LocalConfiguration &funcConf,
                         const ufo::Variables &initialVars) {
  ufo::Variables vars = initialVars;
  vars += Variable(funcName, funcConf);

  Variables res = vars.allFromGroup("ObsValue");
  oops::Log::info() << res << std::endl;

  EXPECT(res.size() == 2);
  EXPECT(hasVariable(res, Variable("eastward_wind@ObsValue")));
  EXPECT(hasVariable(res, Variable("temperature@ObsValue")));
}

// -----------------------------------------------------------------------------
// Test that ufo::Variables::allFromGroup() gets variables from the functions
void testAllFromGroup() {
  Variables vars;
  vars += Variable("GeoVaLs/height");
  vars += Variable("MetaData/latitude");
  vars += Variable("temperature@ObsValue");
  vars += Variable("MetaData/longitude");

  const eckit::Configuration &conf = ::test::TestEnvironment::config();
  testAllFromGroupFor("ObsFunction/Conditional",
                      eckit::LocalConfiguration(conf, "float conditional"),
                      vars);
  testAllFromGroupFor("IntObsFunction/Conditional",
                      eckit::LocalConfiguration(conf, "int conditional"),
                      vars);
  testAllFromGroupFor("StringObsFunction/Conditional",
                      eckit::LocalConfiguration(conf, "string conditional"),
                      vars);
  testAllFromGroupFor("DateTimeObsFunction/Conditional",
                      eckit::LocalConfiguration(conf, "datetime conditional"),
                      vars);
}

// -----------------------------------------------------------------------------
void testHasGroupFor(const std::string &funcName,
                     const eckit::LocalConfiguration &funcConf,
                     const ufo::Variables &initialVars) {
  ufo::Variables vars = initialVars;
  vars += Variable(funcName, funcConf);
  EXPECT(vars.hasGroup("ObsValue"));
}

// Test that ufo::Variables::hasGroup() works for functions
void testHasGroup() {
  ufo::Variables vars;
  vars += Variable("MetaData/latitude");
  vars += Variable("MetaData/longitude");

  EXPECT(vars.hasGroup("MetaData"));
  EXPECT(!vars.hasGroup("ObsValue"));

  const eckit::Configuration &conf = ::test::TestEnvironment::config();
  testHasGroupFor("ObsFunction/Conditional",
                  eckit::LocalConfiguration(conf, "float conditional"),
                  vars);
  testHasGroupFor("IntObsFunction/Conditional",
                  eckit::LocalConfiguration(conf, "int conditional"),
                  vars);
  testHasGroupFor("StringObsFunction/Conditional",
                  eckit::LocalConfiguration(conf, "string conditional"),
                  vars);
  testHasGroupFor("DateTimeObsFunction/Conditional",
                  eckit::LocalConfiguration(conf, "datetime conditional"),
                  vars);
}

// -----------------------------------------------------------------------------

class Variables : public oops::Test {
 public:
  Variables() {}
  virtual ~Variables() {}
 private:
  std::string testid() const override {return "ufo::test::Variables";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Variables/testVariable")
      { testVariable(); });

    ts.emplace_back(CASE("ufo/Variables/testConstructor")
      { testConstructor(); });

    ts.emplace_back(CASE("ufo/Variables/testAllFromGroup")
      { testAllFromGroup(); });

    ts.emplace_back(CASE("ufo/Variables/testHasGroup")
      { testHasGroup(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_VARIABLES_H_
