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
#include "oops/../test/TestEnvironment.h"
#include "oops/base/Variables.h"
#include "oops/runs/Test.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testVariable() {
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("Variables", conf);
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
  }
}

// -----------------------------------------------------------------------------
void testConstructor() {
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("oops variables", conf);
  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    // read variable from config
    oops::Variables oopsvars(conf[jj]);
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
// Test that ufo::Variables::allFromGroup() gets variables from the functions
void testAllFromGroup() {
  Variables vars;
  vars += Variable("height@GeoVaLs");
  vars += Variable("Velocity@ObsFunction");
  vars += Variable("latitude@MetaData");
  vars += Variable("temperature@ObsValue");
  vars += Variable("longitude@MetaData");

  Variables res = vars.allFromGroup("ObsValue");
  oops::Log::info() << res << std::endl;

  EXPECT(res.size() == 3);
  EXPECT(hasVariable(res, Variable("eastward_wind@ObsValue")));
  EXPECT(hasVariable(res, Variable("northward_wind@ObsValue")));
  EXPECT(hasVariable(res, Variable("temperature@ObsValue")));
}

// -----------------------------------------------------------------------------
// Test that ufo::Variables::hasGroup() works for functions
void testHasGroup() {
  ufo::Variables vars;
  vars += Variable("latitude@MetaData");
  vars += Variable("longitude@MetaData");

  EXPECT(vars.hasGroup("MetaData"));
  EXPECT(!vars.hasGroup("ObsValue"));

  vars += Variable("Velocity@ObsFunction");

  EXPECT(vars.hasGroup("ObsValue"));
}


// -----------------------------------------------------------------------------

class Variables : public oops::Test {
 public:
  Variables() {}
  virtual ~Variables() {}
 private:
  std::string testid() const {return "ufo::test::Variables";}

  void register_tests() const {
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
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_VARIABLES_H_
