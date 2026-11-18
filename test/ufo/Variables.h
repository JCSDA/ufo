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
    const std::vector<std::string> expectedDims =
        var.channels().empty() ? std::vector<std::string>{"Location"}
                               : std::vector<std::string>{"Location", "Channel"};
    EXPECT(var.dimList() == expectedDims);
    for (std::size_t jvar = 0; jvar < var.size(); ++jvar) {
      EXPECT(var.variable(jvar) == refvars[jvar]);
      EXPECT(var[jvar].dimList() == expectedDims);
    }
    // test the fullName() method
    const std::string refFullName = conf[jj].getString("reference full name");
    EXPECT_EQUAL(var.fullName(), refFullName);
  }
}

// -----------------------------------------------------------------------------
void testDimensionProvenance() {
  const std::vector<std::string> locationDims{"Location"};
  const std::vector<std::string> channelDims{"Location", "Channel"};
  const Variable parent("ObsValue/brightnessTemperature", std::vector<int>{3, 5});
  const Variable primitive = parent[1];
  EXPECT(parent.dimList() == channelDims);
  EXPECT_EQUAL(primitive.variable(), "brightnessTemperature_5");
  EXPECT_EQUAL(primitive.variable(0), "brightnessTemperature_5");
  EXPECT_EQUAL(primitive.size(), 1);
  EXPECT(primitive.channels().empty());
  EXPECT(primitive.dimList() == channelDims);
  EXPECT_EQUAL(primitive.toOopsObsVariables()[0], "brightnessTemperature_5");
  EXPECT_EQUAL(primitive.toOopsObsVariables().size(), 1);

  const Variable copied(primitive);
  Variable assigned("MetaData/latitude");
  assigned = primitive;
  const Variable regrouped(primitive, "DiagnosticFlags");
  EXPECT(copied.dimList() == channelDims);
  EXPECT(assigned.dimList() == channelDims);
  EXPECT(regrouped.dimList() == channelDims);
  EXPECT_EQUAL(regrouped.fullName(), "DiagnosticFlags/brightnessTemperature_5");
  EXPECT_EQUAL(primitive[0].variable(), "brightnessTemperature_5");
  EXPECT(primitive[0].dimList() == channelDims);
  EXPECT(regrouped[0].dimList() == channelDims);

  const Variable literal("ObsValue/brightnessTemperature_5");
  EXPECT(literal.dimList() == locationDims);
  EXPECT(literal[0].dimList() == locationDims);
  EXPECT(primitive != literal);
  assigned = literal;
  EXPECT(assigned.dimList() == locationDims);

  const Variable noChannels("ObsValue/temperature", std::vector<int>{});
  EXPECT(noChannels.dimList() == locationDims);
  EXPECT(noChannels[0].dimList() == locationDims);
  const Variable oneChannel("ObsValue/brightnessTemperature", std::vector<int>{5});
  EXPECT(oneChannel[0].dimList() == channelDims);

  const ufo::Variables mixed(std::vector<Variable>{literal, parent});
  EXPECT(mixed.variable(0).dimList() == locationDims);
  EXPECT(mixed.variable(1).dimList() == channelDims);
  EXPECT(mixed.variable(2).dimList() == channelDims);
  EXPECT_EQUAL(mixed.variable(1).variable(), "brightnessTemperature_3");
  EXPECT_EQUAL(mixed.variable(2).variable(), "brightnessTemperature_5");
  const ufo::Variables regroupedVars(mixed, "DiagnosticFlags");
  EXPECT(regroupedVars.variable(0).dimList() == locationDims);
  EXPECT(regroupedVars.variable(2).dimList() == channelDims);
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
// Test operator== for ufo::Variables and ufo::Variable
void testCompare() {
  ufo::Variable var1("MetaData/var1", {1, 3});
  ufo::Variable var2("MetaData/var1", {1, 2});
  ufo::Variable var3("MetaData/var2");
  ufo::Variable var4("ObsValue/var1", {1, 3});
  ufo::Variable var5(var1, "MetaData");
  EXPECT(var1 != var2);
  EXPECT(var1 != var3);
  EXPECT(var1 != var4);
  EXPECT(var1 == var5);

  ufo::Variables vars1({var1, var2, var3});
  ufo::Variables vars2({var2, var3, var1});
  ufo::Variables vars3({var1, var2});
  ufo::Variables vars4;
  ufo::Variables vars5 = vars3;
  EXPECT(vars1 == vars2);
  EXPECT(vars1 != vars3);
  EXPECT(vars4 != vars2);
  EXPECT(vars5 == vars3);
  vars5 += var3;
  EXPECT(vars5 != vars3);
  EXPECT(vars5 == vars1);

  EXPECT(var1[1] == var5[1]);
  EXPECT(var1[1] != var4[1]);
  EXPECT(Variable(var4[1], var1.group()) == var1[1]);
  EXPECT(var1[1][0] == var1[1]);

  const Variable literal("MetaData/var1_3");
  EXPECT(var1[1] != literal);
  EXPECT(!(literal == var1[1]));

  const ufo::Variables primitiveVars(std::vector<Variable>{var1[0], var1[1]});
  const ufo::Variables reorderedVars(std::vector<Variable>{var5[1], var5[0]});
  const ufo::Variables literalVars(std::vector<Variable>{var1[0], literal});
  EXPECT(primitiveVars == reorderedVars);
  EXPECT(primitiveVars != literalVars);
}

// -----------------------------------------------------------------------------
// Test ufo::Variables::toOopsVariables and ctor from oops::ObsVariables
void testToFromOopsObsVariables() {
  // single variable with channels
  ufo::Variable uVar("windEastward", {1, 2, 3});
  oops::ObsVariables uVar_oops = uVar.toOopsObsVariables();
  ufo::Variables ufoVars({uVar});
  oops::ObsVariables ufoVars_oops = ufoVars.toOopsObsVariables();
  oops::ObsVariables ufoVars_ref({"windEastward"}, {1, 2, 3});
  ufo::Variables ufoVars_from(ufoVars_ref);
  oops::Log::info() << "One-var channel ufo variable: " << uVar << std::endl;
  oops::Log::info() << "One-var channel ufo variables: " << ufoVars << std::endl;
  oops::Log::info() << "One-var channel ufo variable to oops vars: " <<
                       uVar.toOopsObsVariables() << std::endl;
  oops::Log::info() << "One-var channel ufo variables to oops vars: " <<
                       ufoVars.toOopsObsVariables() << std::endl;
  oops::Log::info() << "Expected oops vars: " << ufoVars_ref << std::endl;
  oops::Log::info() << "One-var channel ufo variables from oops vars: " <<
                       ufoVars_from << std::endl;
  EXPECT(ufoVars_oops == ufoVars_ref);
  EXPECT(ufoVars == ufoVars_from);

  ufo::Variable vVar("windNorthward", {1, 2, 3});
  ufoVars += vVar;
  oops::ObsVariables ufoVars_oops2 = ufoVars.toOopsObsVariables();
  oops::ObsVariables ufoVars_ref2({"windEastward", "windNorthward"}, {1, 2, 3});
  ufo::Variables ufoVars_from2(ufoVars_ref2);
  oops::Log::info() << "Two-var channel ufo variables: " << ufoVars << std::endl;
  oops::Log::info() << "Two-var channel ufo variables to oops vars: " <<
                       ufoVars_oops2 << std::endl;
  oops::Log::info() << "Expected oops vars: " << ufoVars_ref2 << std::endl;
  oops::Log::info() << "Two-var channel ufo variables from oops vars: " <<
                       ufoVars_from2 << std::endl;
  EXPECT(ufoVars_oops2 == ufoVars_ref2);
  EXPECT(ufoVars == ufoVars_from2);
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

    ts.emplace_back(CASE("ufo/Variables/testDimensionProvenance")
      { testDimensionProvenance(); });

    ts.emplace_back(CASE("ufo/Variables/testConstructor")
      { testConstructor(); });

    ts.emplace_back(CASE("ufo/Variables/testAllFromGroup")
      { testAllFromGroup(); });

    ts.emplace_back(CASE("ufo/Variables/testHasGroup")
      { testHasGroup(); });

    ts.emplace_back(CASE("ufo/Variables/testCompare")
      { testCompare(); });

    ts.emplace_back(CASE("ufo/Variables/testToFromOopsObsVariables")
      { testToFromOopsObsVariables(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_VARIABLES_H_
