/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PRIMITIVEVARIABLES_H_
#define TEST_UFO_PRIMITIVEVARIABLES_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/utils/PrimitiveVariables.h"

namespace eckit
{
  // Don't use the contracted output for floats: the current implementation works only
  // with integer types.
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

class TestCaseParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestCaseParameters, Parameters);
 public:
  oops::RequiredParameter<std::vector<Variable>> variables{"variables", this};
  oops::RequiredParameter<std::vector<std::string>> expectedNames{"expected names", this};
  oops::RequiredParameter<std::vector<std::string>> expectedGroups{"expected groups", this};
  oops::RequiredParameter<std::vector<std::vector<float>>> expectedValues{"expected values", this};
};

/// Code shared by all tests
class TestFixture : private boost::noncopyable {
 public:
  static const ioda::ObsSpace & obsspace() { return *getInstance().obsspace_; }
  static const ObsFilterData & data() { return *getInstance().data_; }

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
    data_.reset(new ObsFilterData(*obsspace_));
  }

  std::shared_ptr<ioda::ObsSpace> obsspace_;
  std::unique_ptr<ObsFilterData> data_;
};

/// Test a range-based for loop over the PrimitiveVariables range.
void testPrimitiveVariables(const eckit::LocalConfiguration &conf) {
  TestCaseParameters parameters;
  parameters.validateAndDeserialize(conf);

  Variables variables;
  for (const Variable &variable : parameters.variables.value())
    variables += variable;

  std::vector<std::string> names;
  std::vector<std::string> groups;
  std::vector<Variable> primitiveVariables;
  std::vector<std::vector<float>> values;

  for (PrimitiveVariable pv : PrimitiveVariables(variables, TestFixture::data())) {
    names.push_back(pv.name());
    groups.push_back(pv.group());
    primitiveVariables.push_back(pv.variable());
    values.push_back(pv.values());
  }

  EXPECT_EQUAL(names, parameters.expectedNames.value());
  EXPECT_EQUAL(groups, parameters.expectedGroups.value());
  EXPECT_EQUAL(values.size(), parameters.expectedValues.value().size());
  for (size_t i = 0; i < values.size(); ++i)
    EXPECT(oops::are_all_close_relative(values[i], parameters.expectedValues.value()[i], 1e-6f));

  EXPECT(std::equal(primitiveVariables.begin(), primitiveVariables.end(), names.begin(),
                    [](const Variable &primitiveVariable, const std::string &name)
                    {return primitiveVariable.variable() == name; }));
  EXPECT(std::equal(primitiveVariables.begin(), primitiveVariables.end(), groups.begin(),
                    [](const Variable &primitiveVariable, const std::string &group)
                    {return primitiveVariable.group() == group; }));
}

/// Test explicitly that the * and -> operator overloads in PrimitiveVariablesIterator work
/// correctly.
void testPrimitiveVariablesIterator(const eckit::LocalConfiguration &conf) {
  TestCaseParameters parameters;
  parameters.validateAndDeserialize(conf);

  Variables variables;
  for (const Variable &variable : parameters.variables.value())
    variables += variable;

  std::vector<std::string> names;
  std::vector<std::vector<float>> values;

  {
    PrimitiveVariables pvars(variables, TestFixture::data());
    for (PrimitiveVariablesIterator it = pvars.begin(); it != pvars.end(); ++it) {
      // Use the * operator to get access to the variable name and the -> operator to get access
      // to its values.
      names.push_back((*it).name());
      values.push_back(it->values());
    }
  }

  EXPECT_EQUAL(names, parameters.expectedNames.value());
  EXPECT_EQUAL(values.size(), parameters.expectedValues.value().size());
  for (size_t i = 0; i < values.size(); ++i)
    EXPECT(oops::are_all_close_relative(values[i], parameters.expectedValues.value()[i], 1e-6f));
}

class PrimitiveVariables : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::PrimitiveVariables";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config(), "cases");
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(conf, testCaseName);
      ts.emplace_back(CASE("ufo/PrimitiveVariables/range/" + testCaseName, testCaseConf)
                      {
                        testPrimitiveVariables(testCaseConf);
                      });
      ts.emplace_back(CASE("ufo/PrimitiveVariables/iterator/" + testCaseName, testCaseConf)
                      {
                        testPrimitiveVariablesIterator(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PRIMITIVEVARIABLES_H_
