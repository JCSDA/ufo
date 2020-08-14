/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PARALLELOBSDISTRIBUTION_H_
#define TEST_UFO_PARALLELOBSDISTRIBUTION_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/ParallelObsDistribution.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  // TODO(wsmigaj) Report this (especially for floats) as a bug in eckit?
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::DateTime> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::Duration> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

template <typename T>
class TestParameters : public oops::Parameters {
 public:
  oops::RequiredParameter<Variable> variable{"variable", this};
  oops::RequiredParameter<std::vector<T>> expectedValues{"expectedValues", this};
};

template <typename T>
void ifTIsDoubleCastDoublesToFloats(std::vector<T> &v)
{}

template <>
void ifTIsDoubleCastDoublesToFloats(std::vector<double> &v) {
  for (double &x : v)
    x = static_cast<float>(x);
}

template <typename T>
void testVariable(const std::string &section) {
  const eckit::Configuration &topConf = ::test::TestEnvironment::config();

  util::DateTime bgn(topConf.getString("window begin"));
  util::DateTime end(topConf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(topConf, "obs space");
  ioda::ObsSpace obsSpace(obsSpaceConf, oops::mpi::comm(), bgn, end);

  TestParameters<T> parameters;
  parameters.deserialize(topConf.getSubConfiguration(section));
  std::vector<T> expectedValues = parameters.expectedValues;
  // IODA stores all floating-point variables in single precision, so if the T template
  // parameter is double, we need to cast the expected values, which have been loaded in double
  // precision from the YAML file, to single precision before comparing.
  ifTIsDoubleCastDoublesToFloats(expectedValues);

  ParallelObsDistribution obsDistribution(obsSpace);
  std::vector<T> globalValues = getGlobalVariableValues<T>(
        obsSpace, obsDistribution,
        parameters.variable.value().variable(), parameters.variable.value().group());

  EXPECT_EQUAL(globalValues, expectedValues);
}

CASE("ufo/ParallelObsDistribution/getGlobalIntVariableValues") {
  testVariable<int>("intTest");
}

CASE("ufo/ParallelObsDistribution/getGlobalFloatVariableValues") {
  testVariable<float>("floatOrDoubleTest");
}

CASE("ufo/ParallelObsDistribution/getGlobalDoubleVariableValues") {
  testVariable<double>("floatOrDoubleTest");
}

CASE("ufo/ParallelObsDistribution/getGlobalDateTimeVariableValues") {
  testVariable<util::DateTime>("dateTimeTest");
}

CASE("ufo/ParallelObsDistribution/members") {
  const eckit::Configuration &topConf = ::test::TestEnvironment::config();

  util::DateTime bgn(topConf.getString("window begin"));
  util::DateTime end(topConf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(topConf, "obs space");
  ioda::ObsSpace obsSpace(obsSpaceConf, oops::mpi::comm(), bgn, end);

  ParallelObsDistribution obsDistribution(obsSpace);
  const size_t gnlocs = obsSpace.gnlocs();
  EXPECT_EQUAL(obsDistribution.globalObsCount(), gnlocs);

  const size_t rank = obsSpace.comm().rank();
  const size_t nlocs = obsSpace.nlocs();
  EXPECT_EQUAL(obsDistribution.localObsCounts()[rank], nlocs);
}

class ParallelObsDistribution : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::ParallelObsDistribution";}

  void register_tests() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PARALLELOBSDISTRIBUTION_H_
