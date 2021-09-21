/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_METOFFICESORT_H_
#define TEST_UFO_METOFFICESORT_H_

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/utils/metoffice/MetOfficeSort.h"

namespace ufo {
namespace test {

CASE("ufo/MetOfficeSort/noKey") {
  const eckit::Configuration &topLevelConf = ::test::TestEnvironment::config();
  for (const eckit::LocalConfiguration &conf : topLevelConf.getSubConfigurations("no key")) {
    std::vector<int> input = conf.getIntVector("input");
    const std::vector<int> expectedOutput = conf.getIntVector("output");
    metOfficeSort(input.begin(), input.end());
    EXPECT_EQUAL(input, expectedOutput);
  }
}

CASE("ufo/MetOfficeSort/withKey") {
  const eckit::Configuration &topLevelConf = ::test::TestEnvironment::config();
  for (const eckit::LocalConfiguration &conf : topLevelConf.getSubConfigurations("with key")) {
    const std::vector<std::string> keys = conf.getStringVector("keys");
    std::vector<size_t> index(keys.size());
    std::iota(index.begin(), index.end(), 0);
    const std::vector<size_t> expectedOutput = conf.getUnsignedVector("output");
    metOfficeSort(index.begin(), index.end(), [&keys] (size_t i) { return keys[i]; } );
    EXPECT_EQUAL(index, expectedOutput);
  }
}

class MetOfficeSort : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::MetOfficeSort";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_METOFFICESORT_H_
