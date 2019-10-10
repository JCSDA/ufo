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
#include "oops/util/IntSetParser.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testConfigConstructor() {
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("Variables", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    // read variables from config
    ufo::Variables vars(conf[jj]);
    vars.removeDuplicates();
    // read reference vector of strings
    std::vector<std::string> refvars(conf[jj].getStringVector("reference"));
    std::sort(refvars.begin(), refvars.end());
    // compare the two
    EXPECT(vars.size() == refvars.size());
    for (std::size_t jvar = 0; jvar < vars.size(); ++jvar) {
      EXPECT(vars[jvar] == refvars[jvar]);
    }
  }
}

// -----------------------------------------------------------------------------

void testChannelsConstructor() {
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("Variables with channels", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    // read varname and channels from the config
    const std::string varname = conf[jj].getString("variable name");
    std::string chlist = conf[jj].getString("channels");
    std::set<int> channelset = oops::parseIntSet(chlist);
    std::vector<int> channels;
    std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels));
    ufo::Variables vars(varname, channels);
    // read reference vector of strings
    std::vector<std::string> refvars(conf[jj].getStringVector("reference"));
    std::sort(refvars.begin(), refvars.end());
    // compare the two
    EXPECT(vars.size() == refvars.size());
    for (std::size_t jvar = 0; jvar < vars.size(); ++jvar) {
      EXPECT(vars[jvar] == refvars[jvar]);
    }
  }
}

// -----------------------------------------------------------------------------
// Test that ufo::Variables::allFromGroup() gets variables from the functions
void testAllFromGroup() {
  ufo::Variables vars;
  vars += ("height@GeoVaLs");
  vars += ("Velocity@ObsFunction");
  vars += ("latitude@MetaData");
  vars += ("temperature@ObsValue");
  vars += ("longitude@MetaData");

  oops::Variables res = vars.allFromGroup("ObsValue");
  oops::Log::info() << res << std::endl;

  EXPECT(res.size() == 3);
  EXPECT(res.has("eastward_wind"));
  EXPECT(res.has("northward_wind"));
  EXPECT(res.has("temperature"));
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

    ts.emplace_back(CASE("ufo/Variables/testConfigConstructor")
      { testConfigConstructor(); });

    ts.emplace_back(CASE("ufo/Variables/testChannelsConstructor")
      { testChannelsConstructor(); });

    ts.emplace_back(CASE("ufo/Variables/testAllFromGroup")
      { testAllFromGroup(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_VARIABLES_H_
