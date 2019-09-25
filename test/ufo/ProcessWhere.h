/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_UFO_PROCESSWHERE_H_
#define TEST_UFO_PROCESSWHERE_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testProcessWhere() {
  oops::Log::info() << "weeeell" << std::endl;
  const eckit::LocalConfiguration conf = ::test::TestEnvironment::config();
  oops::Log::info() << conf;

  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));

  eckit::LocalConfiguration obsconf(conf, "ObsSpace");

  ioda::ObsSpace ospace(obsconf, bgn, end);
  ObsFilterData data(ospace);

  oops::Log::info() << "created ObsFilterData" << std::endl;

  const int nlocs = obsconf.getInt("nlocs");
  EXPECT(data.nlocs() == nlocs);

  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ProcessWhere", confs);
  oops::Log::info() << confs.size() << " confs" << std::endl;
  for (size_t jconf = 0; jconf < confs.size(); ++jconf) {
    std::vector<bool> result = processWhere(confs[jconf], data);
    const int size_ref = confs[jconf].getInt("size where true");
    const int size = std::count(result.begin(), result.end(), true);
    oops::Log::info() << "reference: " << size_ref << ", compare with " << size << std::endl;
    EXPECT(size == size_ref);
  }
}

// -----------------------------------------------------------------------------

class ProcessWhere : public oops::Test {
 public:
  ProcessWhere() {}
  virtual ~ProcessWhere() {}
 private:
  std::string testid() const {return "ufo::test::ProcessWhere";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ProcessWhere/testProcessWhere")
      { testProcessWhere(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PROCESSWHERE_H_
