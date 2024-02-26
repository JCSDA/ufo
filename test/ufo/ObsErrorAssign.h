/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSERRORASSIGN_H_
#define TEST_UFO_OBSERRORASSIGN_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/BlackList.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testObsErrorAssign(const eckit::LocalConfiguration &conf) {
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  std::vector<std::string> varnames {"pressure", "stationIdentification", "observationTypeNum",
                                     "latitudeBand", "sensorChannelNumber", "dataProviderOrigin",
                                     "satelliteIdentifier"};
  for (std::string varname : varnames) {
    if (conf.has(varname)) {
      const std::vector<int> dat = conf.getIntVector(varname);
      obsspace.put_db("MetaData", varname, dat);
    }
  }
  if (conf.has("stringVar")) {
      const std::vector<std::string> dat = conf.getStringVector("stringVar");
      obsspace.put_db("MetaData", "stringVar", dat);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "ObsError assign");
  ufo::BlackListParameters bfparam;
  bfparam.deserialize(filterConf);

  ufo::BlackList filter(obsspace, bfparam, qcflags, obserr);
  if (conf.has("expectExceptionWithMessage")) {
    const std::string msg = conf.getString("expectExceptionWithMessage");
    EXPECT_THROWS_MSG(filter.preProcess(), msg.c_str());
    return;
  } else {
    filter.preProcess();
  }

  std::vector<float> expectedObsError;
  if (conf.has("expected_obserror")) {
    expectedObsError = conf.getFloatVector("expected_obserror");
  } else {
    expectedObsError = conf.getFloatVector("expected_obserror_variance");
    for (size_t ind=0; ind < expectedObsError.size(); ind++) {
      expectedObsError[ind] = std::sqrt(expectedObsError[ind]);
    }
  }
  int ind = 0;
  for (size_t varn = 0; varn < obserr->nvars(); ++varn) {
    for (size_t locn = 0; locn < obserr->nlocs(); ++locn) {
      EXPECT(oops::is_close_absolute((*obserr)[varn][locn], expectedObsError[ind], 1e-4f, 0,
                                     oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE));
      ind++;
    }
  }
}

class ObsErrorAssign : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::ObsErrorAssign";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/ObsErrorAssign/" + testCaseName, testCaseConf)
                      {
                        testObsErrorAssign(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSERRORASSIGN_H_
