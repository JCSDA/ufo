/*
 * (C) Crown copyright 2023, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_RECONDITIONING_H_
#define TEST_UFO_RECONDITIONING_H_

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
#include "oops/util/Logger.h"

#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

#include "ufo/errors/ObsErrorBase.h"
#include "ufo/errors/ObsErrorCrossVarCov.h"
#include "ufo/errors/ObsErrorDiagonal.h"
#include "ufo/errors/ObsErrorWithinGroupCov.h"


namespace ufo {
namespace test {

class ReconditioningTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ReconditioningTestParameters, Parameters)
 public:
  oops::Parameter<bool> testReader{"test reader", false, this};
  oops::OptionalParameter<std::vector<double>> refVec{"reference", this};
};

// -----------------------------------------------------------------------------
void testNoReconditioning() {
  eckit::LocalConfiguration Tconf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(Tconf.getSubConfiguration("time window"));

  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    ufo::ObsErrorParametersWrapper Params;
    Params.validateAndDeserialize(rconf);
    if (Params.errorParams().reconditioning.value().ReconMethod.value()
        != ufo::ObsErrorReconditionerMethod::NORECONDITIONING)
      continue;
    const eckit::LocalConfiguration obsSpaceConf(conf[jj], "obs space");
    ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::myself(), timeWindow, oops::mpi::myself());

    ObsErrorBase* R = ufo::ObsErrorFactory::create(Params.errorParams(), obsspace);
    ObsErrorBase* RRecon = ufo::ObsErrorFactory::create(Params.errorParams(), obsspace);
    oops::Log::info() << "Corr before:\n" << R << std::endl;

    ioda::ObsVector mask(obsspace, "ObsError");
    RRecon->update(mask);
    oops::Log::info() << "Corr after:\n" << RRecon << std::endl;

    const double rmseR = R->getRMSE();
    const double rmseRR = RRecon->getRMSE();
    EXPECT(oops::is_close(rmseR, rmseRR, 1e-10));
  }
}

// -----------------------------------------------------------------------------
void compareKnownOutput() {
  eckit::LocalConfiguration Tconf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(Tconf.getSubConfiguration("time window"));

  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    if (!conf[jj].has("obs error test")) {
      continue;
    }
    const eckit::LocalConfiguration obsSpaceConf(conf[jj], "obs space");
    ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::myself(), timeWindow, oops::mpi::myself());

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    const eckit::LocalConfiguration testconf(conf[jj], "obs error test");
    ufo::ObsErrorParametersWrapper Params;
    ReconditioningTestParameters TestParams;
    Params.validateAndDeserialize(rconf);
    TestParams.validateAndDeserialize(testconf);

    ObsErrorBase* R = ufo::ObsErrorFactory::create(Params.errorParams(), obsspace);
    ioda::ObsVector mask(obsspace, "ObsError");
    ioda::ObsVector sample(obsspace, "ObsValue");
    std::vector<double> refVec = TestParams.refVec.value().value();
    // update method reconditions the ObsError matrix
    R->update(mask);
    // multiply method applies reconditioned ObsError matrix to the sample vector
    R->multiply(sample);
    std::vector<double> sampleVec;
    sample.maskAndSerialize(sample, sampleVec);
    oops::Log::info() << "R times sample vector: " << sampleVec << std::endl;
    oops::Log::info() << "Reference vector: " << refVec << std::endl << std::endl;
    ASSERT(sampleVec.size() == refVec.size());

    for (size_t i = 0; i < sampleVec.size(); ++i) {
      EXPECT(oops::is_close(sampleVec[i], refVec[i], 1e-5));
    }
  }
}

// -----------------------------------------------------------------------------
void testNoValidOptionSelected() {
  eckit::LocalConfiguration Tconf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(Tconf.getSubConfiguration("time window"));

  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    if (!conf[jj].has("expectExceptionWithMessage")) {
      continue;
    }
    const eckit::LocalConfiguration obsSpaceConf(conf[jj], "obs space");
    ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::myself(), timeWindow, oops::mpi::myself());

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    ufo::ObsErrorParametersWrapper Params;
    Params.validateAndDeserialize(rconf);

    const std::string msg = conf[jj].getString("expectExceptionWithMessage");
    EXPECT_THROWS_MSG(ufo::ObsErrorFactory::create(Params.errorParams(), obsspace), msg.c_str());
  }
}

// -----------------------------------------------------------------------------
class Reconditioning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::Reconditioning";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Recondtioning/") {
                      testNoReconditioning();
                    });
    ts.emplace_back(CASE("ufo/Recondtioning/") {
                      compareKnownOutput();
                    });
    ts.emplace_back(CASE("ufo/Recondtioning/") {
                      testNoValidOptionSelected();
                    });
  }
  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_RECONDITIONING_H_
