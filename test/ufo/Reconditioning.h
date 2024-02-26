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

#include "oops/base/ObsError.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"

#include "ufo/instantiateObsErrorFactory.h"

#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

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

  ufo::instantiateObsErrorFactory();
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    ufo::ObsErrorCrossVarCovParameters Params;
    Params.validateAndDeserialize(rconf);
    if (Params.reconditioning.value()->ReconMethod.value()
        != ufo::ReconditionMethod::NORECONDITIONING)
      continue;
    const eckit::LocalConfiguration obsSpaceConf(conf[jj], "obs space");
    ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::myself(), timeWindow, oops::mpi::myself());

    ObsErrorCrossVarCov R(Params, obsspace, oops::mpi::myself());
    ObsErrorCrossVarCov RRecon(Params, obsspace, oops::mpi::myself());
    oops::Log::info() << "Corr before:\n" << R << std::endl;

    ioda::ObsVector mask(obsspace, "ObsError");
    RRecon.recondition(Params, mask);
    oops::Log::info() << "Corr after:\n" << RRecon << std::endl;

    const double rmseR = R.getRMSE();
    const double rmseRR = RRecon.getRMSE();
    EXPECT(oops::is_close(rmseR, rmseRR, 1e-10));
  }
}

// -----------------------------------------------------------------------------
void compareKnownOutput() {
  eckit::LocalConfiguration Tconf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(Tconf.getSubConfiguration("time window"));

  ufo::instantiateObsErrorFactory();

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
    ObsErrorCrossVarCovParameters Params;
    ReconditioningTestParameters TestParams;
    Params.validateAndDeserialize(rconf);
    TestParams.validateAndDeserialize(testconf);

    ObsErrorCrossVarCov R(Params, obsspace, oops::mpi::myself());
    ioda::ObsVector mask(obsspace, "ObsError");
    ioda::ObsVector sample(obsspace, "ObsValue");
    std::vector<double> refVec = TestParams.refVec.value().value();
    R.recondition(Params, mask);
    R.multiply(sample);
    Eigen::VectorXd sampleVec = sample.packEigen(sample);
    oops::Log::info() << "R times sample vector: " << sampleVec.transpose() << std::endl;
    oops::Log::info() << "Reference vector: " << refVec << std::endl << std::endl;

    for (size_t i = 0; i < sampleVec.size(); ++i) {
      EXPECT(oops::is_close(sampleVec[i], refVec[i], 1e-5));
    }
  }
}

// -----------------------------------------------------------------------------
void testNoValidOptionSelected() {
  eckit::LocalConfiguration Tconf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(Tconf.getSubConfiguration("time window"));

  ufo::instantiateObsErrorFactory();

  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    if (!conf[jj].has("expectExceptionWithMessage")) {
      continue;
    }
    const eckit::LocalConfiguration obsSpaceConf(conf[jj], "obs space");
    ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::myself(), timeWindow, oops::mpi::myself());

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    ObsErrorCrossVarCovParameters Params;
    Params.validateAndDeserialize(rconf);

    const std::string msg = conf[jj].getString("expectExceptionWithMessage");
    EXPECT_THROWS_MSG(ObsErrorCrossVarCov R(Params, obsspace, oops::mpi::myself()), msg.c_str());
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
