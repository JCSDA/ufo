/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSBIASCOVARIANCEDETAILS_H_
#define TEST_UFO_OBSBIASCOVARIANCEDETAILS_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {
namespace test {
// -----------------------------------------------------------------------------

void testObsBiasCovarianceDetails() {
  eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  //  Setup ObsSpace
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));
  std::vector<eckit::LocalConfiguration> obsconfs
    = conf.getSubConfigurations("observations");

  for (auto & oconf : obsconfs) {
    const double tolerance = oconf.getDouble("tolerance");

    eckit::LocalConfiguration ospconf(oconf, "obs space");
    ioda::ObsSpace odb(ospconf, oops::mpi::world(), timeWindow, oops::mpi::myself());

    // Setup ObsBias
    eckit::LocalConfiguration biasconf = oconf.getSubConfiguration("obs bias");
    ObsBias ybias(odb, biasconf);

    // Setup ObsBiasIncrements
    eckit::LocalConfiguration biaserrconf = biasconf.getSubConfiguration("covariance");
    ObsBiasIncrement ybias_inc(odb, biasconf);

    // Setup ObsBiasCovariance (include reading from file)
    ObsBiasCovariance ybias_cov(odb, biasconf);

    // Randomize increments
    ybias_cov.randomize(ybias_inc);

    ObsBiasIncrement ybias_inc_2(ybias_inc);
    ObsBiasIncrement ybias_inc_3(ybias_inc);

    // linearize for first outer loop - expect to throw a message
    // since QC flags and effective errors weren't set yet
    biaserrconf.set("iteration", 0);
    EXPECT_THROWS(ybias_cov.linearize(ybias, biaserrconf));

    // mimic QC flags from first outer loop
    const std::vector<int> qc_flags(odb.nlocs(), 50);
    const std::vector<std::string> vars = odb.obsvariables().variables();
    for ( const auto & var : vars)
      odb.put_db("EffectiveQC0", var , qc_flags);

    // mimic effective errors
    const std::vector<float> errs(odb.nlocs(), 1.0);
    for ( const auto & var : vars)
      odb.put_db("EffectiveError0", var , errs);

    // mimic predictors
    ioda::ObsVector predx(odb);
    for (std::size_t jj = 0; jj < predx.size(); ++jj)
      predx[jj] = 1.0;
    for (const auto & pred : ybias_cov.predictorNames()) {
      predx.save(pred + "Predictor");
    }

    if (biasconf.has("covariance.output file")) {
      // Need a temporary workaround for issues that surface where multiple file handles
      // are pointing to the same file. This situation needs to be avoided, but will take
      // some refactoring to accomplish. For the mean time we can keep test refernce files
      // of the expected output and use those for the construct call below which
      // will keep the number of file handles at one per file.
      //
      // The idea is to name the test reference file the same as the output file
      // for the write call below, but keep the test reference file in
      // "Data/ufo/testinput_tier_1" while the output file lives in "Data".
      //
      // TODO(SRH) once the refactoring is completed for avoiding the "multiple file handles
      // pointing to the same file" issue, we can restore this code to directly read
      // the file created by the write call below.
      std::string output_file = biasconf.getString("covariance.output file");
      ybias_cov.write(biasconf);

      // replace the prefix "Data" with "Data/ufo/testinput_tier_1"
      auto pos = output_file.find("Data/");
      output_file.replace(pos, 4, "Data/ufo/testinput_tier_1");
      biasconf.set("covariance.prior.input file", output_file);
      ObsBiasCovariance ybias_cov2(odb, biasconf);
      ybias_cov.multiply(ybias_inc, ybias_inc_2);
      ybias_cov2.multiply(ybias_inc, ybias_inc_3);
      EXPECT(ybias_inc_2.norm() - ybias_inc_3.norm() < tolerance);
      oops::Log::test() << "ufo::testObsBiasCovarianceDetails read / write is verified"
                        << std::endl;
    }

    ybias_inc_2.zero();
    ybias_inc_3.zero();

    // Randomize increments again
    ybias_cov.randomize(ybias_inc);

    // linearize for the first outer loop again -- now QC and errors are set
    ybias_cov.linearize(ybias, biaserrconf);

    // delta_bias * B
    ybias_cov.multiply(ybias_inc, ybias_inc_2);

    EXPECT(ybias_inc.norm() != ybias_inc_2.norm());
    oops::Log::test() << "ufo::testObsBiasCovarianceDetails multiply is verified" << std::endl;

    // delta_bias / B
    ybias_cov.inverseMultiply(ybias_inc_2, ybias_inc_3);

    // Verifing the reading is right
    EXPECT(ybias_inc.norm() - ybias_inc_3.norm() < tolerance);
    oops::Log::test() << "ufo::testObsBiasCovarianceDetails inverseMultiply is verified"
                      << std::endl;
  }
}

// -----------------------------------------------------------------------------

class ObsBiasCovarianceDetails : public oops::Test {
 public:
  ObsBiasCovarianceDetails() {}
  virtual ~ObsBiasCovarianceDetails() {}
 private:
  std::string testid() const override {return "ufo::test::ObsBiasCovarianceDetails";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsBias/testObsBiasCovarianceDetails")
      { testObsBiasCovarianceDetails(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSBIASCOVARIANCEDETAILS_H_
