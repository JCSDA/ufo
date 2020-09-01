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
#include "oops/parallel/mpi/mpi.h"
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
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));
  std::vector<eckit::LocalConfiguration> obsconfs
    = conf.getSubConfigurations("observations");

  for (auto & oconf : obsconfs) {
    ioda::ObsSpace odb(oconf.getSubConfiguration("obs space"), oops::mpi::comm(), bgn, end);

    // Setup ObsBias
    ObsBias ybias(odb, oconf);

    // Setup ObsBiasIncrements
    ObsBiasIncrement ybias_inc(odb, oconf);
    ObsBiasIncrement ybias_inc_2(ybias_inc);
    ObsBiasIncrement ybias_inc_3(ybias_inc);
    ybias_inc_2.zero();
    ybias_inc_3.zero();

    // Setup ObsBiasCovariance (include reading from file)
    ObsBiasCovariance ybias_cov(odb, oconf);

    // Randomize increments
    ybias_cov.randomize(ybias_inc);

    // linearize for first outer loop
    oconf.set("iteration", 0);
    ybias_cov.linearize(ybias, oconf);

    // linearize for second outer loop
    oconf.set("iteration", 1);
    EXPECT_THROWS(ybias_cov.linearize(ybias, oconf));

    // mimic QC flags from first outer loop
    const std::vector<int> qc_flags(odb.nlocs(), 50);
    const std::vector<std::string> vars = odb.obsvariables().variables();
    for ( const auto & var : vars)
     odb.put_db("EffectiveQC0", var , qc_flags);

    // Randomize increments again
    ybias_cov.randomize(ybias_inc);

    // linearize for second outer loop
    ybias_cov.linearize(ybias, oconf);

    // delta_bias * B
    ybias_cov.multiply(ybias_inc, ybias_inc_2);

    EXPECT(ybias_inc.norm() != ybias_inc_2.norm());
    oops::Log::test() << "ufo::testObsBiasCovarianceDetails multiply is verified" << std::endl;

    // delta_bias / B
    ybias_cov.inverseMultiply(ybias_inc_2, ybias_inc_3);

    // Verifing the reading is right
    const double tolerance = oconf.getDouble("tolerance");
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
  std::string testid() const {return "ufo::test::ObsBiasCovarianceDetails";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsBias/testObsBiasCovarianceDetails")
      { testObsBiasCovarianceDetails(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSBIASCOVARIANCEDETAILS_H_
