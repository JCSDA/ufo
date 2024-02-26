/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSBIASINCREMENT_H_
#define TEST_UFO_OBSBIASINCREMENT_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "test/TestEnvironment.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {
namespace test {
// -----------------------------------------------------------------------------

/// \brief Tests ObsBiasIncrement::read and ObsBiasIncrement::write methods
/// \details Test that if we write and then read or initialize ObsBiasIncrement the coefficients
/// are correct for the read/initialized ObsBiasIncrement.
void testObsBiasIncrementReadWrite() {
  eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  std::vector<eckit::LocalConfiguration> obsconfs = conf.getSubConfigurations("observations");

  for (auto & oconf : obsconfs) {
    ioda::ObsTopLevelParameters obsparams;
    obsparams.validateAndDeserialize(oconf.getSubConfiguration("obs space"));
    eckit::LocalConfiguration obsconf(oconf, "obs space");
    ioda::ObsSpace odb(obsconf, oops::mpi::world(), timeWindow, oops::mpi::myself());

    // set up ObsBias parameters
    eckit::LocalConfiguration biasconf = oconf.getSubConfiguration("obs bias");

    // set up modified configuration for reading in later
    eckit::LocalConfiguration biasconf_forread = biasconf;
    std::string readTestRefFile(biasconf.getString("increment output file"));
    biasconf_forread.set("test input file", readTestRefFile);

    // initialize ObsBiasIncrement by taking the difference between a non-trivial ObsBias and
    // a zero-valued ObsBias
    ObsBiasIncrement ybiasinc(odb, biasconf);
    ObsBias ybias1(odb, biasconf);
    ObsBias ybias2(odb, biasconf);
    ybias2.zero();
    ybiasinc.diff(ybias1, ybias2);

    // save original coefficients for comparison later
    const Eigen::VectorXd original_coeffs = ybiasinc.data();

    // write out; assign zero; read back in; compare with original coefficients
    eckit::LocalConfiguration conf;
    ybiasinc.write(conf);
    ybiasinc.zero();
    EXPECT_EQUAL(ybiasinc.data().norm(), 0.0);
    ybiasinc.read(biasconf_forread);
    EXPECT_EQUAL((original_coeffs - ybiasinc.data()).norm(), 0.0);
  }
}

// -----------------------------------------------------------------------------

class ObsBiasIncrement : public oops::Test {
 public:
  ObsBiasIncrement() = default;
  virtual ~ObsBiasIncrement() = default;

 private:
  std::string testid() const override {return "ufo::test::ObsBiasIncrement";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsBiasIncrement/testObsBiasIncrementReadWrite")
      { testObsBiasIncrementReadWrite(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSBIASINCREMENT_H_
