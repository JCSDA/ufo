/*
 * (C) Copyright 2021- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSBIAS_H_
#define TEST_UFO_OBSBIAS_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/ObsBias.h"

namespace ufo {
namespace test {
// -----------------------------------------------------------------------------

/// \brief Tests ObsBias::read and ObsBias::write methods
/// \details Test that if we write and then read or initialize ObsBias the coefficients
/// are correct for the read/initialized ObsBias.
void testObsBiasReadWrite() {
  eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));
  std::vector<eckit::LocalConfiguration> obsconfs
    = conf.getSubConfigurations("observations");

  for (auto & oconf : obsconfs) {
    ioda::ObsSpace odb(oconf.getSubConfiguration("obs space"), oops::mpi::world(),
                       bgn, end, oops::mpi::myself());

    // setup ObsBias parameters
    eckit::LocalConfiguration biasconf = oconf.getSubConfiguration("obs bias");
    ObsBiasParameters biasparams;
    biasparams.validateAndDeserialize(biasconf);

    // setup ObsBias parameters with input file == output file from the original yaml
    eckit::LocalConfiguration biasconf_forread = biasconf;
    biasconf_forread.set("input file", biasconf.getString("output file"));
    ObsBiasParameters biasparams_forread;
    biasparams_forread.validateAndDeserialize(biasconf_forread);

    // init ObsBias
    ObsBias ybias(odb, biasparams);
    oops::Log::test() << "Initialized obs bias: " << ybias << std::endl;
    // save original coefficients for comparison later
    const Eigen::VectorXd original_coeffs = ybias.data();
    // write out; assign zero; read back in; compare with original coefficients
    ybias.write(biasparams);
    ybias.zero();
    EXPECT_EQUAL(ybias.data().norm(), 0.0);
    ybias.read(biasparams_forread);
    oops::Log::test() << "Obs bias after write out; assign zero; read back in: "
                      << ybias << std::endl;
    EXPECT_EQUAL((original_coeffs - ybias.data()).norm(), 0.0);
    // create ObsBias from output coefficients; compare with original coefficients
    ObsBias ybias2(odb, biasparams_forread);
    oops::Log::test() << "Obs bias initialized from the written out file: " << ybias2 << std::endl;
    EXPECT_EQUAL((original_coeffs - ybias2.data()).norm(), 0.0);
  }
}

// -----------------------------------------------------------------------------

class ObsBias : public oops::Test {
 public:
  ObsBias() = default;
  virtual ~ObsBias() = default;

 private:
  std::string testid() const override {return "ufo::test::ObsBias";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsBias/testObsBiasReadWrite")
      { testObsBiasReadWrite(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSBIAS_H_
