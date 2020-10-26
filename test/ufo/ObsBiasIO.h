/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSBIASIO_H_
#define TEST_UFO_OBSBIASIO_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <iomanip>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"
#include "ufo/ObsBias.h"

namespace ufo {
namespace test {
// -----------------------------------------------------------------------------

void testObsBiasIO() {
  eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  //  Setup ObsSpace
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  ioda::ObsSpace odb(obsconf, oops::mpi::world(), bgn, end, oops::mpi::myself());
  const size_t nlocs = odb.nlocs();

  // Setup ObsBias (include reading from file)
  ObsBias ybias(odb, conf);

  // Verifing the reading is right
  double norm = ybias.norm();
  double norm_ref = conf.getDouble("obs bias.norm ref");
  double tolerance = conf.getDouble("obs bias.tolerance");

  oops::Log::debug() << "norm = " << std::setw(15) << std::setprecision(8)
                     << norm << std::endl;
  oops::Log::debug() << "norm_ref = " << std::setw(15) << std::setprecision(8)
                     << norm_ref << std::endl;
  EXPECT(std::abs(norm - norm_ref) < tolerance);
  oops::Log::test() << "ufo::testObsBiasIO Reading is verified" << std::endl;

  // Change the values in ObsBias
  const double any_number = -888.0;
  std::fill(ybias.begin(), ybias.end(), any_number);

  // Write out to file
  ybias.write(conf);

  // Read again for verification
  std::string output_filename = conf.getString("obs bias.analysis.dataout");
  conf.set("obs bias.prior.datain", output_filename);

  ybias.read(conf);

  // Verifing the reading is correct
  norm = ybias.norm() / std::sqrt(ybias.size());

  oops::Log::debug() << "norm = " << std::setw(15) << std::setprecision(8)
                     << norm << std::endl;
  oops::Log::debug() << "any_number = " << std::setw(15) << std::setprecision(8)
                     << std::abs(any_number) << std::endl;
  EXPECT(std::abs(norm - std::abs(any_number)) < tolerance);
  oops::Log::test() << "ufo::testObsBiasIO Writting is verified" << std::endl;
}

// -----------------------------------------------------------------------------

class ObsBiasIO : public oops::Test {
 public:
  ObsBiasIO() {}
  virtual ~ObsBiasIO() {}
 private:
  std::string testid() const override {return "ufo::test::ObsBiasIO";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsBias/testObsBiasIO")
      { testObsBiasIO(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSBIASIO_H_
