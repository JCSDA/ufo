/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_DIFFUSIONDIRAC_H_
#define TEST_UFO_DIFFUSIONDIRAC_H_

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"

#include "ufo/errors/ObsErrorDiffusion.h"
#include "ufo/errors/ObsErrorWithinGroupCov.h"

#include "ufo/ObsTraits.h"

#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace ufo {
namespace test {

using ObsSpace_ = ioda::ObsSpace;

// Dirac test for ObsErrorDiffusion: place a unit impulse at one obs, apply R, and log
// the propagated covariance at the peak and two neighbors. Runs on the world
// communicator to exercise the multi-PE path; read-outs are keyed by global obs index
// so the result is identical on any PE count.
void diracDiffusion() {
  // get time window
  eckit::LocalConfiguration Tconf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(Tconf.getSubConfiguration("time window"));

  // get obs configurations
  std::vector<eckit::LocalConfiguration> conf;
  ::test::TestEnvironment::config().get("observations", conf);
  eckit::LocalConfiguration testConf;
  ::test::TestEnvironment::config().get("test settings", testConf);

  oops::Log::info() << " CONFIGURATION [0]:" << std::endl;
  oops::Log::info() << conf[0] << std::endl;

  // get obsSpace confs
  eckit::LocalConfiguration obsSpaceDiff(conf[0], "obs space");

  // World communicator so obs are distributed across ranks and ObsErrorDiffusion runs
  // its multi-PE path (the time communicator -- 4th arg -- stays myself()).
  ioda::ObsSpace odb_geo_diff(obsSpaceDiff, oops::mpi::world(), timeWindow, oops::mpi::myself());

  int nlocs = odb_geo_diff.nlocs();

  oops::Log::info() << " --> number of locations: " << nlocs << std::endl;

  // Source-file obs index: partition-invariant, so the dirac location and diagnostic
  // read-outs are PE-agnostic (identical on 1 PE and N PE).
  const std::vector<std::size_t> & srcIndex = odb_geo_diff.index();
  const eckit::mpi::Comm & comm = odb_geo_diff.comm();

  // Value of dy at global obs index gLoc, on any rank: exactly one rank owns it
  // (non-overlapping dist), so a SUM over (value on owner, 0 elsewhere) broadcasts it.
  auto globalValueAt = [&](const ioda::ObsVector & v, std::size_t gLoc) -> double {
    double val = 0.0;
    for (int i = 0; i < nlocs; ++i)
      if (srcIndex[i] == gLoc) val = v[i];
    comm.allReduceInPlace(val, eckit::mpi::Operation::SUM);
    return val;
  };

  std::vector<float> lons(nlocs);
  std::vector<float> lats(nlocs);

  odb_geo_diff.get_db("MetaData", "longitude", lons);
  odb_geo_diff.get_db("MetaData", "latitude", lats);

  // get obsError conf
  eckit::LocalConfiguration obsErrorDiff(conf[0], "obs error");

  ObsErrorDiffusionParameters diffParams;

  diffParams.validateAndDeserialize(obsErrorDiff);

  oops::Log::info() << " ObsErr Conf [0]:" << std::endl;
  oops::Log::info() << obsErrorDiff << std::endl;

  ioda::ObsVector dy(odb_geo_diff, "ObsError");
  // create diffusion obsError
  ObsErrorDiffusion R_diff(diffParams, odb_geo_diff, oops::mpi::myself());
  R_diff.update(dy);

  R_diff.randomize(dy);
  ioda::ObsVector dy_diff(dy);

  // `dirac location` is a global obs index; set the impulse on whichever rank owns
  // that obs (same physical obs at any PE count).
  int location = testConf.getInt("dirac location");
  if (testConf.getString("test type") == "dirac") {
    dy_diff.zero();
    for (int i = 0; i < nlocs; ++i)
      if (srcIndex[i] == static_cast<std::size_t>(location)) dy_diff[i] = 1;
  }

  oops::Log::info() << "Initial diffusion dy_diff:\n" << dy_diff.data() << std::endl;

  R_diff.multiply(dy_diff);

  oops::Log::info() << "Final diffusion dy:\n" << dy_diff.data() << std::endl;
  oops::Log::info() << "Diffusion Dirac R*dy stats:" << std::endl << dy_diff << std::endl;
  oops::Log::test() << "Diffusion Dirac R*dy stats:" << std::endl << dy_diff << std::endl;
  // Read out by GLOBAL obs index so the values are identical on 1 PE and N PE.
  oops::Log::test() << "Final peak of dirac impulse: "
                    << globalValueAt(dy_diff, location) << std::endl;
  // print out propagated covariance at the two grid locations closest to the peak
  oops::Log::test() << "Propagated covariance at location 38: "
                    << globalValueAt(dy_diff, 38) << std::endl;
  oops::Log::test() << "Propagated covariance at location 51: "
                    << globalValueAt(dy_diff, 51) << std::endl;
}

class DiffusionDirac : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::Diffusion";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/DiffusionDirac/") {
                      diracDiffusion();
                    });
  }
  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_DIFFUSIONDIRAC_H_
