/*
 * (C) Crown copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_LINEAROBSBIASOPERATOR_H_
#define TEST_UFO_LINEAROBSBIASOPERATOR_H_

#include <iomanip>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Locations.h"
#include "oops/interface/SampledLocations.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsBiasOperator.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

// Tests whether application of the linear obs bias operator to an obs bias increment obtained
// by subtracting two obs biases from each other produces the expected results.
CASE("ufo/LinearObsBiasOperator/testLinearObsBiasOperator") {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  std::vector<eckit::LocalConfiguration> typeconfs;
  ::test::TestEnvironment::config().get("observations", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ioda::ObsSpace &odb = Test_::obspace()[jj].obsspace();
    const eckit::Configuration &conf = typeconfs[jj];

    // initialize obs bias objects
    eckit::LocalConfiguration biasConf = conf.getSubConfiguration("obs bias");
    eckit::LocalConfiguration targetBiasConf = conf.getSubConfiguration("target obs bias");
    const ObsBias bias(odb, biasConf);
    const ObsBias targetBias(odb, targetBiasConf);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf, "geovals");
    oops::Variables requiredVars = bias.requiredVars();
    requiredVars += targetBias.requiredVars();
    GeoVaLs geovals(gconf, odb, requiredVars);
    geovals.setDefaultFormat(GeoVaLFormat::REDUCED);

    // set up obs diagnostics
    oops::ObsVariables requiredHdiagnostics;
    requiredHdiagnostics += bias.requiredHdiagnostics();
    std::vector<float> lons(odb.nlocs());
    std::vector<float> lats(odb.nlocs());
    std::vector<util::DateTime> times(odb.nlocs());
    odb.get_db("MetaData", "latitude", lats);
    odb.get_db("MetaData", "longitude", lons);
    odb.get_db("MetaData", "dateTime", times);
    auto locs = std::make_unique<SampledLocations>(lons, lats, times, odb.distribution());
    ObsDiagnostics ydiags(odb, SampledLocations_(std::move(locs)), requiredHdiagnostics);

    // set TL trajectory to the geovals and the bias coeff. from the files
    LinearObsBiasOperator biasOperator(odb);
    biasOperator.setTrajectory(geovals, bias, ydiags);

    // set the bias increment to the difference between two obs biases initialised earlier
    ObsBiasIncrement biasInc(odb, biasConf);
    biasInc.diff(targetBias, bias);

    // apply the linear obs bias operator
    ioda::ObsVector dy(odb);
    biasOperator.computeObsBiasTL(biasInc, dy);

    // verify results
    const double dy_rms = dy.rms();
    const double expected_dy_rms = conf.getDouble("rms ref");
    const double tol = conf.getDouble("relative tolerance");
    EXPECT(oops::is_close_relative(dy_rms, expected_dy_rms, tol));
  }
}

// -----------------------------------------------------------------------------

class LinearObsBiasOperator : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::LinearObsBiasOperator";}

  void register_tests() const override {  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_LINEAROBSBIASOPERATOR_H_
