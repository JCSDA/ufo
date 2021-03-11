/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PREDICTOR_H_
#define TEST_UFO_PREDICTOR_H_

#include <memory>
#include <set>
#include <string>
#include <vector>


#define ECKIT_TESTING_SELF_REGISTER_CASES 0


#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsOperator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/IntSetParser.h"
#include "test/interface/ObsTestsFixture.h"
#include "ufo/ObsTraits.h"


namespace ufo {
namespace test {

// -----------------------------------------------------------------------------
void testPredictor() {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
  typedef oops::ObsDiagnostics<ufo::ObsTraits>    ObsDiags_;

  std::vector<eckit::LocalConfiguration> typeconfs;
  ::test::TestEnvironment::config().get("observations", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ioda::ObsSpace &ospace = Test_::obspace()[jj].obsspace();
    const eckit::Configuration &conf = typeconfs[jj];

    /// initialize bias correction
    eckit::LocalConfiguration bconf(conf, "obs bias");
    const ObsBias ybias(ospace, bconf);
    const std::set<int> jobs = oops::parseIntSet(bconf.getString("jobs"));
    int njobs = jobs.size();
    // get predictor names
    std::vector<std::string> predictor_names = ybias.requiredPredictors();

    /// read geovals from the file
    eckit::LocalConfiguration gconf(conf, "geovals");
    oops::Variables gvars;
    gvars += ybias.requiredVars();
    const GeoVaLs gval(gconf, ospace, gvars);

    /// read Obs diagnostics from file
    oops::Variables diagvars;
    diagvars += ybias.requiredHdiagnostics();
    std::vector<float> lons(ospace.nlocs());
    std::vector<float> lats(ospace.nlocs());
    std::vector<util::DateTime> times(ospace.nlocs());
    ospace.get_db("MetaData", "latitude", lats);
    ospace.get_db("MetaData", "longitude", lons);
    ospace.get_db("MetaData", "datetime", times);
    Locations locs(lons, lats, times, ospace.comm());
    ObsDiagnostics ydiags(ospace, locs, diagvars);

    /// Calculate predictor values
    const std::size_t npreds = predictor_names.size();
    std::vector<ioda::ObsVector> predData(npreds, ioda::ObsVector(ospace));

    const Predictors & predictors = ybias.predictors();
    for (std::size_t p = 0; p < npreds; ++p) {
      predictors[p]->compute(ospace, gval, ydiags, predData[p]);
      predData[p].save(predictors[p]->name() + "Predictor");
    }

    // Read in tolerance from yaml
    const double tol = conf.getDouble("tolerance");

    //  Get output variable names and read test data
    // Note prepend predictor_ to predictor names to distingush
    // predictor reference from obsbias reference
    std::vector<std::string> vars;
    for (std::size_t jp = 0; jp < npreds; ++jp) {
      vars.push_back("predictor_" + predictor_names[jp]);
    }
    std::vector<int> channels(jobs.begin(), jobs.end());
    const oops::Variables testvars(vars, channels);
    ioda::ObsDataVector<float> ref(ospace, testvars, "TestReference");

    /// For each predictor for each channel compare computed predictor values to reference
    for (std::size_t jp = 0; jp < npreds; ++jp) {
      const std::size_t nlocs  = predData[jp].nlocs();
      for (std::size_t jc = 0; jc < njobs; ++jc) {
        std::vector<float> TestData =
                ref["predictor_" + predictor_names[jp] + '_' + std::to_string(channels[jc])];

        // compare test and reference vectors
        double rms = 0.0;
        for (size_t jl = 0; jl < nlocs; jl++) {
          TestData[jl] -= predData[jp][jl*njobs + jc];
          rms += TestData[jl] * TestData[jl];
        }
        rms = std::sqrt(rms / nlocs);

        EXPECT(rms <= tol);
        oops::Log::debug() << "Vector difference between reference and computed for predictor_"
                             + predictor_names[jp] + '_' + std::to_string(channels[jc]) + ": "
                          << TestData << std::endl;
      }
    }
  }
}

class Predictor : public oops::Test {
 public:
  Predictor() {}
  virtual ~Predictor() {}
 private:
  std::string testid() const override {return "ufo::test::Predictor";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/Predictor/testPredictor")
      { testPredictor(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PREDICTOR_H_