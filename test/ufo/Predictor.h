/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PREDICTOR_H_
#define TEST_UFO_PREDICTOR_H_

#include <Eigen/Core>

#include <memory>
#include <set>
#include <string>
#include <utility>
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
#include "oops/util/Expect.h"
#include "oops/util/IntSetParser.h"
#include "test/interface/ObsTestsFixture.h"
#include "ufo/ObsTraits.h"


namespace ufo {
namespace test {

// -----------------------------------------------------------------------------
void testPredictor() {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
  typedef oops::ObsDiagnostics<ufo::ObsTraits>    ObsDiags_;
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  std::vector<eckit::LocalConfiguration> typeconfs;
  ::test::TestEnvironment::config().get("observations", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ioda::ObsSpace &ospace = Test_::obspace()[jj].obsspace();
    const eckit::Configuration &conf = typeconfs[jj];

    // initialize bias correction
    eckit::LocalConfiguration bconf(conf, "obs bias");
    const ObsBias ybias(ospace, bconf);
    // get predictor names
    std::vector<std::string> predictor_names = ybias.requiredPredictors();

    // Initialize GeoVaLs
    oops::Variables gvars;
    gvars += ybias.requiredVars();
    std::unique_ptr<GeoVaLs> gval;
    if (gvars.size() > 0) {
      // Read GeoVaLs from a file
      eckit::LocalConfiguration gconf(conf, "geovals");
      gval.reset(new GeoVaLs(gconf, ospace, gvars));
    } else {
      // Create an empty GeoVaLs object
      gval.reset(new GeoVaLs(ospace.distribution(), oops::Variables()));
    }
    gval->setDefaultFormat(GeoVaLFormat::REDUCED);

    // initialize Obs diagnostics
    oops::ObsVariables diagvars;
    diagvars += ybias.requiredHdiagnostics();
    std::vector<float> lons(ospace.nlocs());
    std::vector<float> lats(ospace.nlocs());
    std::vector<util::DateTime> times(ospace.nlocs());
    ospace.get_db("MetaData", "latitude", lats);
    ospace.get_db("MetaData", "longitude", lons);
    ospace.get_db("MetaData", "dateTime", times);
    auto locs = std::make_unique<SampledLocations>(lons, lats, times, ospace.distribution());
    ObsDiagnostics ydiags(ospace, SampledLocations_(std::move(locs)), diagvars);

    bool expect_error_message = false;

    // Calculate predictor values
    const std::size_t npreds = predictor_names.size();
    EXPECT(npreds > 0);
    std::vector<ioda::ObsVector> predData(npreds, ioda::ObsVector(ospace));

    const Predictors & predictors = ybias.predictors();
    for (std::size_t p = 0; p < npreds; ++p) {
      if (conf.has("expectExceptionWithMessage")) {
        const std::string msg = conf.getString("expectExceptionWithMessage");
        EXPECT_THROWS_MSG(predictors[p]->compute(ospace, *gval, ydiags, ybias,
                                                 predData[p]), msg.c_str());
        expect_error_message = true;
        break;
      }
      predictors[p]->compute(ospace, *gval, ydiags, ybias, predData[p]);
      predData[p].save(predictors[p]->name() + "Predictor");
    }

    if (expect_error_message) {
      continue;
    }

    // Read in tolerance from yaml
    const double tol = conf.getDouble("tolerance");

    // Get output variable names and read test data
    // Note prepend predictor_ to predictor names to distingush
    // predictor reference from obsbias reference
    std::vector<std::string> vars;
    for (std::size_t jp = 0; jp < npreds; ++jp) {
      vars.push_back("predictor_" + predictor_names[jp]);
    }
    const oops::ObsVariables testvars = ospace.assimvariables();
    const std::size_t nvars = testvars.size();
    if (!testvars.channels().empty()) {
      // At present we can label predictors with either the channel number or the variable
      // name, but not both. So make sure there's only one multi-channel variable.
      ASSERT(nvars == testvars.channels().size());
    }

    /// For each predictor for each variable compare computed predictor values to reference
    std::vector<double> testData(ospace.nlocs());
    for (std::size_t jp = 0; jp < npreds; ++jp) {
      const std::size_t nlocs = predData[jp].nlocs();
      for (std::size_t jv = 0; jv < nvars; ++jv) {
        std::string refVarName = "predictor_" + predictor_names[jp] + '_';
        if (testvars.channels().empty())
          refVarName += testvars[jv];
        else
          refVarName += std::to_string(testvars.channels()[jv]);

        ospace.get_db("TestReference", refVarName, testData);

        // compare test and reference vectors
        double rms = 0.0;
        for (size_t jl = 0; jl < nlocs; jl++) {
          testData[jl] -= predData[jp][jl*nvars + jv];
          rms += testData[jl] * testData[jl];
        }
        rms = std::sqrt(rms / nlocs);

        oops::Log::debug() << "Vector difference between reference and computed for "
                           << refVarName << ": " << testData << std::endl;
        EXPECT(rms <= tol);
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
