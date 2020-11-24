/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSFUNCTION_H_
#define TEST_UFO_OBSFUNCTION_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void dataVectorDiff(const ioda::ObsSpace & ospace, ioda::ObsDataVector<float> & vals,
                    const ioda::ObsDataVector<float> & ref, std::vector<float> & rms_out) {
  /// Loop through variables and calculate rms for each variable
  for (size_t ivar = 0; ivar < vals.nvars() ; ++ivar) {
    float rms = 0.0;
    int nobs = 0;
    for (size_t jj = 0; jj < ref.nlocs() ; ++jj) {
      vals[ivar][jj] -= ref[ivar][jj];
      rms += vals[ivar][jj] * vals[ivar][jj];
      nobs++;
    }
    ospace.comm().allReduceInPlace(rms, eckit::mpi::sum());
    ospace.comm().allReduceInPlace(nobs, eckit::mpi::sum());
    if (nobs > 0) rms = sqrt(rms / static_cast<float>(nobs));
    rms_out[ivar] = rms;
  }
}

// -----------------------------------------------------------------------------

void testFunction() {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;

  std::vector<eckit::LocalConfiguration> typeconfs;
  ::test::TestEnvironment::config().get("observations", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ioda::ObsSpace &ospace = Test_::obspace()[jj].obsspace();
    const eckit::Configuration &conf = typeconfs[jj];

///  Setup ObsFilterData
    ObsFilterData inputs(ospace);

///  Get function name and which group to use for H(x)
    const eckit::LocalConfiguration obsfuncconf(conf, "obs function");
    Variable funcname(obsfuncconf);

///  Setup function
    ObsFunction obsfunc(funcname);
    ufo::Variables allfuncvars = obsfunc.requiredVariables();

///  Setup GeoVaLs
    const oops::Variables geovars = allfuncvars.allFromGroup("GeoVaLs").toOopsVariables();
    std::unique_ptr<GeoVaLs> gval;
    if (geovars.size() > 0) {
      const eckit::LocalConfiguration gconf(conf, "geovals");
      gval.reset(new GeoVaLs(gconf, ospace, geovars));
      inputs.associate(*gval);
    }

///  Setup ObsDiags
    const oops::Variables diagvars = allfuncvars.allFromGroup("ObsDiag").toOopsVariables();
    std::unique_ptr<ObsDiagnostics> diags;
    if (diagvars.size() > 0) {
      const eckit::LocalConfiguration diagconf(conf, "obs diagnostics");
      diags.reset(new ObsDiagnostics(diagconf, ospace, diagvars));
      inputs.associate(*diags);
    }

///  Get output variable names
    const oops::Variables outputvars(obsfuncconf, "variables");
///  Compute function result
    ioda::ObsDataVector<float> vals(ospace, outputvars, "ObsFunction", false);
    obsfunc.compute(inputs, vals);
    vals.save("TestResult");
    int nvars = vals.nvars();

///  Compute function result through ObsFilterData
    ioda::ObsDataVector<float> vals_ofd(ospace, outputvars, "ObsFunction", false);
    inputs.get(funcname, vals_ofd);

///  Read reference values from ObsSpace
    ioda::ObsDataVector<float> ref(ospace, outputvars, "TestReference");


    const double tol = obsfuncconf.getDouble("tolerance");

///  Calculate rms(f(x) - ref) and compare to tolerance
    std::vector<float> rms_out(nvars);
    dataVectorDiff(ospace, vals, ref, rms_out);

    oops::Log::info() << "Vector difference between reference and computed: " << std::endl;
    oops::Log::info() << vals << std::endl;
    for (size_t ivar = 0; ivar < nvars; ivar++) {
      EXPECT(rms_out[ivar] < 100*tol);  //  change tol from percent to actual value.
    }

    dataVectorDiff(ospace, vals_ofd, ref, rms_out);
    oops::Log::info() << "Vector difference between reference and computed via ObsFilterData: "
                      << std::endl;
    oops::Log::info() << vals_ofd << std::endl;
    for (size_t ivar = 0; ivar < nvars; ivar++) {
      EXPECT(rms_out[ivar] < 100*tol);  //  change tol from percent to actual value.
    }
  }
}

// -----------------------------------------------------------------------------

class ObsFunction : public oops::Test {
 public:
  ObsFunction() {}
  virtual ~ObsFunction() {}
 private:
  std::string testid() const override {return "ufo::test::ObsFunction";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFunction/testFunction")
      { testFunction(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFUNCTION_H_
