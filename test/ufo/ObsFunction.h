/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_UFO_OBSFUNCTION_H_
#define TEST_UFO_OBSFUNCTION_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/GeoVaLs.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testFunction() {
  const eckit::LocalConfiguration conf = ::test::TestEnvironment::config();
///  Setup ObsSpace
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));
  const eckit::LocalConfiguration obsconf(conf, "ObsSpace");
  ioda::ObsSpace ospace(obsconf, bgn, end);

///  Get function name and which group to use for H(x)
  const eckit::LocalConfiguration obsfuncconf(conf, "ObsFunction");
  std::string funcname = obsfuncconf.getString("name");
  std::string grpname = obsfuncconf.getString("inputGroupName");

///  Setup function
  ObsFunction obsfunc(funcname);

///  Setup GeoVaLs
  oops::Variables geovars = obsfunc.requiredGeoVaLs();
  const eckit::LocalConfiguration gconf(conf, "GeoVaLs");
  const GeoVaLs gval(gconf, ospace, geovars);

///  Compute function result
  ioda::ObsDataVector<float> vals(ospace, funcname, "ObsFunction", false);
  ObsFilterData inputs(ospace);
  inputs.associate(gval);
  obsfunc.compute(inputs, vals);
  vals.save("TestResult");

///  Read reference values from ObsSpace
  std::string varref = obsfuncconf.getString("reference");
  ioda::ObsDataVector<float> ref(ospace, varref, "TestReference");

  const double tol = obsfuncconf.getDouble("tolerance");

///  Calculate rms(f(x) - ref) and compare to tolerance
  double zrms = 0.0;
  int nobs = 0;
  for (size_t jj = 0; jj < ref.nlocs() ; ++jj) {
    vals[0][jj] -= ref[0][jj];
    zrms += vals[0][jj] * vals[0][jj];
    nobs++;
  }
  ospace.comm().allReduceInPlace(zrms, eckit::mpi::sum());
  ospace.comm().allReduceInPlace(nobs, eckit::mpi::sum());
  if (nobs > 0) zrms = sqrt(zrms / static_cast<double>(nobs));
  oops::Log::info() << "Vector difference between reference and computed: " << vals << std::endl;
  EXPECT(zrms < 100*tol);  //  change tol from percent to actual value.
}

// -----------------------------------------------------------------------------

class ObsFunction : public oops::Test {
 public:
  ObsFunction() {}
  virtual ~ObsFunction() {}
 private:
  std::string testid() const {return "ufo::test::ObsFunction";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFunction/testFunction")
      { testFunction(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFUNCTION_H_
