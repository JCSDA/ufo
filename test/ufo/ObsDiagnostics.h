/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSDIAGNOSTICS_H_
#define TEST_UFO_OBSDIAGNOSTICS_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"

namespace eckit
{
  // Don't use the contracted output for this type: the current implementation works only
  // with integer types.
  // TODO(wsmigaj) Report this as a bug in eckit.
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testObsDiagnostics() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

  //  Setup ObsSpace
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));
  const eckit::LocalConfiguration obsconf(conf, "obs space");
  ioda::ObsTopLevelParameters obsparams;
  obsparams.validateAndDeserialize(obsconf);
  ioda::ObsSpace ospace(obsparams, oops::mpi::world(), bgn, end, oops::mpi::myself());
  const size_t nlocs = ospace.nlocs();

  // initialize observation operator (set variables requested from the model,
  // variables simulated by the observation operator, other init)
  eckit::LocalConfiguration obsopconf(conf, "obs operator");
  ObsOperatorParametersWrapper obsopparams;
  obsopparams.validateAndDeserialize(obsopconf);
  ObsOperator hop(ospace, obsopparams);

  // read geovals from the file
  eckit::LocalConfiguration gconf(conf, "geovals");
  GeoVaLsParameters geovalsparams;
  geovalsparams.validateAndDeserialize(gconf);
  const GeoVaLs gval(geovalsparams, ospace, hop.requiredVars());

  // initialize bias correction
  eckit::LocalConfiguration biasconf = conf.getSubConfiguration("obs bias");
  ObsBiasParameters biasparams;
  biasparams.validateAndDeserialize(biasconf);
  const ObsBias ybias(ospace, biasparams);

  // create obsvector to hold H(x)
  ioda::ObsVector hofx(ospace);

  // create obsvector to hold bias
  ioda::ObsVector bias(ospace);
  bias.zero();

  // create diagnostics to hold HofX diags
  eckit::LocalConfiguration diagconf(conf, "obs diagnostics");
  oops::Variables diagvars(diagconf, "variables");
  EXPECT(diagvars.size() > 0);
  std::unique_ptr<Locations> locs(hop.locations());
  ObsDiagnostics diags(ospace, *(locs.get()), diagvars);

  // call H(x) to compute diagnostics
  hop.simulateObs(gval, hofx, ybias, bias, diags);

  // read tolerance and reference Diagnostics
  const double tol = conf.getDouble("tolerance");
  eckit::LocalConfiguration diagrefconf(conf, "reference obs diagnostics");
  GeoVaLsParameters diagrefparams;
  diagrefparams.validateAndDeserialize(diagrefconf);
  ObsDiagnostics diagref(diagrefparams, ospace, diagvars);

  // loop over all diag variables and levels and compare with reference
  for (size_t ivar = 0; ivar < diagvars.size(); ivar++) {
    const size_t nlevs = diags.nlevs(diagvars[ivar]);
    EXPECT(nlevs == diagref.nlevs(diagvars[ivar]));
    for (size_t ilev = 0; ilev < nlevs; ilev++) {
      std::vector<float> ref(nlocs);
      std::vector<float> computed(nlocs);
      diags.get(computed, diagvars[ivar], ilev);
      diagref.get(ref, diagvars[ivar], ilev);

      float rms = 0.0;
      for (size_t iloc = 0; iloc < nlocs; iloc++) {
        ref[iloc] -= computed[iloc];
        rms += ref[iloc] * ref[iloc];
      }
      rms = sqrt(rms / nlocs);

      EXPECT(rms < 100*tol);
      oops::Log::info() << diagvars[ivar] << ", level " << ilev <<
          ": difference between reference and computed: " << ref << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

class ObsDiagnostics : public oops::Test {
 public:
  ObsDiagnostics() {}
  virtual ~ObsDiagnostics() {}
 private:
  std::string testid() const override {return "ufo::test::ObsDiagnostics";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsDiagnostics/testObsDiagnostics")
      { testObsDiagnostics(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSDIAGNOSTICS_H_
