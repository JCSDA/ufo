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
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Locations.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"
#include "ufo/ObsTraits.h"

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

class ObsDiagnosticsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsDiagnosticsParameters, Parameters)

 public:
  /// \brief List of obs diagnostics to be calculated.
  oops::RequiredParameter<oops::ObsVariables> variables{"variables", this};
};

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ObsDiagnostics test.
class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> timeWindow{"time window", this};
  oops::RequiredParameter<eckit::LocalConfiguration> obsSpace{"obs space", this};
  oops::RequiredParameter<eckit::LocalConfiguration> obsOperator{"obs operator", this};
  oops::RequiredParameter<eckit::LocalConfiguration> geovals{"geovals", this};
  oops::Parameter<eckit::LocalConfiguration> obsBias{"obs bias", eckit::LocalConfiguration(), this};
  oops::RequiredParameter<ObsDiagnosticsParameters> obsDiagnostics{"obs diagnostics", this};
  oops::RequiredParameter<eckit::LocalConfiguration> referenceObsDiagnostics{
    "reference obs diagnostics", this};
  oops::RequiredParameter<double> tolerance{"tolerance", this};
};

// -----------------------------------------------------------------------------

void testObsDiagnostics() {
  const eckit::Configuration & conf = ::test::TestEnvironment::config();
  TestParameters params;
  params.validateAndDeserialize(conf);

  //  Setup ObsSpace
  ioda::ObsSpace ospace(params.obsSpace, oops::mpi::world(),
                        util::TimeWindow(conf.getSubConfiguration("time window")),
                        oops::mpi::myself());
  const size_t nlocs = ospace.nlocs();

  // initialize observation operator (set variables requested from the model,
  // variables simulated by the observation operator, other init)
  ObsOperator hop(ospace, params.obsOperator);

  // initialize bias correction
  const ObsBias ybias(ospace, params.obsBias);

  // initialize geovals
  oops::Variables hopvars = hop.requiredVars();
  oops::Variables reducedHopvars = ybias.requiredVars();
  hopvars += reducedHopvars;  // the reduced format is derived from the sampled format
  // read geovals from the file
  GeoVaLs gval(params.geovals, ospace, hopvars);
  // convert geovals to the reduced format
  hop.computeReducedVars(reducedHopvars, gval);

  // create obsvector to hold H(x)
  ioda::ObsVector hofx(ospace);

  // create obsvector to hold bias
  ioda::ObsVector bias(ospace);
  bias.zero();

  // create diagnostics to hold HofX diags
  const oops::ObsVariables &diagvars = params.obsDiagnostics.value().variables;
  EXPECT(diagvars.size() > 0);
  ObsDiagnostics diags(ospace, hop.locations(), diagvars);
  typedef ioda::ObsDataVector<int> QCFlags_t;
  QCFlags_t qc_flags(ospace, ospace.obsvariables(), std::string());
  // call H(x) to compute diagnostics
  hop.simulateObs(gval, hofx, ybias, qc_flags, bias, diags);

  // read tolerance and reference Diagnostics
  const double tol = params.tolerance;
  ObsDiagnostics diagref(params.referenceObsDiagnostics, ospace, diagvars);

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
