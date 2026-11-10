/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSOPERATORSHARED_H_
#define TEST_UFO_OBSOPERATORSHARED_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/ObsTraits.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

/// Unlike test::ObsOperator, which builds and deletes one operator at a time, this test
/// keeps the operators of all obs spaces alive together, as an application does. Operators
/// that share resources across obs spaces (e.g. a shared CRTM coefficient load) then see
/// each other. H(x) of each operator is checked against its reference after all operators
/// are built, and again for the remaining operators after the first one is deleted.
class SharedOperators {
  typedef ::test::ObsTestsFixture<ObsTraits>  Test_;
  typedef oops::GeoVaLs<ObsTraits>            GeoVaLs_;
  typedef oops::ObsAuxControl<ObsTraits>      ObsAuxCtrl_;
  typedef oops::ObsDiagnostics<ObsTraits>     ObsDiags_;
  typedef oops::ObsOperator<ObsTraits>        ObsOperator_;
  typedef oops::ObsVector<ObsTraits>          ObsVector_;

 public:
  SharedOperators()
    : obsConfs_(::test::TestEnvironment::config().getSubConfigurations("observations")) {
    for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
      const eckit::LocalConfiguration obsOpConf(obsConfs_[jj], "obs operator");
      hops_.push_back(std::make_unique<ObsOperator_>(Test_::obspace()[jj], obsOpConf));
    }
  }

  void deleteOperator(std::size_t jj) {hops_[jj].reset();}

  /// Simulate obs space jj with its operator and compare with the reference.
  void checkSimulateObs(std::size_t jj) const {
    const eckit::Configuration & obsConf = obsConfs_[jj];
    ObsOperator_ & hop = *hops_[jj];
    oops::Log::info() << "Simulating " << Test_::obspace()[jj].obsname() << std::endl;

    oops::ObsDataVector<ObsTraits, int> qc_flags(Test_::obspace()[jj],
                                                 Test_::obspace()[jj].obsvariables(),
                                                 std::string());
    const eckit::LocalConfiguration bconf = obsConf.getSubConfiguration("obs bias");
    ObsAuxCtrl_ ybias(Test_::obspace()[jj], bconf);

    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias.requiredVars();
    hopvars += reducedHopvars;
    GeoVaLs_ gval(eckit::LocalConfiguration(obsConf, "geovals"), Test_::obspace()[jj], hopvars);
    hop.computeReducedVars(reducedHopvars, gval);

    ObsVector_ hofx(Test_::obspace()[jj]);
    ObsVector_ bias(Test_::obspace()[jj]);
    bias.zero();
    oops::ObsVariables diagvars;
    diagvars += ybias.requiredHdiagnostics();
    ObsDiags_ diags(Test_::obspace()[jj], hop.locations(), diagvars);

    hop.simulateObs(gval, hofx, ybias, qc_flags, bias, diags);

    const double tol = obsConf.getDouble("tolerance");
    if (obsConf.has("vector ref")) {
      ObsVector_ obsref(Test_::obspace()[jj], obsConf.getString("vector ref"));
      obsref -= hofx;
      oops::Log::info() << "Vector difference between reference and computed: " << obsref;
      EXPECT(obsref.rms() < 100*tol);
    } else {
      EXPECT(oops::is_close(obsConf.getDouble("rms ref"), hofx.rms(), tol));
    }
  }

  std::size_t size() const {return hops_.size();}

 private:
  std::vector<eckit::LocalConfiguration> obsConfs_;
  std::vector<std::unique_ptr<ObsOperator_>> hops_;
};

// -----------------------------------------------------------------------------

void testSimulateObsShared() {
  SharedOperators ops;
  for (std::size_t jj = 0; jj < ops.size(); ++jj) ops.checkSimulateObs(jj);
  // The remaining operators must work after another user of their resources is deleted
  ops.deleteOperator(0);
  for (std::size_t jj = 1; jj < ops.size(); ++jj) ops.checkSimulateObs(jj);
}

// -----------------------------------------------------------------------------

class ObsOperatorShared : public oops::Test {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
 public:
  using oops::Test::Test;
  virtual ~ObsOperatorShared() {}
 private:
  std::string testid() const override {return "ufo::test::ObsOperatorShared";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsOperatorShared/testSimulateObsShared")
      { testSimulateObsShared(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSOPERATORSHARED_H_
