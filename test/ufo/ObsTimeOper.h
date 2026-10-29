/*
 * (C) Copyright 2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSTIMEOPER_H_
#define TEST_UFO_OBSTIMEOPER_H_

#include <algorithm>
#include <cmath>
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
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Range.h"
#include "oops/util/TimeWindow.h"
#include "test/TestEnvironment.h"
#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {
namespace test {

typedef oops::Locations<ObsTraits> Locations_;
typedef ioda::ObsDataVector<int> QCFlags_t;

// -----------------------------------------------------------------------------
/// Fill `gv` the way oops::GetValues does by default (nearest state in time): each path takes
/// the value of the model state whose time slot (state time +/- half the state spacing,
/// clipped to the window) contains the path time. State k at location i holds values[k][i].
void fillLikeGetValues(GeoVaLs & gv, const Locations_ & locs, const util::TimeWindow & window,
                       const std::vector<util::DateTime> & stateTimes,
                       const std::vector<std::vector<double>> & values) {
  const SampledLocations & paths = locs.samplingMethod(0).sampledLocations();
  const oops::Variable var = gv.getVars()[0];
  const size_t nlocs = paths.nlocs();
  std::vector<size_t> locOfPath(paths.size());
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (paths.pathsGroupedByLocation().empty()) {
      locOfPath[jloc] = jloc;
    } else {
      const util::Range<size_t> & r = paths.pathsGroupedByLocation()[jloc];
      for (size_t jp = r.begin; jp < r.end; ++jp) locOfPath[jp] = jloc;
    }
  }
  const util::Duration halfWidth = (stateTimes[1] - stateTimes[0]) / 2;
  std::vector<double> filled(paths.size(), util::missingValue<double>());
  for (size_t k = 0; k < stateTimes.size(); ++k) {
    const std::vector<bool> mask =
        window.createSubWindow(stateTimes[k], halfWidth).createTimeMask(paths.times());
    for (size_t jp = 0; jp < paths.size(); ++jp) {
      if (mask[jp]) filled[jp] = values[k][locOfPath[jp]];
    }
  }
  for (size_t jp = 0; jp < paths.size(); ++jp) {
    gv.putProfile(std::vector<double>{filled[jp]}, var, jp, GeoVaLFormat::SAMPLED);
  }
}

// -----------------------------------------------------------------------------

double dotSampled(const GeoVaLs & gv1, const GeoVaLs & gv2) {
  const oops::Variable var = gv1.getVars()[0];
  double zz = 0.0;
  std::vector<double> v1(gv1.nlevs(var, GeoVaLFormat::SAMPLED)), v2(v1.size());
  for (size_t jp = 0; jp < gv1.nprofiles(var, GeoVaLFormat::SAMPLED); ++jp) {
    gv1.getProfile(v1, var, jp, GeoVaLFormat::SAMPLED);
    gv2.getProfile(v2, var, jp, GeoVaLFormat::SAMPLED);
    for (size_t jl = 0; jl < v1.size(); ++jl) zz += v1[jl] * v2[jl];
  }
  return zz;
}

// -----------------------------------------------------------------------------

struct TimeOperTestSetup {
  TimeOperTestSetup()
    : conf(::test::TestEnvironment::config()),
      window(conf.getSubConfiguration("time window")),
      ospace(conf.getSubConfiguration("obs space"), oops::mpi::world(), window,
             oops::mpi::myself()) {
    for (const std::string & t : conf.getStringVector("state times"))
      stateTimes.push_back(util::DateTime(t));
    for (const eckit::LocalConfiguration & c : conf.getSubConfigurations("state values"))
      values.push_back(c.getDoubleVector("values"));
    reference = conf.getDoubleVector("reference hofx");
    tolerance = conf.getDouble("tolerance");
  }
  const eckit::LocalConfiguration conf;
  util::TimeWindow window;
  ioda::ObsSpace ospace;
  std::vector<util::DateTime> stateTimes;
  std::vector<std::vector<double>> values;
  std::vector<double> reference;
  double tolerance;
};

// -----------------------------------------------------------------------------
/// Each observation must be sampled at the two state times bracketing it (nudged inside the
/// window where a bound is excluded).
void testLocations() {
  TimeOperTestSetup s;
  ObsOperator hop(s.ospace, s.conf.getSubConfiguration("obs operator"));
  const Locations_ locs = hop.locations();
  const SampledLocations & paths = locs.samplingMethod(0).sampledLocations();
  const std::vector<std::string> ref = s.conf.getStringVector("reference path times");
  EXPECT_EQUAL(paths.nlocs(), s.ospace.nlocs());
  EXPECT_EQUAL(paths.size(), ref.size());
  for (size_t jp = 0; jp < std::min(paths.size(), ref.size()); ++jp)
    EXPECT_EQUAL(paths.times()[jp], util::DateTime(ref[jp]));
}

// -----------------------------------------------------------------------------
/// H(x) must be the linear time interpolation of the bracketing states, not either state.
void testSimulateObs() {
  TimeOperTestSetup s;
  ObsOperator hop(s.ospace, s.conf.getSubConfiguration("obs operator"));
  ObsBias ybias(s.ospace, eckit::LocalConfiguration());
  const Locations_ locs = hop.locations();
  GeoVaLs gv(locs, hop.requiredVars(), std::vector<size_t>{1});
  fillLikeGetValues(gv, locs, s.window, s.stateTimes, s.values);
  hop.computeReducedVars(ybias.requiredVars(), gv);

  ioda::ObsVector hofx(s.ospace), bias(s.ospace);
  bias.zero();
  ObsDiagnostics diags(s.ospace, hop.locations(), oops::ObsVariables());
  QCFlags_t qc(s.ospace, s.ospace.obsvariables(), std::string());
  hop.simulateObs(gv, hofx, ybias, qc, bias, diags);

  EXPECT_EQUAL(hofx.nlocs(), s.reference.size());
  for (size_t jo = 0; jo < hofx.nlocs(); ++jo) {
    oops::Log::info() << "obs " << jo << " hofx = " << hofx[jo]
                      << " reference = " << s.reference[jo] << std::endl;
    EXPECT(std::abs(hofx[jo] - s.reference[jo]) < s.tolerance);
  }
}

// -----------------------------------------------------------------------------
/// The TL must interpolate like the forward operator and the AD must be its adjoint.
void testTangentLinearAndAdjoint() {
  TimeOperTestSetup s;
  ObsOperator hop(s.ospace, s.conf.getSubConfiguration("obs operator"));
  LinearObsOperator hoptl(s.ospace, s.conf.getSubConfiguration("obs operator"));
  ObsBias ybias(s.ospace, eckit::LocalConfiguration());
  ObsBiasIncrement ybinc(s.ospace, eckit::LocalConfiguration());
  QCFlags_t qc(s.ospace, s.ospace.obsvariables(), std::string());
  const Locations_ locs = hop.locations();

  GeoVaLs traj(locs, hop.requiredVars(), std::vector<size_t>{1});
  fillLikeGetValues(traj, locs, s.window, s.stateTimes, s.values);
  hoptl.setTrajectory(traj, ybias, qc);

  // TL: the operator is linear, so H'(dx) with dx = x reproduces the reference
  GeoVaLs dx(locs, hop.requiredVars(), std::vector<size_t>{1});
  fillLikeGetValues(dx, locs, s.window, s.stateTimes, s.values);
  ioda::ObsVector dy(s.ospace);
  hoptl.simulateObsTL(dx, dy, ybinc);
  for (size_t jo = 0; jo < dy.nlocs(); ++jo)
    EXPECT(std::abs(dy[jo] - s.reference[jo]) < s.tolerance);

  // AD: <H' dx, dy2> == <dx, H'^T dy2>
  ioda::ObsVector dy2(s.ospace);
  for (size_t jo = 0; jo < dy2.nlocs(); ++jo) dy2[jo] = 1.0 + jo;
  GeoVaLs dxad(locs, hop.requiredVars(), std::vector<size_t>{1});
  dxad.zero();
  hoptl.simulateObsAD(dxad, dy2, ybinc);
  const double zy = dy.dot_product_with(dy2);
  const double zx = dotSampled(dx, dxad);
  oops::Log::info() << "<H'dx, dy> = " << zy << " <dx, H'^T dy> = " << zx << std::endl;
  EXPECT(std::abs(zy - zx) < 1.0e-10 * std::abs(zy));
}

// -----------------------------------------------------------------------------

class ObsTimeOper : public oops::Test {
 public:
  ObsTimeOper() = default;
  virtual ~ObsTimeOper() = default;

 private:
  std::string testid() const override {return "ufo::test::ObsTimeOper";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsTimeOper/testLocations")
                    { testLocations(); });
    ts.emplace_back(CASE("ufo/ObsTimeOper/testSimulateObs")
                    { testSimulateObs(); });
    ts.emplace_back(CASE("ufo/ObsTimeOper/testTangentLinearAndAdjoint")
                    { testTangentLinearAndAdjoint(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSTIMEOPER_H_
