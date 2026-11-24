/*
 * (C) Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_RTTOVEMPTYINPUT_H_
#define TEST_UFO_RTTOVEMPTYINPUT_H_

#include <string>
#include <string_view>  // NOLINT(build/include_order)
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/ObsVariables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/TimeWindow.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/QCflags.h"
#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"
#include "ufo/ObsTraits.h"

namespace ufo {
namespace test {

// Roundoff differences between batch layouts can exceed kStateTolerance for the Jacobian
// diagnostic and the TL/AD outputs derived from it, so those comparisons use kRelaxedTolerance.
constexpr double kStateTolerance = 1.0e-10;
constexpr double kRelaxedTolerance = 1.0e-8;
constexpr const char* kDiagnostic = "brightness_temperature_jacobian_air_temperature_1";

// Choose the rejected locations so that RTTOV batches are empty at the start, in the middle or at
// the end of the data, and so that accepted profiles are left with sparse profile IDs.
bool isRejected(const std::string_view pattern, const size_t loc, const size_t nlocs) {
  if (pattern == "all") return true;
  if (pattern == "first" || pattern == "qc") return loc < 2;
  if (pattern == "last") return loc >= nlocs - 2;
  if (pattern == "pairs") return (loc / 2) % 2 == 0;
  return loc % 2 == 1;  // "odd"
}

// Rejected profiles must give no output, and accepted profiles must give the same output as when
// every profile is accepted.
CASE("ufo/RTTOV/emptyInput") {
  const eckit::Configuration& conf = ::test::TestEnvironment::config();
  ioda::ObsSpace obsspace(conf.getSubConfiguration("obs space"), oops::mpi::world(),
                          util::TimeWindow(conf.getSubConfiguration("time window")),
                          oops::mpi::myself());
  const size_t nlocs = obsspace.nlocs();
  EXPECT(nlocs > 4);
  const eckit::LocalConfiguration gconf(conf, "geovals");
  const eckit::LocalConfiguration bconf;
  ObsBias bias(obsspace, bconf);
  ObsBiasIncrement biasInc(obsspace, bconf);
  biasInc.zero();
  ioda::ObsVector biasVector(obsspace);
  biasVector.zero();
  ioda::ObsDataVector<int> qc(obsspace, obsspace.obsvariables(), std::string());
  const oops::ObsVariables diagnosticVars(std::vector<std::string>{kDiagnostic});
  std::vector<float> validAzimuth;
  obsspace.get_db("MetaData", "sensorAzimuthAngle", validAzimuth);

  auto setEmissivity = [&](const std::vector<float>& emissivity) {
    for (const int channel : obsspace.assimvariables().channels()) {
      obsspace.put_db("TestEmissivity", "emissivity_" + std::to_string(channel), emissivity,
                      {"Location", "Channel"});
    }
  };
  auto compare = [&](const std::vector<bool>& rejected, const auto& actual, const auto& expected,
                     const size_t nvars, const double rejectedValue, const double tolerance) {
    for (size_t i = 0; i < nlocs * nvars; ++i) {
      if (rejected[i / nvars]) {
        EXPECT_EQUAL(actual[i], rejectedValue);
      } else {
        EXPECT(oops::is_close_relative(actual[i], expected[i], tolerance));
      }
    }
  };

  // Profiles are rejected either by an invalid sensor azimuth or by an invalid supplied
  // emissivity, and the batch size sets how many profiles RTTOV sees in each call.
  for (const bool suppliedEmissivity : {false, true}) {
    for (const int batchSize : {22, 44, 10000}) {
      eckit::LocalConfiguration opconf(conf, "obs operator");
      eckit::LocalConfiguration options(opconf, "obs options");
      options.set("max_channels_per_batch", batchSize);
      if (suppliedEmissivity) options.set("surface emissivity group", "TestEmissivity");
      opconf.set("obs options", options);
      obsspace.put_db("MetaData", "sensorAzimuthAngle", validAzimuth);
      if (suppliedEmissivity) setEmissivity(std::vector<float>(nlocs, 0.95f));
      qc.zero();

      // Reference output with every profile accepted.
      ObsOperator referenceOp(obsspace, opconf);
      GeoVaLs trajectory(gconf, obsspace, referenceOp.requiredVars());
      LinearObsOperator referenceLinearOp(obsspace, opconf);
      referenceLinearOp.setTrajectory(trajectory, bias, qc);
      GeoVaLs increment(gconf, obsspace, referenceLinearOp.requiredVars());
      increment.random();
      ioda::ObsVector referenceTL(obsspace);
      referenceLinearOp.simulateObsTL(increment, referenceTL, biasInc);
      EXPECT(referenceTL.rms() > 0.0);
      ioda::ObsVector referenceHofx(obsspace);
      ObsDiagnostics referenceDiags(obsspace, referenceOp.locations(), diagnosticVars);
      referenceOp.simulateObs(trajectory, referenceHofx, bias, qc, biasVector, referenceDiags);
      EXPECT(referenceDiags.nlevs(kDiagnostic) > 0);

      std::vector<std::string> patterns{"all", "first", "last", "pairs", "odd"};
      // The QC pattern rejects the second profile by QC flag, after an empty batch.
      if (!suppliedEmissivity) patterns.push_back("qc");
      for (const std::string& pattern : patterns) {
        oops::Log::info() << "RTTOV empty-input case: batch=" << batchSize
                          << ", supplied emissivity=" << suppliedEmissivity
                          << ", pattern=" << pattern << std::endl;
        qc.zero();
        eckit::LocalConfiguration scenarioConf(opconf);
        std::vector<bool> rejected(nlocs);
        std::vector<float> azimuth(validAzimuth);
        std::vector<float> emissivity(nlocs, 0.95f);
        for (size_t loc = 0; loc < nlocs; ++loc) {
          rejected[loc] = isRejected(pattern, loc, nlocs);
          if (!rejected[loc]) continue;
          if (pattern == "qc" && loc == 1) {
            for (size_t chan = 0; chan < qc.nvars(); ++chan) qc[chan][loc] = QCflags::black;
          } else if (suppliedEmissivity) {
            emissivity[loc] = 1.5f;
          } else {
            azimuth[loc] = 1000.0f;
          }
        }
        if (pattern == "qc") {
          eckit::LocalConfiguration scenarioOptions(options);
          scenarioOptions.set("UseQCFlagsToSkipHofX", true);
          scenarioConf.set("obs options", scenarioOptions);
        }
        obsspace.put_db("MetaData", "sensorAzimuthAngle", azimuth);
        if (suppliedEmissivity) setEmissivity(emissivity);

        ObsOperator op(obsspace, scenarioConf);
        ObsDiagnostics diags(obsspace, op.locations(), diagnosticVars);
        ioda::ObsVector hofx(obsspace);
        op.simulateObs(trajectory, hofx, bias, qc, biasVector, diags);
        compare(rejected, hofx, referenceHofx, hofx.nvars(), util::missingValue<double>(),
                kStateTolerance);
        if (pattern != "all") {
          EXPECT_EQUAL(diags.nlevs(kDiagnostic), referenceDiags.nlevs(kDiagnostic));
          for (size_t lev = 0; lev < diags.nlevs(kDiagnostic); ++lev) {
            std::vector<double> diagnostic(nlocs);
            std::vector<double> referenceDiagnostic(nlocs);
            diags.get(diagnostic, kDiagnostic, lev);
            referenceDiags.get(referenceDiagnostic, kDiagnostic, lev);
            EXPECT_EQUAL(diagnostic.size(), nlocs);
            compare(rejected, diagnostic, referenceDiagnostic, 1, util::missingValue<double>(),
                    kRelaxedTolerance);
          }
        }

        LinearObsOperator linearOp(obsspace, scenarioConf);
        linearOp.setTrajectory(trajectory, bias, qc);
        ioda::ObsVector tl(obsspace);
        linearOp.simulateObsTL(increment, tl, biasInc);
        compare(rejected, tl, referenceTL, tl.nvars(), 0.0, kRelaxedTolerance);

        // The reference AD sees no forcing for the rejected profiles.
        ioda::ObsVector forcing(obsspace);
        forcing.random();
        ioda::ObsVector maskedForcing(forcing);
        for (size_t i = 0; i < nlocs * tl.nvars(); ++i) {
          if (rejected[i / tl.nvars()]) maskedForcing[i] = 0.0;
        }
        GeoVaLs ad(increment);
        GeoVaLs referenceAD(increment);
        linearOp.simulateObsAD(ad, forcing, biasInc);
        referenceLinearOp.simulateObsAD(referenceAD, maskedForcing, biasInc);
        const double referenceRms = referenceAD.rms();
        ad -= referenceAD;
        oops::Log::info() << "AD difference rms " << ad.rms() << ", reference rms "
                          << referenceRms << std::endl;
        EXPECT(ad.rms() <= kStateTolerance + kRelaxedTolerance * referenceRms);
      }
    }
  }
}

class RTTOVEmptyInput : public oops::Test {
 private:
  std::string testid() const override { return "ufo::test::RTTOVEmptyInput"; }
  void register_tests() const override {
    // The CASE above registers itself with eckit, so there is nothing to add.
  }
  void clear() const override {
    // The test creates its own objects inside the CASE, so there is nothing to clean up.
  }
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_RTTOVEMPTYINPUT_H_
