/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSFILTERS_H_
#define TEST_UFO_OBSFILTERS_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsFilters.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/UfoTrait.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have passed quality control in
//! at least one variable.
//!
std::vector<size_t> getPassedObservationIndices(const UfoTrait::ObsDataVector<int> &qcFlags) {
  std::vector<size_t> indices;
  for (size_t locIndex = 0; locIndex < qcFlags.nlocs(); ++locIndex) {
    bool passed = false;
    for (size_t varIndex = 0; varIndex < qcFlags.nvars(); ++varIndex) {
      if (qcFlags[varIndex][locIndex] == 0) {
        passed = true;
        break;
      }
    }
    if (passed) {
      indices.push_back(locIndex);
    }
  }
  return indices;
}

void testFilters() {
  typedef ::test::ObsTestsFixture<UfoTrait> Test_;
  typedef oops::GeoVaLs<ufo::UfoTrait>           GeoVaLs_;
  typedef oops::ObsDiagnostics<ufo::UfoTrait>    ObsDiags_;
  typedef oops::ObsAuxControl<ufo::UfoTrait>     ObsAuxCtrl_;
  typedef oops::ObsFilters<ufo::UfoTrait>        ObsFilters_;
  typedef oops::ObsOperator<ufo::UfoTrait>       ObsOperator_;
  typedef oops::ObsVector<ufo::UfoTrait>         ObsVector_;

  const eckit::LocalConfiguration obsconf(::test::TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> typeconfs;
  obsconf.get("ObsTypes", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
/// init QC and error
    boost::shared_ptr<oops::ObsDataVector<ufo::UfoTrait, float> > obserr
      (new oops::ObsDataVector<ufo::UfoTrait, float>(Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables(), "ObsError"));
    boost::shared_ptr<oops::ObsDataVector<ufo::UfoTrait, int> >
      qcflags(new oops::ObsDataVector<ufo::UfoTrait, int>  (Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables()));

//  Create filters and run preProcess
    ObsFilters_ filters(Test_::obspace()[jj], typeconfs[jj], qcflags, obserr);
    filters.preProcess();

/// call priorFilter and postFilter if hofx is available
    oops::Variables geovars = filters.requiredGeoVaLs();
    oops::Variables diagvars = filters.requiredHdiagnostics();
    if (typeconfs[jj].has("HofX")) {
///   read GeoVaLs from file if required
      if (geovars.size() > 0) {
        const eckit::LocalConfiguration gconf(typeconfs[jj], "GeoVaLs");
        const GeoVaLs_ gval(gconf, Test_::obspace()[jj], geovars);
        filters.priorFilter(gval);
      } else {
        oops::Log::info() << "Filters don't require geovals, priorFilter not called" << std::endl;
      }
///   read H(x) and ObsDiags from file
      oops::Log::info() << "HofX section specified, reading HofX from file" << std::endl;
      const std::string hofxgroup = typeconfs[jj].getString("HofX");
      ObsVector_ hofx(Test_::obspace()[jj], hofxgroup);
      eckit::LocalConfiguration obsdiagconf;
      if (diagvars.size() > 0) {
        obsdiagconf = eckit::LocalConfiguration(typeconfs[jj], "ObsDiag");
        oops::Log::info() << "ObsDiag section speciifed, reading ObsDiag from file" << std::endl;
      }
      const ObsDiags_ diags(obsdiagconf, Test_::obspace()[jj], diagvars);
      filters.postFilter(hofx, diags);
    } else if (typeconfs[jj].has("ObsOperator")) {
///   read GeoVaLs, compute H(x) and ObsDiags
      oops::Log::info() << "ObsOperator section specified, computing HofX" << std::endl;
      const eckit::LocalConfiguration obsopconf(typeconfs[jj], "ObsOperator");
      ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
      const ObsAuxCtrl_ ybias(Test_::obspace()[jj], typeconfs[jj]);
      ObsVector_ hofx(Test_::obspace()[jj]);
      oops::Variables vars;
      vars += hop.variables();
      vars += filters.requiredGeoVaLs();
      const eckit::LocalConfiguration gconf(typeconfs[jj], "GeoVaLs");
      const GeoVaLs_ gval(gconf, Test_::obspace()[jj], vars);
      ObsDiags_ diags(Test_::obspace()[jj],
                      hop.locations(Test_::obspace()[jj].windowStart(),
                                    Test_::obspace()[jj].windowEnd()),
                      filters.requiredHdiagnostics());
      filters.priorFilter(gval);
      hop.simulateObs(gval, hofx, ybias, diags);
      filters.postFilter(hofx, diags);
    } else if (geovars.size() > 0) {
///   Only call priorFilter
      const eckit::LocalConfiguration gconf(typeconfs[jj], "GeoVaLs");
      const GeoVaLs_ gval(gconf, Test_::obspace()[jj], geovars);
      filters.priorFilter(gval);
      oops::Log::info() << "HofX or ObsOperator sections not provided for filters, " <<
                           "postFilter not called" << std::endl;
    } else {
///   no need to run priorFilter or postFilter
      oops::Log::info() << "GeoVaLs not required, HofX or ObsOperator sections not " <<
                           "provided for filters, only preProcess was called" << std::endl;
    }

    qcflags->save("EffectiveQC");
    const std::string errname = "EffectiveError";
    obserr->save(errname);

//  Compare with known results
    if (typeconfs[jj].has("passedObservationsBenchmark")) {
      const std::vector<size_t> passedObsBenchmark =
          typeconfs[jj].getUnsignedVector("passedObservationsBenchmark");
      const std::vector<size_t> passedObs = getPassedObservationIndices(qcflags->obsdatavector());
      EXPECT_EQUAL(passedObs, passedObsBenchmark);
    }

    const int passedBenchmark = typeconfs[jj].getInt("passedBenchmark");
    const int passed = numZero(*qcflags);
    EXPECT_EQUAL(passed, passedBenchmark);
  }
}

// -----------------------------------------------------------------------------

class ObsFilters : public oops::Test {
 public:
  ObsFilters() {}
  virtual ~ObsFilters() {}
 private:
  std::string testid() const {return "test::ObsFilters";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFilters/testFilters")
      { testFilters(); });
  }
};

// =============================================================================

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFILTERS_H_
