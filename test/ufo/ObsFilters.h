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
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/Variable.h"
#include "ufo/UfoTrait.h"

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  // TODO(wsmigaj) Report this (especially for floats) as a bug in eckit?
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::DateTime> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::Duration> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flags satisfy the
//! predicate in at least one variable.
//!
//! \param qcFlags
//!   Vector of quality control flags for all observations.
//! \param predicate
//!   A function object taking an argument of type int and returning bool.
//!
template <typename Predicate>
std::vector<size_t> getObservationIndicesWhere(
    const UfoTrait::ObsDataVector<int> &qcFlags, const Predicate &predicate) {
  std::vector<size_t> indices;
  for (size_t locIndex = 0; locIndex < qcFlags.nlocs(); ++locIndex) {
    bool satisfied = false;
    for (size_t varIndex = 0; varIndex < qcFlags.nvars(); ++varIndex) {
      if (predicate(qcFlags[varIndex][locIndex])) {
        satisfied = true;
        break;
      }
    }
    if (satisfied) {
      indices.push_back(locIndex);
    }
  }
  return indices;
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have passed quality control in
//! at least one variable.
//!
std::vector<size_t> getPassedObservationIndices(const UfoTrait::ObsDataVector<int> &qcFlags) {
  return getObservationIndicesWhere(qcFlags, [](int qcFlag) { return qcFlag == 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have failed quality control in
//! at least one variable.
//!
std::vector<size_t> getFailedObservationIndices(const UfoTrait::ObsDataVector<int> &qcFlags) {
  return getObservationIndicesWhere(qcFlags, [](int qcFlag) { return qcFlag != 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flag is set to \p flag in
//! at least one variable.
//!
std::vector<size_t> getFlaggedObservationIndices(const UfoTrait::ObsDataVector<int> &qcFlags,
                                                 int flag) {
  return getObservationIndicesWhere(qcFlags, [flag](int qcFlag) { return qcFlag == flag; });
}

// -----------------------------------------------------------------------------

//!
//! Return the number of elements of \p data with at least one nonzero component.
//!
size_t numNonzero(const UfoTrait::ObsDataVector<int> & data) {
  size_t result = 0;
  for (size_t locIndex = 0; locIndex < data.nlocs(); ++locIndex) {
    for (size_t varIndex = 0; varIndex < data.nvars(); ++varIndex) {
      if (data[varIndex][locIndex] != 0)
        ++result;
    }
  }
  return result;
}

// -----------------------------------------------------------------------------

//!
//! Return the number of elements of \p data with at least one component equal to \p value.
//!
size_t numEqualTo(const UfoTrait::ObsDataVector<int> & data, int value) {
  size_t result = 0;
  for (size_t locIndex = 0; locIndex < data.nlocs(); ++locIndex) {
    for (size_t varIndex = 0; varIndex < data.nvars(); ++varIndex) {
      if (data[varIndex][locIndex] == value)
        ++result;
    }
  }
  return result;
}

// -----------------------------------------------------------------------------

template <typename T>
void expectVariablesEqual(const UfoTrait::ObsSpace &obsspace,
                          const ufo::Variable &referenceVariable,
                          const ufo::Variable &testVariable)
{
  std::vector<T> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<T> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT_EQUAL(reference, test);
}

// -----------------------------------------------------------------------------

void expectVariablesApproximatelyEqual(const UfoTrait::ObsSpace &obsspace,
                                       const ufo::Variable &referenceVariable,
                                       const ufo::Variable &testVariable,
                                       float absTol)
{
  std::vector<float> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<float> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT(oops::are_all_close_absolute(reference, test, absTol));
}

// -----------------------------------------------------------------------------

void testFilters() {
  typedef ::test::ObsTestsFixture<UfoTrait> Test_;
  typedef oops::GeoVaLs<ufo::UfoTrait>           GeoVaLs_;
  typedef oops::ObsDiagnostics<ufo::UfoTrait>    ObsDiags_;
  typedef oops::ObsAuxControl<ufo::UfoTrait>     ObsAuxCtrl_;
  typedef oops::ObsFilters<ufo::UfoTrait>        ObsFilters_;
  typedef oops::ObsOperator<ufo::UfoTrait>       ObsOperator_;
  typedef oops::ObsVector<ufo::UfoTrait>         ObsVector_;
  typedef oops::ObsSpace<ufo::UfoTrait>          ObsSpace_;

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
      if (typeconfs[jj].has("ObsBias")) vars += ybias.requiredGeoVaLs();
      const eckit::LocalConfiguration gconf(typeconfs[jj], "GeoVaLs");
      const GeoVaLs_ gval(gconf, Test_::obspace()[jj], vars);
      oops::Variables diagvars;
      diagvars += filters.requiredHdiagnostics();
      if (typeconfs[jj].has("ObsBias")) diagvars += ybias.requiredHdiagnostics();
      ObsDiags_ diags(Test_::obspace()[jj],
                      hop.locations(Test_::obspace()[jj].windowStart(),
                                    Test_::obspace()[jj].windowEnd()),
                                    diagvars);
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
    bool atLeastOneBenchmarkFound = false;

    if (typeconfs[jj].has("passedObservationsBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> passedObsBenchmark =
          typeconfs[jj].getUnsignedVector("passedObservationsBenchmark");
      const std::vector<size_t> passedObs = getPassedObservationIndices(qcflags->obsdatavector());
      EXPECT_EQUAL(passedObs, passedObsBenchmark);
    }

    if (typeconfs[jj].has("passedBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const int passedBenchmark = typeconfs[jj].getInt("passedBenchmark");
      const int passed = numZero(*qcflags);
      EXPECT_EQUAL(passed, passedBenchmark);
    }

    if (typeconfs[jj].has("failedObservationsBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> failedObsBenchmark =
          typeconfs[jj].getUnsignedVector("failedObservationsBenchmark");
      const std::vector<size_t> failedObs = getFailedObservationIndices(qcflags->obsdatavector());
      EXPECT_EQUAL(failedObs, failedObsBenchmark);
    }

    if (typeconfs[jj].has("failedBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const int failedBenchmark = typeconfs[jj].getInt("failedBenchmark");
      const int failed = numNonzero(qcflags->obsdatavector());
      EXPECT_EQUAL(failed, failedBenchmark);
    }

    if (typeconfs[jj].has("benchmarkFlag")) {
      const int flag = typeconfs[jj].getInt("benchmarkFlag");

      if (typeconfs[jj].has("flaggedObservationsBenchmark")) {
        atLeastOneBenchmarkFound = true;
        const std::vector<size_t> flaggedObsBenchmark =
            typeconfs[jj].getUnsignedVector("flaggedObservationsBenchmark");
        const std::vector<size_t> flaggedObs =
            getFlaggedObservationIndices(qcflags->obsdatavector(), flag);
        EXPECT_EQUAL(flaggedObsBenchmark, flaggedObsBenchmark);
      }

      if (typeconfs[jj].has("flaggedBenchmark")) {
        atLeastOneBenchmarkFound = true;
        const int flaggedBenchmark = typeconfs[jj].getInt("flaggedBenchmark");
        const int flagged = numEqualTo(qcflags->obsdatavector(), flag);
        EXPECT_EQUAL(flagged, flaggedBenchmark);
      }
    }

    if (typeconfs[jj].has("compareVariables")) {
      for (const eckit::LocalConfiguration &compareVariablesConf :
           typeconfs[jj].getSubConfigurations("compareVariables")) {
        atLeastOneBenchmarkFound = true;

        ufo::Variable referenceVariable(compareVariablesConf.getSubConfiguration("reference"));
        ufo::Variable testVariable(compareVariablesConf.getSubConfiguration("test"));

        const UfoTrait::ObsSpace &obsspace = Test_::obspace()[jj].obsspace();
        switch (obsspace.dtype(referenceVariable.group(), referenceVariable.variable())) {
        case ioda::ObsDtype::Integer:
          expectVariablesEqual<int>(obsspace, referenceVariable, testVariable);
          break;
        case ioda::ObsDtype::String:
          expectVariablesEqual<std::string>(obsspace, referenceVariable, testVariable);
          break;
        case ioda::ObsDtype::DateTime:
          expectVariablesEqual<util::DateTime>(obsspace, referenceVariable, testVariable);
          break;
        case ioda::ObsDtype::Float:
          if (!compareVariablesConf.has("absTol")) {
            expectVariablesEqual<float>(obsspace, referenceVariable, testVariable);
          } else {
            const float tol = compareVariablesConf.getFloat("absTol");
            expectVariablesApproximatelyEqual(obsspace, referenceVariable, testVariable, tol);
          }
          break;
        case ioda::ObsDtype::None:
          ASSERT_MSG(false, "Reference variable not found in observation space");
        }
      }
    }

    EXPECT(atLeastOneBenchmarkFound);
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
