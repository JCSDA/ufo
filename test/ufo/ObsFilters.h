/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSFILTERS_H_
#define TEST_UFO_OBSFILTERS_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

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
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsTraits.h"

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
//! \brief Convert indices of observations held by this process to global observation indices.
//!
//! It is assumed that observations are distributed to processes in a round-robin fashion.
//! For example, 8 observations are mapped to 3 processes in the following way:
//!
//! Global obs. index | Process index | Local obs. index
//! ----------------- | ------------- | ----------------
//! 0                 | 0             | 0
//! 1                 | 1             | 0
//! 2                 | 2             | 0
//! 3                 | 0             | 1
//! 4                 | 1             | 1
//! 5                 | 2             | 1
//! 6                 | 0             | 2
//! 7                 | 1             | 2
//!
void convertLocalObsIndicesToGlobal(const eckit::mpi::Comm &comm, std::vector<size_t> &indices) {
  const size_t rank = comm.rank();
  const size_t size = comm.size();
  for (size_t &index : indices)
    index = index * size + rank;
}

// -----------------------------------------------------------------------------

///
/// \brief Gather data from all tasks and deliver the combined data to all tasks.
///
/// \returns A vector that contains the elements of \p v from process 0 followed by the elements
/// of \p v from process 1 etc.
///
template <typename T>
std::vector<T> allGatherv(const eckit::mpi::Comm &comm, const std::vector<T> &v) {
  eckit::mpi::Buffer<T> buffer(comm.size());
  comm.allGatherv(v.begin(), v.end(), buffer);
  return buffer.buffer;
}

//!
//! Return the indices of observations whose quality control flags satisfy the
//! predicate in at least one variable.
//!
//! \param qcFlags
//!   Vector of quality control flags for all observations.
//! \param predicate
//!   A function object taking an argument of type int and returning bool.
//! \param comm
//!   The MPI communicator used by the ObsSpace.
//!
template <typename Predicate>
std::vector<size_t> getObservationIndicesWhere(
    const eckit::mpi::Comm &comm,
    const ObsTraits::ObsDataVector<int> &qcFlags,
    const Predicate &predicate) {
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

  convertLocalObsIndicesToGlobal(comm, indices);
  indices = allGatherv(comm, indices);
  std::sort(indices.begin(), indices.end());
  return indices;
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have passed quality control in
//! at least one variable.
//!
std::vector<size_t> getPassedObservationIndices(const eckit::mpi::Comm &comm,
                                                const ObsTraits::ObsDataVector<int> &qcFlags) {
  return getObservationIndicesWhere(comm, qcFlags, [](int qcFlag) { return qcFlag == 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have failed quality control in
//! at least one variable.
//!
std::vector<size_t> getFailedObservationIndices(const eckit::mpi::Comm &comm,
                                                const ObsTraits::ObsDataVector<int> &qcFlags) {
  return getObservationIndicesWhere(comm, qcFlags, [](int qcFlag) { return qcFlag != 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flag is set to \p flag in
//! at least one variable.
//!
std::vector<size_t> getFlaggedObservationIndices(const eckit::mpi::Comm &comm,
                                                 const ObsTraits::ObsDataVector<int> &qcFlags,
                                                 int flag) {
  return getObservationIndicesWhere(comm, qcFlags, [flag](int qcFlag) { return qcFlag == flag; });
}

// -----------------------------------------------------------------------------

//!
//! Return the number of elements of \p data with at least one nonzero component.
//!
size_t numNonzero(const ObsTraits::ObsDataVector<int> & data) {
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
size_t numEqualTo(const ObsTraits::ObsDataVector<int> & data, int value) {
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
void expectVariablesEqual(const ObsTraits::ObsSpace &obsspace,
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

void expectVariablesApproximatelyEqual(const ObsTraits::ObsSpace &obsspace,
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
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
  typedef oops::GeoVaLs<ufo::ObsTraits>           GeoVaLs_;
  typedef oops::ObsDiagnostics<ufo::ObsTraits>    ObsDiags_;
  typedef oops::ObsAuxControl<ufo::ObsTraits>     ObsAuxCtrl_;
  typedef oops::ObsFilters<ufo::ObsTraits>        ObsFilters_;
  typedef oops::ObsOperator<ufo::ObsTraits>       ObsOperator_;
  typedef oops::ObsVector<ufo::ObsTraits>         ObsVector_;
  typedef oops::ObsSpace<ufo::ObsTraits>          ObsSpace_;

  std::vector<eckit::LocalConfiguration> typeconfs;
  ::test::TestEnvironment::config().get("observations", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
/// init QC and error
    std::shared_ptr<oops::ObsDataVector<ufo::ObsTraits, float> > obserr
      (new oops::ObsDataVector<ufo::ObsTraits, float>(Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables(), "ObsError"));
    std::shared_ptr<oops::ObsDataVector<ufo::ObsTraits, int> >
      qcflags(new oops::ObsDataVector<ufo::ObsTraits, int>  (Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables()));

//  Create filters and run preProcess
    ObsFilters_ filters(Test_::obspace()[jj], typeconfs[jj], qcflags, obserr);
    filters.preProcess();

/// call priorFilter and postFilter if hofx is available
    oops::Variables geovars = filters.requiredVars();
    oops::Variables diagvars = filters.requiredHdiagnostics();
    if (typeconfs[jj].has("HofX")) {
///   read GeoVaLs from file if required
      std::unique_ptr<const GeoVaLs_> gval;
      if (geovars.size() > 0) {
        const eckit::LocalConfiguration gconf(typeconfs[jj], "geovals");
        gval.reset(new GeoVaLs_(gconf, Test_::obspace()[jj], geovars));
        filters.priorFilter(*gval);
      } else {
        oops::Log::info() << "Filters don't require geovals, priorFilter not called" << std::endl;
      }
///   read H(x) and ObsDiags from file
      oops::Log::info() << "HofX section specified, reading HofX from file" << std::endl;
      const std::string hofxgroup = typeconfs[jj].getString("HofX");
      ObsVector_ hofx(Test_::obspace()[jj], hofxgroup);
      eckit::LocalConfiguration obsdiagconf;
      if (diagvars.size() > 0) {
        obsdiagconf = eckit::LocalConfiguration(typeconfs[jj], "obs diagnostics");
        oops::Log::info() << "Obs diagnostics section specified, reading obs diagnostics from file"
                          << std::endl;
      }
      const ObsDiags_ diags(obsdiagconf, Test_::obspace()[jj], diagvars);
      filters.postFilter(hofx, diags);
    } else if (typeconfs[jj].has("obs operator")) {
///   read GeoVaLs, compute H(x) and ObsDiags
      oops::Log::info() << "ObsOperator section specified, computing HofX" << std::endl;
      const eckit::LocalConfiguration obsopconf(typeconfs[jj], "obs operator");
      ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
      const ObsAuxCtrl_ ybias(Test_::obspace()[jj], typeconfs[jj]);
      ObsVector_ hofx(Test_::obspace()[jj]);
      oops::Variables vars;
      vars += hop.requiredVars();
      vars += filters.requiredVars();
      if (typeconfs[jj].has("obs bias")) vars += ybias.requiredVars();
      const eckit::LocalConfiguration gconf(typeconfs[jj], "geovals");
      const GeoVaLs_ gval(gconf, Test_::obspace()[jj], vars);
      oops::Variables diagvars;
      diagvars += filters.requiredHdiagnostics();
      if (typeconfs[jj].has("obs bias")) diagvars += ybias.requiredHdiagnostics();
      ObsDiags_ diags(Test_::obspace()[jj],
                      hop.locations(Test_::obspace()[jj].windowStart(),
                                    Test_::obspace()[jj].windowEnd()),
                                    diagvars);
      filters.priorFilter(gval);
      hop.simulateObs(gval, hofx, ybias, diags);
      filters.postFilter(hofx, diags);
    } else if (geovars.size() > 0) {
///   Only call priorFilter
      const eckit::LocalConfiguration gconf(typeconfs[jj], "geovals");
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
    const ObsTraits::ObsSpace &obsspace = Test_::obspace()[jj].obsspace();

    if (typeconfs[jj].has("passedObservationsBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> passedObsBenchmark =
          typeconfs[jj].getUnsignedVector("passedObservationsBenchmark");
      const std::vector<size_t> passedObs = getPassedObservationIndices(
            obsspace.comm(), qcflags->obsdatavector());
      EXPECT_EQUAL(passedObs, passedObsBenchmark);
    }

    if (typeconfs[jj].has("passedBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const int passedBenchmark = typeconfs[jj].getInt("passedBenchmark");
      int passed = numEqualTo(qcflags->obsdatavector(), ufo::QCflags::pass);
      obsspace.comm().allReduceInPlace(passed, eckit::mpi::Operation::SUM);
      EXPECT_EQUAL(passed, passedBenchmark);
    }

    if (typeconfs[jj].has("failedObservationsBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> failedObsBenchmark =
          typeconfs[jj].getUnsignedVector("failedObservationsBenchmark");
      const std::vector<size_t> failedObs = getFailedObservationIndices(
            obsspace.comm(), qcflags->obsdatavector());
      EXPECT_EQUAL(failedObs, failedObsBenchmark);
    }

    if (typeconfs[jj].has("failedBenchmark")) {
      atLeastOneBenchmarkFound = true;
      const int failedBenchmark = typeconfs[jj].getInt("failedBenchmark");
      int failed = numNonzero(qcflags->obsdatavector());
      obsspace.comm().allReduceInPlace(failed, eckit::mpi::Operation::SUM);
      EXPECT_EQUAL(failed, failedBenchmark);
    }

    if (typeconfs[jj].has("benchmarkFlag")) {
      const int flag = typeconfs[jj].getInt("benchmarkFlag");

      if (typeconfs[jj].has("flaggedObservationsBenchmark")) {
        atLeastOneBenchmarkFound = true;
        const std::vector<size_t> flaggedObsBenchmark =
            typeconfs[jj].getUnsignedVector("flaggedObservationsBenchmark");
        const std::vector<size_t> flaggedObs =
            getFlaggedObservationIndices(obsspace.comm(), qcflags->obsdatavector(), flag);
        EXPECT_EQUAL(flaggedObsBenchmark, flaggedObsBenchmark);
      }

      if (typeconfs[jj].has("flaggedBenchmark")) {
        atLeastOneBenchmarkFound = true;
        const int flaggedBenchmark = typeconfs[jj].getInt("flaggedBenchmark");
        int flagged = numEqualTo(qcflags->obsdatavector(), flag);
        obsspace.comm().allReduceInPlace(flagged, eckit::mpi::Operation::SUM);
        EXPECT_EQUAL(flagged, flaggedBenchmark);
      }
    }

    if (typeconfs[jj].has("compareVariables")) {
      for (const eckit::LocalConfiguration &compareVariablesConf :
           typeconfs[jj].getSubConfigurations("compareVariables")) {
        atLeastOneBenchmarkFound = true;

        ufo::Variable referenceVariable(compareVariablesConf.getSubConfiguration("reference"));
        ufo::Variable testVariable(compareVariablesConf.getSubConfiguration("test"));

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
