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
#include "ioda/distribution/DistributionUtils.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/missingValues.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  template <> struct VectorPrintSelector<util::DateTime> { typedef VectorPrintSimple selector; };
}  // namespace eckit

// -----------------------------------------------------------------------------

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

const char *expectComputeToThrow = "expect compute to throw exception with message";
const char *expectConstructorToThrow = "expect constructor to throw exception with message";

// -----------------------------------------------------------------------------

/// \param[out] num_missing_mismatches
///   Number of locations containing a missing value in `vals` but not in `ref`, or in `ref` but not
///   in `vals`.
void dataVectorDiff(const ioda::ObsSpace & ospace, ioda::ObsDataVector<float> & vals,
                    const ioda::ObsDataVector<float> & ref, std::vector<float> & rms_out,
                    size_t &num_missing_mismatches) {
  const float missing = util::missingValue<float>();
  num_missing_mismatches = 0;
  /// Loop through variables and calculate rms for each variable
  for (size_t ivar = 0; ivar < vals.nvars() ; ++ivar) {
    for (size_t jj = 0; jj < ref.nlocs() ; ++jj) {
      if (vals[ivar][jj] != missing && ref[ivar][jj] != missing) {
        vals[ivar][jj] -= ref[ivar][jj];
      } else {
        vals[ivar][jj] = missing;
      }
      if ((vals[ivar][jj] != missing) ^ (ref[ivar][jj] != missing)) {
        ++num_missing_mismatches;
      }
    }
    int nobs = globalNumNonMissingObs(*ospace.distribution(), 1, vals[ivar]);
    float rms = dotProduct(*ospace.distribution(), 1, vals[ivar], vals[ivar]);
    if (nobs > 0) rms = sqrt(rms / static_cast<float>(nobs));
    rms_out[ivar] = rms;
  }
}

// -----------------------------------------------------------------------------

/// Checks that `vals` and `vals_ofd` are equal to `ref`.
///
/// \param vals
///   ObsFunction values calculated by calling ObsFunction::compute() directly.
/// \param vals_ofd
///   ObsFunction values calculated by calling ObsFilterData::get().
/// \param ref
///   Reference values.
template <typename T>
void checkResults(const ioda::ObsSpace &/*ospace*/,
                  const eckit::Configuration &/*obsfuncconf*/,
                  const ioda::ObsDataVector<T> & vals,
                  const ioda::ObsDataVector<T> & vals_ofd,
                  const ioda::ObsDataVector<T> & ref) {
  const size_t nvars = ref.nvars();
  EXPECT_EQUAL(vals.nvars(), nvars);

  for (size_t i = 0; i < nvars; ++i)
    EXPECT_EQUAL(vals[i], ref[i]);

  for (size_t i = 0; i < nvars; ++i)
    EXPECT_EQUAL(vals_ofd[i], ref[i]);
}

/// Specialization for floats: expects the RMSs of the differences between `vals` and `ref` and
/// between `vals_ofd` and `ref` to be within a specified tolerance.
///
/// The contents of `vals` and `vals_ofd` are destroyed.
void checkResults(const ioda::ObsSpace &ospace,
                  const eckit::Configuration &obsfuncconf,
                  ioda::ObsDataVector<float> & vals,
                  ioda::ObsDataVector<float> & vals_ofd,
                  const ioda::ObsDataVector<float> & ref) {
  const double tol = obsfuncconf.getDouble("tolerance");
  const bool expectMissingToMatch =
      obsfuncconf.getBool("expect missing value locations to match", false);
  const size_t nvars = ref.nvars();
  EXPECT_EQUAL(vals.nvars(), nvars);

  ///  Calculate rms(f(x) - ref) and compare to tolerance
  std::vector<float> rms_out(nvars);
  size_t numMissingMismatches;
  dataVectorDiff(ospace, vals, ref, rms_out, numMissingMismatches);

  oops::Log::info() << "Vector difference between reference and computed: " << std::endl;
  oops::Log::info() << vals << std::endl;
  for (size_t ivar = 0; ivar < nvars; ivar++) {
    // FIXME(someone): whatever this does, it certainly doesn't convert percentages to fractions
    EXPECT(rms_out[ivar] < 100*tol);  //  change tol from percent to actual value.
  }
  if (expectMissingToMatch)
    EXPECT_EQUAL(numMissingMismatches, 0);

  dataVectorDiff(ospace, vals_ofd, ref, rms_out, numMissingMismatches);
  oops::Log::info() << "Vector difference between reference and computed via ObsFilterData: "
                        << std::endl;
  oops::Log::info() << vals_ofd << std::endl;
  for (size_t ivar = 0; ivar < nvars; ivar++) {
    // FIXME(someone): whatever this does, it certainly doesn't convert percentages to fractions
    EXPECT(rms_out[ivar] < 100*tol);  //  change tol from percent to actual value.
  }
  if (expectMissingToMatch)
    EXPECT_EQUAL(numMissingMismatches, 0);
}

// -----------------------------------------------------------------------------

template <typename T>
void doTestFunction(ioda::ObsSpace &ospace, const eckit::Configuration &conf) {
///  Get function name
  const eckit::LocalConfiguration obsfuncconf(conf, "obs function");
  Variable funcname(obsfuncconf);

///  Setup function
  if (conf.has(expectConstructorToThrow)) {
    // The constructor is expected to throw an exception containing the specified string.
    const std::string expectedMessage = conf.getString(expectConstructorToThrow);
    EXPECT_THROWS_MSG(ObsFunction<T>{funcname}, expectedMessage.c_str());
    return;
  }
  ObsFunction<T> obsfunc(funcname);
  ufo::Variables allfuncvars = obsfunc.requiredVariables();

///  Setup ObsFilterData
  ObsFilterData inputs(ospace);

///  Setup GeoVaLs
  const oops::Variables geovars = allfuncvars.allFromGroup("GeoVaLs").toOopsVariables();
  std::unique_ptr<GeoVaLs> gval;
  if (geovars.size() > 0) {
    const eckit::LocalConfiguration gconf(conf, "geovals");
    gval.reset(new GeoVaLs(gconf, ospace, geovars));
    gval->setDefaultFormat(GeoVaLFormat::REDUCED);
    inputs.associate(*gval);
  }

///  Setup zero ObsBias
  ioda::ObsVector bias(ospace);
  bias.zero();
  inputs.associate(bias, "ObsBiasData");

///  Setup ObsDiags
  const oops::ObsVariables diagvars = allfuncvars.allFromGroup("ObsDiag").toOopsObsVariables();
  std::unique_ptr<ObsDiagnostics> diags;
  if (diagvars.size() > 0) {
    const eckit::LocalConfiguration diagconf(conf, "obs diagnostics");
    diags.reset(new ObsDiagnostics(diagconf, ospace, diagvars));
    inputs.associate(*diags);
  }

///  Get output variable names
  const oops::ObsVariables outputvars(obsfuncconf, "variables");
///  Compute function result
  ioda::ObsDataVector<T> vals(ospace, outputvars, funcname.group(), false);
  if (!conf.has(expectComputeToThrow)) {
    obsfunc.compute(inputs, vals);
  } else {
    // The call to compute() is expected to throw an exception containing the specified string.
    const std::string expectedMessage = conf.getString(expectComputeToThrow);
    EXPECT_THROWS_MSG(obsfunc.compute(inputs, vals), expectedMessage.c_str());
    return;
  }
  vals.save("TestResult");

///  Compute function result through ObsFilterData
  ioda::ObsDataVector<T> vals_ofd(ospace, outputvars, funcname.group(), false);
  inputs.get(funcname, vals_ofd);

///  Read reference values from ObsSpace
  ioda::ObsDataVector<T> ref(ospace, outputvars, "TestReference");

  checkResults(ospace, obsfuncconf, vals, vals_ofd, ref);
}

// -----------------------------------------------------------------------------

void testFunction() {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;

  std::vector<eckit::LocalConfiguration> typeconfs;
  ::test::TestEnvironment::config().get("observations", typeconfs);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ioda::ObsSpace &ospace = Test_::obspace()[jj].obsspace();
    const eckit::Configuration &conf = typeconfs[jj];

/// Use the function group to determine the type of values produced by the function and thus
/// select the right specialization of doTestFunction().
    Variable funcname(conf.getSubConfiguration("obs function"));
    if (funcname.group() == ObsFunctionTraits<float>::groupName)
      doTestFunction<float>(ospace, conf);
    else if (funcname.group() == ObsFunctionTraits<int>::groupName)
      doTestFunction<int>(ospace, conf);
    else if (funcname.group() == ObsFunctionTraits<std::string>::groupName)
      doTestFunction<std::string>(ospace, conf);
    else if (funcname.group() == ObsFunctionTraits<util::DateTime>::groupName)
      doTestFunction<util::DateTime>(ospace, conf);
    else
      throw eckit::testing::TestException("Variable " + funcname.fullName() +
                                          " is not an ObsFunction", Here());
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

  void clear() const override {
    ::test::ObsTestsFixture<ObsTraits>::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFUNCTION_H_
