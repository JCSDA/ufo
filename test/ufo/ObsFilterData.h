/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSFILTERDATA_H_
#define TEST_UFO_OBSFILTERDATA_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

template <typename T>
void testHasDtypeAndGet(const ufo::ObsFilterData& data, ioda::ObsSpace &ospace,
                        const ufo::Variable &var,
                        ioda::ObsDtype expectedDtype, const std::vector<T> &expectedValues,
                        const bool skipDerived = false) {
  EXPECT(data.has(var));
  EXPECT(data.dtype(var) == expectedDtype);
  // Extract variable into an std::vector
  {
    std::vector<T> vec;
    data.get(var, vec, skipDerived);
    EXPECT(vec == expectedValues);
  }

  // Extract variable into an ObsDataVector
  {
    ioda::ObsDataVector<T> vec(ospace, var.variable());
    data.get(var, vec, skipDerived);
    EXPECT_EQUAL(vec.nvars(), 1);
    EXPECT(vec[0] == expectedValues);
  }
}

// -----------------------------------------------------------------------------

template <typename T>
void testHasDtypeAndGetVector(const ufo::ObsFilterData& data, ioda::ObsSpace &ospace,
                              const ufo::Variable &var,
                              ioda::ObsDtype expectedDtype, const std::vector<T> &expectedValues,
                              const bool skipDerived = false) {
  EXPECT(data.has(var));
  EXPECT(data.dtype(var) == expectedDtype);

  // Extract variable into an std::vector
  std::vector<T> vec;
  data.get(var, vec, skipDerived);
  EXPECT(vec == expectedValues);
}

// -----------------------------------------------------------------------------

void compareMissingValues(const std::vector<float> & float_ref,
                          const std::vector<double> & double_ref) {
  EXPECT(float_ref.size() == double_ref.size());
///  Get the data in float and double types, and check the missing data values
///  match
  const float missingFloat = util::missingValue<float>();
  const float missingDouble = util::missingValue<double>();
  for (int ivar=0; ivar < float_ref.size(); ivar++) {
    if ((float_ref[ivar] == missingFloat) != (double_ref[ivar] == missingDouble)) {
      std::cout << "Missings not equal: " << ivar << "  " << float_ref[ivar] <<
        "  " << double_ref[ivar] << std::endl;
      throw std::runtime_error("Missing data does not match between precisions");
    }
  }
}

// -----------------------------------------------------------------------------
// Check the different methods to read in geovals, and check that they have
// missing data in the same locations

void checkGeoVaLsGet(const ufo::ObsFilterData& data,
                     const GeoVaLs & gval,
                     ufo::Variable variable,
                     int vectorLength,
                     int nlevs) {
  std::vector<float> vec;
  std::vector<float> ref(vectorLength);
  std::vector<double> ref2(vectorLength);
///  nlevs == 1: 2D geovals, could be retrieved with get(var)
  if (nlevs == 1) {
    data.get(variable, vec);
    gval.get(ref, variable.toOopsVariable());
    EXPECT(vec == ref);
    gval.get(ref2, variable.toOopsVariable());
    compareMissingValues(ref, ref2);
///  otherwise need get(var, level) to retrieve
  } else {
    data.get(variable, nlevs - 1, vec);
    gval.getAtLevel(ref, variable.toOopsVariable(), nlevs - 1);
    EXPECT(vec == ref);
    gval.getAtLevel(ref2, variable.toOopsVariable(), nlevs - 1);
    compareMissingValues(ref, ref2);
  }
}

// -----------------------------------------------------------------------------
// Check the different methods to read in obs diagnostics, and check that they
// contain missing data at the same locations

void checkObsDiagsGet(const ufo::ObsFilterData & data,
                      const ObsDiagnostics & obsDiags,
                      ufo::Variable variable,
                      int vectorLength,
                      int nlevs) {
  std::vector<float> vec;
  std::vector<float> ref(vectorLength);
  std::vector<double> ref2(vectorLength);
///  nlevs == 1: 2D obsdiags, could be retrieved with get(var)
  if (nlevs == 1) {
    data.get(variable, vec);
    obsDiags.get(ref, variable.variable());
    EXPECT(vec == ref);
    obsDiags.get(ref2, variable.variable());
    compareMissingValues(ref, ref2);
///  otherwise need get(var, level) to retrieve
  } else {
    data.get(variable, nlevs - 1, vec);
    obsDiags.get(ref, variable.variable(), nlevs - 1);
    EXPECT(vec == ref);
    obsDiags.get(ref2, variable.variable(), nlevs - 1);
    compareMissingValues(ref, ref2);
  }
}

// -----------------------------------------------------------------------------

void testSkipDerived(const eckit::LocalConfiguration & conf,
                     const util::TimeWindow & timeWindow,
                     const bool skipDerived) {
///  Setup ObsSpace
    const eckit::LocalConfiguration obsconf(conf, "obs space");
    ioda::ObsSpace ospace(obsconf, oops::mpi::world(), timeWindow, oops::mpi::myself());

/// Setup ObsFilterData
    ObsFilterData data(ospace);

///  Check that get() works with the skipDerived option:
    const eckit::LocalConfiguration dataconf(conf, "test data");
    std::vector<eckit::LocalConfiguration> varconfs;
    dataconf.get("float variables", varconfs);
    ufo::Variables obsvars(varconfs);
    for (size_t jvar = 0; jvar < obsvars.nvars(); ++jvar) {
      const ufo::Variable &var = obsvars.variable(jvar);
      std::vector<float> ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), ref, {}, skipDerived);
      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Float, ref, skipDerived);
    }
}

// -----------------------------------------------------------------------------

void testObsFilterData() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));
  std::vector<eckit::LocalConfiguration> confs;
  conf.get("obs filter data", confs);
  for (size_t jconf = 0; jconf < confs.size(); ++jconf) {
///  Test the skipDerived option if it is present. After doing so move to the next configuration.
    if (confs[jconf].has("skipDerived")) {
      testSkipDerived(confs[jconf], timeWindow, confs[jconf].getBool("skipDerived"));
      continue;
    }

///  Setup ObsSpace
    const eckit::LocalConfiguration obsconf(confs[jconf], "obs space");
    ioda::ObsSpace ospace(obsconf, oops::mpi::world(), timeWindow, oops::mpi::myself());

///  Setup GeoVaLs
    const eckit::LocalConfiguration gconf(confs[jconf], "geovals");
    std::vector<eckit::LocalConfiguration> varconfs;
    confs[jconf].get("geovals variables", varconfs);
    const Variables geovars(varconfs);
    GeoVaLs gval(gconf, ospace, geovars.toOopsVariables());
    gval.setDefaultFormat(GeoVaLFormat::REDUCED);

///  Setup ObsDiags
    const eckit::LocalConfiguration obsdiagconf(confs[jconf], "obs diagnostics");
    varconfs.clear();
    confs[jconf].get("obs diagnostics variables", varconfs);
    const Variables diagvars(varconfs);
    const ObsDiagnostics obsdiags(obsdiagconf, ospace, diagvars.toOopsObsVariables());

///  Setup H(x)
    const eckit::LocalConfiguration hofxconf(confs[jconf], "HofX");
    const std::string hofxgroup = hofxconf.getString("group");
    ioda::ObsVector hofx(ospace, hofxgroup);

///  Setup ObsErrors
    const eckit::LocalConfiguration obserrorconf(confs[jconf], "ObsError");
    varconfs.clear();
    obserrorconf.get("variables", varconfs);
    const Variables obserrorvars(varconfs);
    ioda::ObsDataVector<float> obserrors(ospace, obserrorvars.toOopsObsVariables());
    for (size_t i = 0; i < obserrors.nvars(); ++i)
      std::fill(obserrors[i].begin(), obserrors[i].end(), 0.1f * i);

///  Setup QCFlags
    const eckit::LocalConfiguration qcflagsconf(confs[jconf], "QCFlags");
    varconfs.clear();
    qcflagsconf.get("variables", varconfs);
    const Variables qcflagsvars(varconfs);
    ioda::ObsDataVector<int> qcflags(ospace, qcflagsvars.toOopsObsVariables());
    for (size_t i = 0; i < qcflags.nvars(); ++i)
      std::fill(qcflags[i].begin(), qcflags[i].end(), i);

///  Setup ObsFilterData and test nlocs
    ObsFilterData data(ospace);
    EXPECT(data.nlocs() == ospace.nlocs());

///  Check that has(), get() and dtype() works on ObsSpace:
    const eckit::LocalConfiguration dataconf(confs[jconf], "test data");
    varconfs.clear();
    dataconf.get("float variables", varconfs);
    ufo::Variables obsvars(varconfs);
    for (size_t jvar = 0; jvar < obsvars.nvars(); ++jvar) {
      const ufo::Variable &var = obsvars.variable(jvar);

      std::vector<float> ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), ref);

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Float, ref);
    }

    for (size_t jvar = 0; jvar < obsvars.nvars(); ++jvar) {
      const ufo::Variable &var = obsvars.variable(jvar);
      std::vector<float> float_ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), float_ref);

      std::vector<double> double_ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), double_ref);
    }

///  Check that has(), get() and dtype() work on variables in ObsSpace
///  for vectors only (e.g. metadata on Channels dimension)
    varconfs.clear();
    dataconf.get("float variables vector only", varconfs);
    ufo::Variables floatvecvars(varconfs);
    for (size_t jvar = 0; jvar < floatvecvars.nvars(); ++jvar) {
      const ufo::Variable &var = floatvecvars.variable(jvar);

      std::vector<float> ref;
      ospace.get_db(var.group(), var.variable(), ref);

      testHasDtypeAndGetVector(data, ospace, var, ioda::ObsDtype::Float, ref);
    }

///  Check that has(), get() and dtype() work on integer variables in ObsSpace:
    varconfs.clear();
    dataconf.get("integer variables", varconfs);
    ufo::Variables intvars(varconfs);
    for (size_t jvar = 0; jvar < intvars.nvars(); ++jvar) {
      const ufo::Variable &var = intvars.variable(jvar);

      std::vector<int> ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), ref);

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Integer, ref);
    }

///  Check that has(), get() and dtype() work on string variables in ObsSpace:
    varconfs.clear();
    dataconf.get("string variables", varconfs);
    ufo::Variables strvars(varconfs);
    for (size_t jvar = 0; jvar < strvars.nvars(); ++jvar) {
      const ufo::Variable &var = strvars.variable(jvar);

      std::vector<std::string> ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), ref);

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::String, ref);
    }

///  Check that has(), get() and dtype() work on datetime variables in ObsSpace:
    varconfs.clear();
    dataconf.get("datetime variables", varconfs);
    ufo::Variables dtvars(varconfs);
    for (size_t jvar = 0; jvar < dtvars.nvars(); ++jvar) {
      const ufo::Variable &var = dtvars.variable(jvar);

      std::vector<util::DateTime> ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), ref);

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::DateTime, ref);
    }

///  Check that has(), get() and dtype() work on bool variables in ObsSpace:
    varconfs.clear();
    dataconf.get("bool variables", varconfs);
    ufo::Variables boolvars(varconfs);
    for (size_t jvar = 0; jvar < boolvars.nvars(); ++jvar) {
      const ufo::Variable &var = boolvars.variable(jvar);

      std::vector<bool> ref(ospace.nlocs());
      ospace.get_db(var.group(), var.variable(), ref);

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Bool, ref);
    }

///  Check that associate(), has(), get() and dtype() work on ObsVector:
///  H(x) not associated yet
///  The important aspect of dtype handling here is it not being in ObsSpace so we
///  needn't test for other containers.
    varconfs.clear();
    hofxconf.get("variables", varconfs);
    ufo::Variables hofxvars(varconfs);
    for (size_t jvar = 0; jvar < hofxvars.nvars(); ++jvar) {
      EXPECT(!data.has(hofxvars.variable(jvar)));
    }
    data.associate(hofx, "HofX");
///  H(x) associated now
    for (size_t jvar = 0; jvar < hofxvars.nvars(); ++jvar) {
      const ufo::Variable &var = hofxvars.variable(jvar);

      std::vector<float> ref(hofx.nlocs());
      for (size_t jloc = 0; jloc < hofx.nlocs(); jloc++) {
        ref[jloc] = hofx[hofxvars.nvars() * jloc + jvar];
      }

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Float, ref);
    }

///  Check that associate(), has(), get() and dtype() work on ObsDataVector<float>:
    for (size_t jvar = 0; jvar < obserrorvars.size(); ++jvar) {
      EXPECT(!data.has(obserrorvars[jvar]));
    }
    data.associate(obserrors, "ObsErrorData");
///  ObsErrorData associated now
    for (size_t jvar = 0; jvar < obserrorvars.nvars(); ++jvar) {
      const ufo::Variable &var = obserrorvars.variable(jvar);

      const std::vector<float> &ref = obserrors[var.variable()];

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Float, ref);
    }

///  Check that associate(), has(), get() and dtype() work on ObsDataVector<int>:
    for (size_t jvar = 0; jvar < qcflagsvars.size(); ++jvar) {
      EXPECT(!data.has(qcflagsvars[jvar]));
    }
    data.associate(qcflags, "QCFlagsData");
///  QCFlags associated now
    for (size_t jvar = 0; jvar < qcflagsvars.nvars(); ++jvar) {
      const ufo::Variable &var = qcflagsvars.variable(jvar);

      const std::vector<int> &ref = qcflags[var.variable()];

      testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Integer, ref);
    }

///  Check that associate(), has() and get() work on GeoVaLs:
///  GeoVaLs not associated yet
    for (size_t jvar = 0; jvar < geovars.nvars(); ++jvar) {
      EXPECT(!data.has(geovars.variable(jvar)));
    }
    data.associate(gval);
///  GeoVaLs associated now
    for (size_t jvar = 0; jvar < geovars.nvars(); ++jvar) {
      EXPECT(data.has(geovars.variable(jvar)));
      int nlevs = data.nlevs(geovars.variable(jvar));
      int nlevs_ref = gval.nlevs(geovars.variable(jvar).toOopsVariable());
      EXPECT(nlevs == nlevs_ref);
      checkGeoVaLsGet(data, gval, geovars.variable(jvar), ospace.nlocs(), nlevs);
    }

///  Check that associate(), has() and get() work on ObsDiags:
///  ObsDiags not associated yet
    for (size_t jvar = 0; jvar < diagvars.nvars(); ++jvar) {
      EXPECT(!data.has(diagvars.variable(jvar)));
    }
    data.associate(obsdiags);
///  ObsDiags associated now
    for (size_t jvar = 0; jvar < diagvars.nvars(); ++jvar) {
      EXPECT(data.has(diagvars.variable(jvar)));
      int nlevs = data.nlevs(diagvars.variable(jvar));
      int nlevs_ref = obsdiags.nlevs(diagvars.variable(jvar).variable());
      EXPECT(nlevs == nlevs_ref);
      checkObsDiagsGet(data, obsdiags, diagvars.variable(jvar), ospace.nlocs(), nlevs);
    }

///  Check that has(), get() and dtype() work on obs functions returning floats:
    if (confs[jconf].has("float obs functions")) {
      for (const eckit::LocalConfiguration &funcconf :
           confs[jconf].getSubConfigurations("float obs functions")) {
        const eckit::LocalConfiguration varconf(funcconf, "variable");
        const ufo::Variable var(varconf);
        const float expectedValue = funcconf.getFloat("expected value");
        const std::vector<float> ref(ospace.nlocs(), expectedValue);

        testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Float, ref);
      }
    }

///  Check that has(), get() and dtype() work on obs functions returning ints:
    if (confs[jconf].has("int obs functions")) {
      for (const eckit::LocalConfiguration &funcconf :
           confs[jconf].getSubConfigurations("int obs functions")) {
        const eckit::LocalConfiguration varconf(funcconf, "variable");
        const ufo::Variable var(varconf);
        const int expectedValue = funcconf.getInt("expected value");
        const std::vector<int> ref(ospace.nlocs(), expectedValue);

        testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::Integer, ref);
      }
    }

///  Check that has(), get() and dtype() work on obs functions returning strings:
    if (confs[jconf].has("string obs functions")) {
      for (const eckit::LocalConfiguration &funcconf :
           confs[jconf].getSubConfigurations("string obs functions")) {
        const eckit::LocalConfiguration varconf(funcconf, "variable");
        const ufo::Variable var(varconf);
        const std::string expectedValue = funcconf.getString("expected value");
        const std::vector<std::string> ref(ospace.nlocs(), expectedValue);

        testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::String, ref);
      }
    }

///  Check that has(), get() and dtype() work on obs functions returning datetimes:
    if (confs[jconf].has("datetime obs functions")) {
      for (const eckit::LocalConfiguration &funcconf :
           confs[jconf].getSubConfigurations("datetime obs functions")) {
        const eckit::LocalConfiguration varconf(funcconf, "variable");
        const ufo::Variable var(varconf);
        const util::DateTime expectedValue(funcconf.getString("expected value"));
        const std::vector<util::DateTime> ref(ospace.nlocs(), expectedValue);

        testHasDtypeAndGet(data, ospace, var, ioda::ObsDtype::DateTime, ref);
      }
    }
  }
}

// -----------------------------------------------------------------------------

class ObsFilterData : public oops::Test {
 public:
  ObsFilterData() {}
  virtual ~ObsFilterData() {}
 private:
  std::string testid() const override {return "ufo::test::ObsFilterData";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFilterData/testObsFilterData")
      { testObsFilterData(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFILTERDATA_H_
