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
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

void testObsFilterData() {
  const eckit::LocalConfiguration conf = ::test::TestEnvironment::config();
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));

  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ObsFilterData", confs);
  for (size_t jconf = 0; jconf < confs.size(); ++jconf) {
///  Setup ObsSpace
    const eckit::LocalConfiguration obsconf(confs[jconf], "ObsSpace");
    ioda::ObsSpace ospace(obsconf, oops::mpi::comm(), bgn, end);

///  Setup GeoVaLs
    const eckit::LocalConfiguration gconf(confs[jconf], "GeoVaLs");
    std::vector<eckit::LocalConfiguration> varconfs;
    gconf.get("variables", varconfs);
    const Variables geovars(varconfs);
    const GeoVaLs gval(gconf, ospace, geovars.toOopsVariables());

///  Setup ObsDiags
    const eckit::LocalConfiguration obsdiagconf(confs[jconf], "ObsDiag");
    obsdiagconf.get("variables", varconfs);
    const Variables diagvars(varconfs);
    const ObsDiagnostics obsdiags(obsdiagconf, ospace, diagvars.toOopsVariables());

///  Setup H(x)
    const eckit::LocalConfiguration hofxconf(confs[jconf], "HofX");
    const std::string hofxgroup = hofxconf.getString("group");
    ioda::ObsVector hofx(ospace, hofxgroup);

///  Setup ObsFilterData and test nlocs
    ObsFilterData data(ospace);
    EXPECT(data.nlocs() == ospace.nlocs());

///  Check that has() and get() works on ObsSpace:
    obsconf.get("float variables", varconfs);
    ufo::Variables obsvars(varconfs);
    for (size_t jvar = 0; jvar < obsvars.nvars(); ++jvar) {
      EXPECT(data.has(obsvars.variable(jvar)));
      std::vector<float> vec;
      data.get(obsvars.variable(jvar), vec);
      std::vector<float> ref(ospace.nlocs());
      ospace.get_db(obsvars.variable(jvar).group(), obsvars.variable(jvar).variable(), ref);
      EXPECT(vec == ref);
    }
///  Check that has() and get() work on integer variables in ObsSpace:
    obsconf.get("integer variables", varconfs);
    ufo::Variables intvars(varconfs);
    for (size_t jvar = 0; jvar < intvars.nvars(); ++jvar) {
      EXPECT(data.has(intvars.variable(jvar)));
      std::vector<int> vec;
      data.get(intvars.variable(jvar), vec);
      std::vector<int> ref(ospace.nlocs());
      ospace.get_db(intvars.variable(jvar).group(), intvars.variable(jvar).variable(), ref);
      EXPECT(vec == ref);
    }

///  Check that associate(), has() and get() work on ObsVector:
///  H(x) not associated yet
    hofxconf.get("variables", varconfs);
    ufo::Variables hofxvars(varconfs);
    for (size_t jvar = 0; jvar < hofxvars.nvars(); ++jvar) {
      EXPECT(!data.has(hofxvars.variable(jvar)));
    }
    data.associate(hofx);
///  H(x) associated now
    for (size_t jvar = 0; jvar < hofxvars.nvars(); ++jvar) {
      EXPECT(data.has(hofxvars.variable(jvar)));
      std::vector<float> vec;
      data.get(hofxvars.variable(jvar), vec);
      std::vector<float> ref(hofx.nlocs());
      for (size_t jloc = 0; jloc < hofx.nlocs(); jloc++) {
        ref[jloc] = hofx[hofxvars.nvars() * jloc + jvar];
      }
      EXPECT(vec == ref);
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
      int nlevs_ref = gval.nlevs(geovars.variable(jvar).variable());
      EXPECT(nlevs == nlevs_ref);
///  nlevs == 1: 2D geovals, could be retrieved with get(var)
      if (nlevs == 1) {
        std::vector<float> vec;
        data.get(geovars.variable(jvar), vec);
        std::vector<float> ref(ospace.nlocs());
        gval.get(ref, geovars.variable(jvar).variable());
        EXPECT(vec == ref);
///  otherwise need get(var, level) to retrieve
      } else {
        std::vector<float> vec;
        data.get(geovars.variable(jvar), nlevs, vec);
        std::vector<float> ref(ospace.nlocs());
        gval.get(ref, geovars.variable(jvar).variable(), nlevs);
        EXPECT(vec == ref);
      }
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
///  nlevs == 1: 2D obsdiags, could be retrieved with get(var)
      if (nlevs == 1) {
        std::vector<float> vec;
        data.get(diagvars.variable(jvar), vec);
        std::vector<float> ref(ospace.nlocs());
        obsdiags.get(ref, diagvars.variable(jvar).variable());
        EXPECT(vec == ref);
///  otherwise need get(var, level) to retrieve
      } else {
        std::vector<float> vec;
        data.get(diagvars.variable(jvar), nlevs, vec);
        std::vector<float> ref(ospace.nlocs());
        obsdiags.get(ref, diagvars.variable(jvar).variable(), nlevs);
        EXPECT(vec == ref);
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
  std::string testid() const {return "ufo::test::ObsFilterData";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFilterData/testObsFilterData")
      { testObsFilterData(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFILTERDATA_H_
