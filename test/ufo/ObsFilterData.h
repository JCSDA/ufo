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
    const eckit::LocalConfiguration obsvarconf(obsconf, "simulate");
    ioda::ObsSpace ospace(obsconf, oops::mpi::comm(), bgn, end);

///  Setup GeoVaLs
    const eckit::LocalConfiguration gconf(confs[jconf], "GeoVaLs");
    const oops::Variables ingeovars(gconf);
    const GeoVaLs gval(gconf, ospace, ingeovars);

///  Setup ObsDiags
    const eckit::LocalConfiguration obsdiagconf(confs[jconf], "ObsDiag");
    const oops::Variables inobsdiagvars(obsdiagconf);
    const ObsDiagnostics obsdiags(obsdiagconf, ospace, inobsdiagvars);

///  Setup H(x)
    const std::string hofxgroup = confs[jconf].getString("HofX");
    ioda::ObsVector hofx(ospace, hofxgroup);

///  Setup ObsFilterData and test nlocs
    ObsFilterData data(ospace);
    EXPECT(data.nlocs() == ospace.nlocs());

///  Check that has() and get() works on ObsSpace:
    ufo::Variables obsvars(obsvarconf, "ObsValue");
    for (size_t jvar = 0; jvar < obsvars.size(); ++jvar) {
      EXPECT(data.has(obsvars[jvar]));
      std::vector<float> vec;
      data.get(obsvars[jvar], vec);
      std::vector<float> ref(ospace.nlocs());
      ospace.get_db(obsvars.group(jvar), obsvars.variable(jvar), ref.size(), ref.data());
      EXPECT(vec == ref);
    }
///  Check that has() and get() work on integer variables in ObsSpace:
    const eckit::LocalConfiguration intvarconf(obsconf, "integer data");
    ufo::Variables intvars(intvarconf);
    for (size_t jvar = 0; jvar < intvars.size(); ++jvar) {
      EXPECT(data.has(intvars[jvar]));
      std::vector<int> vec;
      data.get(intvars[jvar], vec);
      std::vector<int> ref(ospace.nlocs());
      ospace.get_db(intvars.group(jvar), intvars.variable(jvar), ref.size(), ref.data());
      EXPECT(vec == ref);
    }

///  Check that associate(), has() and get() work on ObsVector:
    ufo::Variables hofxvars(obsvarconf, "HofX");
///  H(x) not associated yet
    for (size_t jvar = 0; jvar < hofxvars.size(); ++jvar) {
      EXPECT(!data.has(hofxvars[jvar]));
    }
    data.associate(hofx);
///  H(x) associated now
    for (size_t jvar = 0; jvar < hofxvars.size(); ++jvar) {
      EXPECT(data.has(hofxvars[jvar]));
      std::vector<float> vec;
      data.get(hofxvars[jvar], vec);
      std::vector<float> ref(hofx.nlocs());
      for (size_t jloc = 0; jloc < hofx.nlocs(); jloc++) {
        ref[jloc] = hofx[hofxvars.size() * jloc + jvar];
      }
      EXPECT(vec == ref);
    }

///  Check that associate(), has() and get() work on GeoVaLs:
    ufo::Variables geovars(gconf, "GeoVaLs");
///  GeoVaLs not associated yet
    for (size_t jvar = 0; jvar < geovars.size(); ++jvar) {
      EXPECT(!data.has(geovars[jvar]));
    }
    data.associate(gval);
///  GeoVaLs associated now
    for (size_t jvar = 0; jvar < geovars.size(); ++jvar) {
      EXPECT(data.has(geovars[jvar]));
      int nlevs = data.nlevs(geovars[jvar]);
      int nlevs_ref = gval.nlevs(geovars.variable(jvar));
      EXPECT(nlevs == nlevs_ref);
///  nlevs == 1: 2D geovals, could be retrieved with get(var)
      if (nlevs == 1) {
        std::vector<float> vec;
        data.get(geovars[jvar], vec);
        std::vector<float> ref(ospace.nlocs());
        gval.get(ref, geovars.variable(jvar));
        EXPECT(vec == ref);
///  otherwise need get(var, level) to retrieve
      } else {
        std::vector<float> vec;
        data.get(geovars[jvar], nlevs, vec);
        std::vector<float> ref(ospace.nlocs());
        gval.get(ref, geovars.variable(jvar), nlevs);
        EXPECT(vec == ref);
      }
    }

///  Check that associate(), has() and get() work on ObsDiags:
    ufo::Variables diagvars(obsdiagconf, "ObsDiag");
///  ObsDiags not associated yet
    for (size_t jvar = 0; jvar < diagvars.size(); ++jvar) {
      EXPECT(!data.has(diagvars[jvar]));
    }
    data.associate(obsdiags);
///  ObsDiags associated now
    for (size_t jvar = 0; jvar < diagvars.size(); ++jvar) {
      EXPECT(data.has(diagvars[jvar]));
      int nlevs = data.nlevs(diagvars[jvar]);
      int nlevs_ref = obsdiags.nlevs(diagvars.variable(jvar));
      EXPECT(nlevs == nlevs_ref);
///  nlevs == 1: 2D obsdiags, could be retrieved with get(var)
      if (nlevs == 1) {
        std::vector<float> vec;
        data.get(diagvars[jvar], vec);
        std::vector<float> ref(ospace.nlocs());
        obsdiags.get(ref, diagvars.variable(jvar));
        EXPECT(vec == ref);
///  otherwise need get(var, level) to retrieve
      } else {
        std::vector<float> vec;
        data.get(diagvars[jvar], nlevs, vec);
        std::vector<float> ref(ospace.nlocs());
        obsdiags.get(ref, diagvars.variable(jvar), nlevs);
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
