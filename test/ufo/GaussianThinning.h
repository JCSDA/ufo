/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_GAUSSIANTHINNING_H_
#define TEST_UFO_GAUSSIANTHINNING_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_shared.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "ufo/filters/Gaussian_Thinning.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testGaussianThinning(const eckit::LocalConfiguration &conf) {
  util::DateTime bgn(conf.getString("window_begin"));
  util::DateTime end(conf.getString("window_end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "ObsSpace");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::comm(), bgn, end);

  if (conf.has("air_pressures")) {
    const std::vector<float> air_pressures = conf.getFloatVector("air_pressures");
    obsspace.put_db("MetaData", "air_pressure", air_pressures);
    const std::vector<float> air_pressure_obserrors(air_pressures.size(), 1.0f);
    obsspace.put_db("ObsError", "air_pressure", air_pressure_obserrors);
  }

  if (conf.has("category")) {
    const std::vector<int> categories = conf.getIntVector("category");
    obsspace.put_db("MetaData", "category", categories);
  }

  if (conf.has("priority")) {
    const std::vector<int> priorities = conf.getIntVector("priority");
    obsspace.put_db("MetaData", "priority", priorities);
  }

  auto obserr = boost::make_shared<ioda::ObsDataVector<float>>(
      obsspace, obsspace.obsvariables(), "ObsError");
  auto qcflags = boost::make_shared<ioda::ObsDataVector<int>>(
      obsspace, obsspace.obsvariables());

  const eckit::LocalConfiguration filterConf(conf, "GaussianThinning");
  ufo::Gaussian_Thinning filter(obsspace, filterConf, qcflags, obserr);
  filter.preProcess();

  const std::vector<size_t> expectedThinnedObsIndices =
      conf.getUnsignedVector("expected_thinned_obs_indices");
  std::vector<size_t> thinnedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::thinned)
      thinnedObsIndices.push_back(i);
  EXPECT_EQUAL(thinnedObsIndices.size(), expectedThinnedObsIndices.size());
  const bool equal = std::equal(thinnedObsIndices.begin(), thinnedObsIndices.end(),
                                expectedThinnedObsIndices.begin());
  EXPECT(equal);
}

CASE("ufo/GaussianThinning/Horizontal mesh 20000") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal mesh 20000"));
}

CASE("ufo/GaussianThinning/Horizontal mesh 20000, extreme longitudes") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal mesh 20000, extreme longitudes"));
}

CASE("ufo/GaussianThinning/Horizontal mesh 20000, extreme latitudes") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Horizontal mesh 20000, extreme latitudes"));
}

CASE("ufo/GaussianThinning/Vertical mesh, single bin") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, single bin"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins, all observations in single bin") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins, "
                                                 "all observations in single bin"));
}

CASE("ufo/GaussianThinning/Thinning in time, single bin") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Thinning in time, single bin"));
}

CASE("ufo/GaussianThinning/Thinning in time, two bins") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Thinning in time, two bins"));
}

CASE("ufo/GaussianThinning/Thinning in time, two bins, all observations in single bin") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Thinning in time, two bins, "
                                                 "all observations in single bin"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins, single category") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins, single category"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins, two categories") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins, two categories"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins, two categories, where clause") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins, "
                                                 "two categories, where clause"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins, equal priorities") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins, equal priorities"));
}

CASE("ufo/GaussianThinning/Vertical mesh, two bins, nonequal priorities") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, two bins, nonequal priorities"));
}

CASE("ufo/GaussianThinning/Vertical mesh, single bin, nonequal priorities, where clause") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, single bin, nonequal priorities, "
                                                 "where clause"));
}

CASE("ufo/GaussianThinning/Vertical mesh, single bin, two categories, nonequal priorities, "
     "where clause") {
  testGaussianThinning(eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                                 "Vertical mesh, single bin, two categories, "
                                                 "nonequal priorities, where clause"));
}

class GaussianThinning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::GaussianThinning";}

  void register_tests() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_GAUSSIANTHINNING_H_
