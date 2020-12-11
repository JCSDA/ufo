/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_PARALLELPOISSONDISKTHINNING_H_
#define TEST_UFO_PARALLELPOISSONDISKTHINNING_H_

#include <algorithm>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/PoissonDiskThinning.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/StringUtils.h"  // for splitVarGroup

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

void testPoissonDiskThinning(const eckit::LocalConfiguration &conf,
                             bool expectValidationError = false) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), bgn, end, oops::mpi::myself());

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "Poisson Disk Thinning");
  ufo::PoissonDiskThinning filter(obsspace, filterConf, qcflags, obserr);

  // Stagger the start of the thinning process on different ranks to detect problems related to
  // initialisation of the random generator seed from system time.
  std::this_thread::sleep_for(obsspace.comm().rank() * std::chrono::milliseconds(1100));

  filter.preProcess();

  std::vector<size_t> retainedObsIndices;
  for (size_t i = 0; i < qcflags->nlocs(); ++i)
    if ((*qcflags)[0][i] == ufo::QCflags::pass)
      retainedObsIndices.push_back(i);

  if (obsSpaceConf.getString("distribution", "RoundRobin") == "InefficientDistribution") {
    // PART 1: Verify that all ranks have retained the same observations.
    const size_t rootRank = 0;
    size_t numObsRetainedOnRoot = retainedObsIndices.size();
    obsspace.comm().broadcast(numObsRetainedOnRoot, rootRank);

    std::vector<size_t> retainedObsIndicesOnRoot(numObsRetainedOnRoot);
    if (obsspace.comm().rank() == rootRank) {
      retainedObsIndicesOnRoot = retainedObsIndices;
    }
    obsspace.comm().broadcast(retainedObsIndicesOnRoot, rootRank);

    EXPECT_EQUAL(retainedObsIndices, retainedObsIndicesOnRoot);
  }

  // PART 2: Verify that all retained observations are sufficiently far from each other
  // and that all rejected observations so close to a retained observation that they can't
  // themselves be retained.

  // Collect indices of observations retained by all processes
  size_t offset = obsspace.nlocs();
  obsspace.distribution().exclusiveScan(offset);

  std::vector<size_t> vecRetainedGlobalObsIndices = retainedObsIndices;
  for (size_t &index : vecRetainedGlobalObsIndices)
    index += offset;
  obsspace.distribution().allGatherv(vecRetainedGlobalObsIndices);
  std::set<size_t> retainedGlobalObsIndices(vecRetainedGlobalObsIndices.begin(),
                                            vecRetainedGlobalObsIndices.end());

  // Collect pressures from all processes
  std::vector<float> pressures(obsspace.nlocs());
  obsspace.get_db("MetaData", "air_pressure", pressures);
  obsspace.distribution().allGatherv(pressures);

  // Collect categories from all processes
  std::vector<int> categories(obsspace.nlocs(), 0);
  if (filterConf.has("category_variable"))
  {
    std::string name, group;
    splitVarGroup(filterConf.getString("category_variable.name"), name, group);
    obsspace.get_db(group, name, categories);
  }
  obsspace.distribution().allGatherv(categories);

  // Check distances between observations
  const float minAllowedDistance = filterConf.getFloat("min_vertical_spacing");
  std::vector<float> minDistanceToRetainedObs(pressures.size(),
                                              std::numeric_limits<float>::max());
  for (size_t i : retainedGlobalObsIndices) {
    for (size_t j = 0; j < pressures.size(); ++j) {
      if (j == i || categories[i] != categories[j])
        continue;
      const float distance = std::abs(pressures[i] - pressures[j]);
      if (oops::contains(retainedGlobalObsIndices, j)) {
        // Retained observations should be far from other retained observation
        EXPECT(distance >= minAllowedDistance);
      } else {
        minDistanceToRetainedObs[j] = std::min(minDistanceToRetainedObs[j], distance);
      }
    }
  }

  for (size_t j = 0; j < pressures.size(); ++j) {
    if (!oops::contains(retainedGlobalObsIndices, j)) {
      // Rejected observations should be close to some retained observation
      EXPECT(minDistanceToRetainedObs[j] < minAllowedDistance);
    }
  }
}

class ParallelPoissonDiskThinning : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::ParallelPoissonDiskThinning";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/ParallelPoissonDiskThinning/" + testCaseName, testCaseConf)
                      {
                        testPoissonDiskThinning(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_PARALLELPOISSONDISKTHINNING_H_
