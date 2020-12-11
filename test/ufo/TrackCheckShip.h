/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_TRACKCHECKSHIP_H_
#define TEST_UFO_TRACKCHECKSHIP_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <boost/none.hpp>
#include <boost/optional.hpp>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/TrackCheckShip.h"
#include "ufo/filters/TrackCheckShipDiagnostics.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

const boost::optional<ufo::TrackCheckShipDiagnostics> setupRunFilter(
    const eckit::LocalConfiguration &conf, std::vector<size_t> *rejectedObsIndices = nullptr) {
  util::DateTime bgn(conf.getString("window begin"));
  util::DateTime end(conf.getString("window end"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), bgn, end, oops::mpi::myself());

  if (conf.has("station_ids")) {
    const std::vector<int> stationIds = conf.getIntVector("station_ids");
    obsspace.put_db("MetaData", "station_id", stationIds);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  const eckit::LocalConfiguration filterConf(conf, "Ship Track Check");
  ufo::TrackCheckShip filter(obsspace, filterConf, qcflags, obserr);
  filter.preProcess();
  if (filterConf.getBool("comparison test", false) && rejectedObsIndices) {
    for (size_t i = 0; i < qcflags->nlocs(); ++i)
      if ((*qcflags)[0][i] == ufo::QCflags::track)
        rejectedObsIndices->push_back(i);
    return boost::none;
  } else if (filterConf.getBool("unit testing mode", false)) {
    return *filter.diagnostics();
  }
  return boost::none;
}

void testTrackCheckShipInitialCalculations(const eckit::LocalConfiguration &conf) {
  auto diagnostics = setupRunFilter(conf);
  if (!diagnostics)
    return;
  const std::vector<double> expectedDistances = conf.getDoubleVector("expected distance");
  const std::vector<double> expectedSpeeds = conf.getDoubleVector("expected speed");
  const std::vector<double> expectedDistancesAveraged =
      conf.getDoubleVector("expected distance averaged");
  const std::vector<double> expectedSpeedsAveraged =
      conf.getDoubleVector("expected speed averaged");
  const std::vector<float> expectedAngles = conf.getFloatVector("expected angle");
  const std::vector<int> expectedShort = conf.getIntVector("expected short");
  const std::vector<int> expectedFast = conf.getIntVector("expected fast");
  const std::vector<int> expectedBends = conf.getIntVector("expected bends");
  const std::vector<double> expectedSumSpeeds = conf.getDoubleVector("expected sum speed");
  const std::vector<double> expectedMeanSpeeds = conf.getDoubleVector("expected mean speed");
  std::vector<double> calculatedDistances, calculatedSpeeds, calculatedDistancesAveraged,
      calculatedSpeedsAveraged, calculatedSumSpeeds, calculatedMeanSpeeds;
  std::vector<float> calculatedAngles;
  std::vector<int> calculatedShort, calculatedFast, calculatedBends;
  for (auto const& tracks : diagnostics->getInitialCalculationResults()) {
    for (auto const& obsStats : tracks.first) {
      calculatedDistances.push_back(obsStats.distance);
      calculatedSpeeds.push_back(obsStats.speed);
      calculatedDistancesAveraged.push_back(obsStats.distanceAveraged);
      calculatedSpeedsAveraged.push_back(obsStats.speedAveraged);
      calculatedAngles.push_back(obsStats.angle);
    }
    calculatedShort.push_back(tracks.second.numShort_);
    calculatedFast.push_back(tracks.second.numFast_);
    calculatedBends.push_back(tracks.second.numBends_);
    calculatedSumSpeeds.push_back(tracks.second.sumSpeed_);
    calculatedMeanSpeeds.push_back(tracks.second.meanSpeed_);
  }
  EXPECT(oops::are_all_close_relative(calculatedDistances, expectedDistances,
                                      .05));
  EXPECT(oops::are_all_close_relative(calculatedSpeeds, expectedSpeeds,
                                      .05));
  EXPECT(oops::are_all_close_relative(
           calculatedDistancesAveraged, expectedDistancesAveraged,
           .05));
  EXPECT(oops::are_all_close_relative(
           calculatedSpeedsAveraged, expectedSpeedsAveraged,
           .05));
  EXPECT(oops::are_all_close_absolute(
           calculatedAngles, expectedAngles,
           5.0f));
  EXPECT_EQUAL(calculatedShort, expectedShort);
  EXPECT_EQUAL(calculatedFast, expectedFast);
  EXPECT_EQUAL(calculatedBends, expectedBends);
  EXPECT(oops::are_all_close_relative(
           calculatedSumSpeeds, expectedSumSpeeds,
           .05));
  EXPECT(oops::are_all_close_relative
         (calculatedMeanSpeeds, expectedMeanSpeeds,
          .05));
}

void testEarlyBreakCondition(const eckit::LocalConfiguration &conf) {
  const eckit::LocalConfiguration filterConf(conf, "Ship Track Check");
  if (!filterConf.getBool("early break check", true)) {
    return;
  }

  auto diagnostics = setupRunFilter(conf);
  if (!diagnostics)
    return;

  std::vector<int> expectedEarlyBreaks(conf.getIntVector("expected early breaks"));
  std::vector<int> calculatedEarlyBreaks;

  for (auto const& earlyBreakResults : diagnostics->getEarlyBreaks()) {
    calculatedEarlyBreaks.push_back(earlyBreakResults);
  }
  EXPECT_EQUAL(expectedEarlyBreaks, calculatedEarlyBreaks);
}

void testRejectedObservations(const eckit::LocalConfiguration &conf) {
  std::vector<size_t> rejectedObsIndices;
  if (setupRunFilter(conf, &rejectedObsIndices))
    return;

  const std::vector<size_t> expectedRejectedObsIndices =
      conf.getUnsignedVector("expected rejected obs indices");
  EXPECT_EQUAL(rejectedObsIndices, expectedRejectedObsIndices);
}

class TrackCheckShip : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::TrackCheckShip";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/TrackCheckShip/" + testCaseName, testCaseConf)
      {
                        testTrackCheckShipInitialCalculations(testCaseConf);
                        testEarlyBreakCondition(testCaseConf);
                        testRejectedObservations(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_TRACKCHECKSHIP_H_
