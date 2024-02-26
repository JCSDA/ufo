/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_TRACKCHECKSHIPMAINLOOP_H_
#define TEST_UFO_TRACKCHECKSHIPMAINLOOP_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/TrackCheckShip.h"
#include "ufo/filters/TrackCheckShipDiagnostics.h"
#include "ufo/filters/Variables.h"

namespace ufo {
namespace test {

void testFirstRejectionSimultaneousIncluded(const eckit::LocalConfiguration &conf) {
  util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsspace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  if (conf.has("station_ids")) {
    const std::vector<int> stationIds = conf.getIntVector("station_ids");
    obsspace.put_db("MetaData", "stationIdentification", stationIds);
  }

  std::shared_ptr<ioda::ObsDataVector<float>> obserr(new ioda::ObsDataVector<float>(
      obsspace, obsspace.obsvariables(), "ObsError"));
  std::shared_ptr<ioda::ObsDataVector<int>> qcflags(new ioda::ObsDataVector<int>(
      obsspace, obsspace.obsvariables()));

  eckit::LocalConfiguration filterConf(conf, "Ship Track Check");
  ufo::TrackCheckShipParameters filterParameters;
  filterParameters.validateAndDeserialize(filterConf);
  ufo::TrackCheckShip filter(obsspace, filterParameters, qcflags, obserr);
  filter.preProcess();
  const ufo::TrackCheckShipDiagnostics diagnostics = *filter.diagnostics();

  const std::vector<size_t> expectedFirstRejectionIndices =
      conf.getUnsignedVector("expected first rejection index");
  const std::vector<int> expectedCategories = conf.getIntVector("expected error category");
  std::vector<size_t> testedFirstRejectionIndices;
  std::vector<int> testedErrorCategories;
  for (auto const& trackFirstRemovalsInfo : diagnostics.getFirstIterativeRemovalInfo()) {
    for (auto const& indexRejected : trackFirstRemovalsInfo.first)
      testedFirstRejectionIndices.push_back(indexRejected);
    testedErrorCategories.push_back(trackFirstRemovalsInfo.second);
  }
  EXPECT_EQUAL(testedFirstRejectionIndices, expectedFirstRejectionIndices);
  EXPECT_EQUAL(testedErrorCategories, expectedCategories);
  auto const expectedDistanceSum = conf.getDoubleVector("expected distance sum", {-1});
  if (*expectedDistanceSum.begin() != -1) {
    auto calculatedDistanceSum = diagnostics.getDistanceSum();
    EXPECT(oops::are_all_close_relative(calculatedDistanceSum, expectedDistanceSum, .05));
  }
  auto const expectedDistancePreviousObservationOmitted =
      conf.getDoubleVector("expected distance previous observation omitted", {-1});
  if (*expectedDistancePreviousObservationOmitted.begin() != -1) {
    auto calculatedDistancePreviousObservationOmitted = diagnostics.getDistancePrevObsOmitted();
    EXPECT(oops::are_all_close_relative(calculatedDistancePreviousObservationOmitted,
                                        expectedDistancePreviousObservationOmitted, .05));
  }
  auto const expectedDistanceCurrentObservationOmitted =
      conf.getDoubleVector("expected distance current observation omitted", {-1});
  if (*expectedDistanceCurrentObservationOmitted.begin() != -1) {
    auto calculatedDistanceCurrentObservationOmitted = diagnostics.getDistanceCurrentObsOmitted();
    EXPECT(oops::are_all_close_relative(calculatedDistanceCurrentObservationOmitted,
                                        expectedDistanceCurrentObservationOmitted, .05));
  }
  auto const expectedTimeSum = conf.getDoubleVector("expected time sum", {-1});
  if (*expectedTimeSum.begin() != -1) {
    auto calculatedTimeSum = diagnostics.getTimeSum();
    EXPECT(oops::are_all_close_relative(calculatedTimeSum, expectedTimeSum, .05));
  }
  auto const expectedPreviousSegmentDistanceProportion = conf.getDoubleVector(
        "expected previous segment distance proportion", {-1});
  if (*expectedPreviousSegmentDistanceProportion.begin() != -1) {
    auto calculatedPreviousSegmentDistanceProportion =
        diagnostics.getPreviousSegmentDistanceProportion();
    EXPECT(oops::are_all_close_relative(calculatedPreviousSegmentDistanceProportion,
                                        expectedPreviousSegmentDistanceProportion, .05));
  }
  auto const expectedPreviousObservationDistanceAveragedProportion = conf.getDoubleVector(
        "expected previous observation distance averaged proportion", {-1});
  if (*expectedPreviousObservationDistanceAveragedProportion.begin() != -1) {
    auto calculatedPreviousObservationDistanceAveragedProportion =
        diagnostics.getPreviousObservationDistanceAveragedProportion();
    EXPECT(oops::are_all_close_relative(
             calculatedPreviousObservationDistanceAveragedProportion,
             expectedPreviousObservationDistanceAveragedProportion, .05));
  }
  auto const expectedPreviousSegmentTimeProportion = conf.getDoubleVector(
        "expected previous segment time proportion", {-1});
  if (*expectedPreviousSegmentTimeProportion.begin() != -1) {
    auto calculatedPreviousSegmentTimeProportion =
        diagnostics.getPreviousSegmentTimeProportion();
    EXPECT(oops::are_all_close_relative(calculatedPreviousSegmentTimeProportion,
                                        expectedPreviousSegmentTimeProportion, .05));
  }
  auto const expectedPreviousAndFastestSegmentTimeProportion = conf.getDoubleVector(
        "expected previous and fastest segment time proportion", {-1});
  if (*expectedPreviousAndFastestSegmentTimeProportion.begin() != -1) {
    auto calculatedPreviousAndFastestSegmentTimeProportion =
        diagnostics.getPreviousAndFastestSegmentTimeProportion();
    EXPECT(oops::are_all_close_relative(calculatedPreviousAndFastestSegmentTimeProportion,
                                        expectedPreviousAndFastestSegmentTimeProportion, .05));
  }

  for (auto const& i : expectedFirstRejectionIndices) {
    EXPECT_EQUAL((*qcflags)[0][i], ufo::QCflags::track);
    // tests that the rejected observations are changing the qc flags correctly
  }
}


class TrackCheckShipMainLoop : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::TrackCheckShipMainLoop";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
    for (const std::string & testCaseName : conf.keys())
    {
      const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(), testCaseName);
      ts.emplace_back(CASE("ufo/TrackCheckShipMainLoop/" + testCaseName, testCaseConf)
                      {
                        testFirstRejectionSimultaneousIncluded(testCaseConf);
                      });
    }
  }

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_TRACKCHECKSHIPMAINLOOP_H_
