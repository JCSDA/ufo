/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_METOFFICEBUDDYPAIRFINDER_H_
#define TEST_UFO_METOFFICEBUDDYPAIRFINDER_H_

#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/MetOfficeBuddyCheckParameters.h"
#include "ufo/filters/MetOfficeBuddyPairFinder.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {
namespace test {

typedef std::pair<int, int> ObsPair;

bool duplicatePairsPresent(const std::vector<ObsPair> &pairs) {
  bool duplicatesFound = false;
  std::set<ObsPair> pairSet;

  // Test for exact duplicates
  for (const ObsPair &pair : pairs) {
    bool inserted = pairSet.insert(pair).second;
    if (!inserted) {
      duplicatesFound = true;
      oops::Log::trace() << "Duplicate pair found: (" << pair.first << ", " << pair.second << ")\n";
    }
  }

  // Test for duplicates in reverse order
  for (const ObsPair &pair : pairs) {
    ObsPair reversePair(pair.second, pair.first);
    bool inserted = pairSet.insert(reversePair).second;
    if (!inserted) {
      duplicatesFound = true;
      oops::Log::trace() << "Pair (" << pair.first << ", " << pair.second
                         << ") present also in reverse order\n";
    }
  }

  oops::Log::trace().flush();

  return duplicatesFound;
}

template <typename Key, typename Value>
Value maxValue(const std::map<Key, Value> &map)
{
  if (map.empty())
    return Value();

  typedef typename std::map<Key, Value>::value_type ValueType;
  return std::max_element(map.begin(), map.end(),
                          [](const ValueType &a, const ValueType &b)
                          { return a.second < b.second; })->second;
}

int maxTotalNumBuddies(const std::vector<ObsPair> &pairs) {
  std::map<int, int> numBuddiesByObsId;
  for (const ObsPair & pair : pairs)
    ++numBuddiesByObsId[pair.first];
  return maxValue(numBuddiesByObsId);
}

int maxNumBuddiesWithSameStationId(const std::vector<ObsPair> &pairs,
                                        const std::vector<int> &stationIds) {
  std::map<int, int> numBuddiesWithSameStationIdByObsId;
  for (const ObsPair & pair : pairs)
    if (stationIds[pair.first] == stationIds[pair.second])
      ++numBuddiesWithSameStationIdByObsId[pair.first];
  return maxValue(numBuddiesWithSameStationIdByObsId);
}

void testDuplicatesAndBuddyCountConstraints(const eckit::LocalConfiguration &conf) {
  const util::TimeWindow timeWindow(conf.getSubConfiguration("time window"));

  const eckit::LocalConfiguration obsSpaceConf(conf, "obs space");
  ioda::ObsSpace obsSpace(obsSpaceConf, oops::mpi::world(), timeWindow, oops::mpi::myself());

  boost::optional<std::vector<float>> airPressures;
  if (obsSpace.has("MetaData", "pressure")) {
    airPressures = std::vector<float>(obsSpace.nlocs());
    obsSpace.get_db("MetaData", "pressure", *airPressures);
  }

  std::vector<float> latitudes(obsSpace.nlocs());
  obsSpace.get_db("MetaData", "latitude", latitudes);

  std::vector<float> longitudes(obsSpace.nlocs());
  obsSpace.get_db("MetaData", "longitude", longitudes);

  std::vector<util::DateTime> datetimes(obsSpace.nlocs());
  obsSpace.get_db("MetaData", "dateTime", datetimes);

  std::vector<int> stationIds(obsSpace.recnum().begin(), obsSpace.recnum().end());

  std::vector<size_t> validObsIds;
  if (conf.has("valid_obs_ids")) {
    validObsIds = conf.getUnsignedVector("valid_obs_ids");
  } else {
    validObsIds.resize(obsSpace.nlocs());
    std::iota(validObsIds.begin(), validObsIds.end(), 0);
  }

  const eckit::LocalConfiguration modernFilterConf(conf, "Met Office Buddy Check, modern");
  MetOfficeBuddyCheckParameters modernOptions;
  modernOptions.deserialize(modernFilterConf);

  std::vector<ObsPair> modernPairs;
  {
    MetOfficeBuddyPairFinder finder(modernOptions, latitudes, longitudes, datetimes,
                                    airPressures.get_ptr(), stationIds);
    const std::vector<MetOfficeBuddyPair> buddyPairs = finder.findBuddyPairs(validObsIds);

    for (const MetOfficeBuddyPair & pair : buddyPairs) {
      modernPairs.push_back(ObsPair(pair.obsIdA, pair.obsIdB));
    }
  }

  EXPECT(maxTotalNumBuddies(modernPairs) <= modernOptions.maxTotalNumBuddies);
  EXPECT(maxNumBuddiesWithSameStationId(modernPairs, stationIds) <=
         modernOptions.maxNumBuddiesWithSameStationId);
  EXPECT_NOT(duplicatePairsPresent(modernPairs));

  const eckit::LocalConfiguration legacyFilterConf(conf, "Met Office Buddy Check, legacy");
  MetOfficeBuddyCheckParameters legacyOptions;
  legacyOptions.deserialize(legacyFilterConf);

  std::vector<ObsPair> legacyPairs;
  {
    MetOfficeBuddyPairFinder finder(legacyOptions, latitudes, longitudes, datetimes,
                                    airPressures.get_ptr(), stationIds);
    const std::vector<MetOfficeBuddyPair> buddyPairs = finder.findBuddyPairs(validObsIds);

    for (const MetOfficeBuddyPair & pair : buddyPairs) {
      legacyPairs.push_back(ObsPair(pair.obsIdA, pair.obsIdB));
    }
  }

  EXPECT(maxTotalNumBuddies(legacyPairs) <= legacyOptions.maxTotalNumBuddies);
  EXPECT(maxNumBuddiesWithSameStationId(legacyPairs, stationIds) <=
         legacyOptions.maxNumBuddiesWithSameStationId);
  EXPECT_NOT(duplicatePairsPresent(legacyPairs));

  // The legacy buddy collector sometimes fails to find all observations that can be classified
  // as buddies while respecting constraints on their maximum number.
  EXPECT(modernPairs.size() > legacyPairs.size());
}

CASE("ufo/MetOfficeBuddyPairFinder/"
     "Duplicates, constraints on buddy counts, legacy pair collector") {
  testDuplicatesAndBuddyCountConstraints(eckit::LocalConfiguration(
                                           ::test::TestEnvironment::config(),
                                           "Duplicates, constraints on buddy counts, "
                                           "legacy pair collector"));
}

std::vector<MetOfficeBuddyPair> findBuddyPairs(const MetOfficeBuddyCheckParameters &options,
                                               const std::vector<float> &latitudes,
                                               const std::vector<float> &longitudes,
                                               const std::vector<util::DateTime> &datetimes,
                                               const std::vector<float> *pressures,
                                               const std::vector<int> &stationIds,
                                               const std::vector<size_t> &validObsIds) {
  MetOfficeBuddyPairFinder finder(options, latitudes, longitudes, datetimes,
                                  pressures, stationIds);
  return finder.findBuddyPairs(validObsIds);
}

void testInvarianceToLongitude(const eckit::LocalConfiguration &conf) {
  const float searchRadius = 100;  // km

  std::vector<float> referenceLatitudes = conf.getFloatVector("reference.latitudes");
  std::vector<float> referenceLongitudes = conf.getFloatVector("reference.longitudes");

  std::vector<util::DateTime> datetimes(2, util::DateTime());
  std::vector<int> stationIds{1, 2};

  std::vector<size_t> validObsIds{0, 1};

  const eckit::LocalConfiguration filterConf(conf, "Met Office Buddy Check");
  MetOfficeBuddyCheckParameters options;
  options.deserialize(filterConf);

  const std::vector<MetOfficeBuddyPair> referenceBuddyPairs = findBuddyPairs(
        options, referenceLatitudes, referenceLongitudes,
        datetimes, nullptr, stationIds, validObsIds);
  EXPECT_EQUAL(referenceBuddyPairs.size(), 1);
  const MetOfficeBuddyPair &refPair = referenceBuddyPairs[0];

  for (const eckit::LocalConfiguration &testConf : conf.getSubConfigurations("test")) {
    std::vector<float> latitudes = testConf.getFloatVector("latitudes");
    std::vector<float> longitudes = testConf.getFloatVector("longitudes");
    const std::vector<MetOfficeBuddyPair> buddyPairs = findBuddyPairs(
          options, latitudes, longitudes, datetimes, nullptr, stationIds, validObsIds);
    EXPECT_EQUAL(buddyPairs.size(), 1);
    const MetOfficeBuddyPair &pair = buddyPairs[0];

    if (pair.obsIdA == refPair.obsIdA) {
      EXPECT_EQUAL(pair.obsIdB, refPair.obsIdB);
      EXPECT(oops::is_close_relative(pair.distanceInKm, refPair.distanceInKm, 1e-6));
      if (latitudes[pair.obsIdB] - latitudes[pair.obsIdA] ==
          referenceLatitudes[refPair.obsIdB] - referenceLatitudes[refPair.obsIdA] &&
          longitudes[pair.obsIdB] - longitudes[pair.obsIdA] ==
          referenceLongitudes[refPair.obsIdB] - referenceLongitudes[refPair.obsIdA]) {
        EXPECT(oops::is_close_relative(pair.rotationAInRad, refPair.rotationAInRad, 1e-6));
        EXPECT(oops::is_close_relative(pair.rotationBInRad, refPair.rotationBInRad, 1e-6));
      }
    } else {
      EXPECT_EQUAL(pair.obsIdA, refPair.obsIdB);
      EXPECT_EQUAL(pair.obsIdB, refPair.obsIdA);
      EXPECT(oops::is_close_relative(pair.distanceInKm, refPair.distanceInKm, 1e-6));
      if (latitudes[pair.obsIdB] - latitudes[pair.obsIdA] ==
          referenceLatitudes[refPair.obsIdA] - referenceLatitudes[refPair.obsIdB] &&
          longitudes[pair.obsIdB] - longitudes[pair.obsIdA] ==
          referenceLongitudes[refPair.obsIdA] - referenceLongitudes[refPair.obsIdB]) {
        EXPECT(oops::is_close_relative(pair.rotationAInRad, refPair.rotationBInRad, 1e-6));
        EXPECT(oops::is_close_relative(pair.rotationBInRad, refPair.rotationAInRad, 1e-6));
      }
    }
  }
}

CASE("ufo/TemporalThinning/Invariance to longitude, different zonal bands") {
  testInvarianceToLongitude(eckit::LocalConfiguration(
                              ::test::TestEnvironment::config(),
                              "Invariance to longitude, different zonal bands"));
}

CASE("ufo/TemporalThinning/Invariance to longitude, same zonal band") {
  testInvarianceToLongitude(eckit::LocalConfiguration(
                              ::test::TestEnvironment::config(),
                              "Invariance to longitude, same zonal band"));
}

void findEndpoint(double startLat, double startLon, double distance, double bearing,
                  double &endLat, double &endLon)
{
    startLat *= ufo::Constants::deg2rad;
    startLon *= ufo::Constants::deg2rad;
    bearing *= ufo::Constants::deg2rad;
    distance /= ufo::Constants::mean_earth_rad;
    endLat = std::asin(std::sin(startLat) * std::cos(distance) +
                       std::cos(startLat) * std::sin(distance) * std::cos(bearing));
    endLon = startLon + std::atan2(std::sin(bearing) * std::sin(distance) * std::cos(startLat),
                               std::cos(distance) - std::sin(startLat) * std::sin(endLat));
    endLat *= ufo::Constants::rad2deg;
    endLon *= ufo::Constants::rad2deg;
}

void findEndpoint(float startLat, float startLon, float distance, float bearing,
                  float &endLat, float &endLon)
{
  double dEndLat, dEndLon;
  findEndpoint(startLat, startLon, distance, bearing, dEndLat, dEndLon);
  endLon = dEndLon;
  endLat = dEndLat;
}

template <typename T>
bool contains(const std::set<T> &set, const T &element) {
  return set.find(element) != set.end();
}

ObsPair reverse(const ObsPair &pair) {
  return ObsPair(pair.second, pair.first);
}

void testSearchRadius(const eckit::LocalConfiguration &conf) {
  const eckit::LocalConfiguration filterConf(conf, "Met Office Buddy Check");
  MetOfficeBuddyCheckParameters options;
  options.deserialize(filterConf);

  const float searchRadius = options.searchRadius;

  for (const eckit::LocalConfiguration &testConf : conf.getSubConfigurations("test")) {
    float centralLatitude = testConf.getFloat("central_latitude");
    float centralLongitude = testConf.getFloat("central_longitude");

    std::vector<float> latitudes{centralLatitude};
    std::vector<float> longitudes{centralLongitude};
    std::vector<ObsPair> pairsExpectedToBePresent, pairsExpectedToBeAbsent;

    const int numBearings = 8;
    for (float i = 0; i < numBearings; ++i) {
      float latitude, longitude;

      findEndpoint(centralLatitude, centralLongitude,
                   0.99f * searchRadius, i * 360.0f / numBearings, latitude, longitude);
      latitudes.push_back(latitude);
      longitudes.push_back(longitude);
      pairsExpectedToBePresent.push_back(ObsPair(0, latitudes.size() - 1));

      findEndpoint(centralLatitude, centralLongitude,
                   1.01f * searchRadius, i * 360.0f / numBearings, latitude, longitude);
      latitudes.push_back(latitude);
      longitudes.push_back(longitude);
      pairsExpectedToBeAbsent.push_back(ObsPair(0, latitudes.size() - 1));
    }

    const std::vector<util::DateTime> datetimes(latitudes.size(), util::DateTime());
    const std::vector<int> stationIds(latitudes.size(), 0);

    std::vector<size_t> validObsIds(latitudes.size());
    std::iota(validObsIds.begin(), validObsIds.end(), 0);

    const std::vector<MetOfficeBuddyPair> buddyPairs = findBuddyPairs(
          options, latitudes, longitudes, datetimes, nullptr, stationIds, validObsIds);

    std::set<ObsPair> pairs;
    for (const MetOfficeBuddyPair & pair : buddyPairs) {
      pairs.insert(ObsPair(pair.obsIdA, pair.obsIdB));
    }

    for (const ObsPair &pair : pairsExpectedToBePresent)
      EXPECT(contains(pairs, pair) || contains(pairs, reverse(pair)));

    for (const ObsPair &pair : pairsExpectedToBeAbsent)
      EXPECT(!contains(pairs, pair) && !contains(pairs, reverse(pair)));
  }
}

CASE("ufo/TemporalThinning/Search radius") {
  testSearchRadius(eckit::LocalConfiguration(::test::TestEnvironment::config(), "Search radius"));
}

class MetOfficeBuddyPairFinder : public oops::Test {
 private:
  std::string testid() const override {return "ufo::test::MetOfficeBuddyPairFinder";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_METOFFICEBUDDYPAIRFINDER_H_
