/*
 * (C) British Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/FindNearestNeighbors.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/GeometryCalculations.h"

namespace ufo {

// Custom hash function for pairs - used for quickly finding unique lat/lon
// pairs and for caching query points
struct PairHash {
  template <typename T1, typename T2>
  std::size_t operator()(const std::pair<T1, T2> &p) const {
    // boost::hash_combine style mix (see
    // https://www.boost.org/doc/libs/1_34_1/doc/html/boost/hash_combine.html)
    const auto h1 = std::hash<T1>{}(p.first);
    const auto h2 = std::hash<T2>{}(p.second);
    auto seed = h1;
    seed ^= h2 + 0x9e3779b9UL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

// -----------------------------------------------------------------------------

FindNearestNeighbors::FindNearestNeighbors(ioda::ObsSpace &obsdb,
                                           const Parameters_ &parameters,
                                           ioda::ObsDataVector<int> &flags,
                                           ioda::ObsDataVector<float> &obserr)
    : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters) {
  oops::Log::trace() << "FindNearestNeighbors constructor" << std::endl;
}

// -----------------------------------------------------------------------------

FindNearestNeighbors::~FindNearestNeighbors() {
  oops::Log::trace() << "FindNearestNeighbors destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void FindNearestNeighbors::applyFilter(
    const std::vector<bool> &apply, const Variables & /* filtervars */,
    std::vector<std::vector<bool>> & /* flagged */) const {
  oops::Log::trace() << "FindNearestNeighbors applyFilter start" << std::endl;
  // Handle output-assignment variable data types
  const std::string group = parameters_.outputAssignment.value().group();
  const std::string variable = parameters_.outputAssignment.value().variable();
  const ioda::ObsDtype dtype = obsdb_.dtype(group, variable);
  if (dtype == ioda::ObsDtype::Float) {
    std::vector<float> outputAssignmentData;
    obsdb_.get_db(group, variable, outputAssignmentData);
    process<float>(outputAssignmentData, apply);
  } else if (dtype == ioda::ObsDtype::Integer) {
    std::vector<int> outputAssignmentData;
    obsdb_.get_db(group, variable, outputAssignmentData);
    process<int>(outputAssignmentData, apply);
  } else if (dtype == ioda::ObsDtype::Integer_64) {
    std::vector<int64_t> outputAssignmentData;
    obsdb_.get_db(group, variable, outputAssignmentData);
    process<int64_t>(outputAssignmentData, apply);
  } else if (dtype == ioda::ObsDtype::String) {
    std::vector<std::string> outputAssignmentData;
    obsdb_.get_db(group, variable, outputAssignmentData);
    process<std::string>(outputAssignmentData, apply);
  } else if (dtype == ioda::ObsDtype::DateTime) {
    std::vector<util::DateTime> outputAssignmentData;
    obsdb_.get_db(group, variable, outputAssignmentData);
    process<util::DateTime>(outputAssignmentData, apply);
  } else if (dtype == ioda::ObsDtype::Bool) {
    std::vector<bool> outputAssignmentData;
    obsdb_.get_db(group, variable, outputAssignmentData);
    process<bool>(outputAssignmentData, apply);
  } else {
    const std::string dtypeStr =
        (dtype == ioda::ObsDtype::None) ? "None" : "Empty";
    throw eckit::UserError(
        "Unsupported type for output assignment variable: " + dtypeStr, Here());
  }
  oops::Log::trace() << "FindNearestNeighbors applyFilter complete"
                     << std::endl;
}

template <typename T>
void FindNearestNeighbors::process(const std::vector<T> &outputAssignmentData,
                                   const std::vector<bool> &apply) const {
  validateParameterValues();
  const size_t nLocsThisRank = obsdb_.nlocs();
  const size_t numNeighbors = parameters_.outputVariables.value().size();

  std::vector<float> lats;
  obsdb_.get_db("MetaData", "latitude", lats);
  std::vector<float> lons;
  obsdb_.get_db("MetaData", "longitude", lons);
  assert(lats.size() == nLocsThisRank);
  assert(lons.size() == nLocsThisRank);

  // The patchObs method tells us if an observation is owned by this MPI rank.
  // REMINDER: For every variable in the globalobs space, there are locations
  // that will be owned by a particular MPI rank, and those locations are the
  // same for EVERY variable (obs vector) in the rank-local obs space.
  std::vector<bool> isOwnedByThisRank(nLocsThisRank);
  obsdb_.distribution()->patchObs(isOwnedByThisRank);

  // Filter out values which are missing or aren't meant to be applied from
  // the query points
  const auto[validQueryLats, validQueryLons, validQueryIndices] =
      extractValidQueryPoints(apply, nLocsThisRank, lats, lons);

  // Filter out values which are not owned by this rank, are missing or aren't
  // meant to be applied from the reference point locations
  const auto[validReferencePointsOutputAssignmentData, validReferenceLats,
              validReferenceLons] =
      extractValidReferencePoints<T>(apply, nLocsThisRank, isOwnedByThisRank,
                                     lats, lons, outputAssignmentData);

  // Gather valid rank-local reference points from all MPI ranks.
  const auto[globalValidReferencePointsOutputAssignmentData,
              globalValidReferenceLats, globalValidReferenceLons] =
      gatherGlobalReferencePoints<T>(validReferencePointsOutputAssignmentData,
                                     validReferenceLats, validReferenceLons);

  // Avoid unnecessary extra computation by removing any duplicate locations
  const auto[uniqueGlobalValidReferencePointsOutputAssignmentData,
              uniqueGlobalValidReferenceLats, uniqueGlobalValidReferenceLons] =
      removeDuplicateReferencePoints<T>(
          globalValidReferencePointsOutputAssignmentData,
          globalValidReferenceLats, globalValidReferenceLons);

  // Now find the nearest neighbors for each query point
  const auto nearestNeighbors = computeNearestNeighbors<T>(
      {validQueryLats, validQueryLons},
      {uniqueGlobalValidReferenceLats, uniqueGlobalValidReferenceLons},
      uniqueGlobalValidReferencePointsOutputAssignmentData, numNeighbors,
      parameters_.distanceMethod.value());

  // reconstruct vectors of distances and reference values with missing values
  // where the query points were missing or the filter wasn't applied
  const auto[distances, neighborValues] =
      reconstructNeighborValuesAndDistances<T>(
          nearestNeighbors, validQueryIndices, nLocsThisRank, numNeighbors);

  // Output distances and nearest neighbor values
  const std::string units = parameters_.distanceUnits.value();
  const std::vector<Variable> &distanceOutputVars =
      parameters_.distanceOutputVariables.value();
  const std::vector<Variable> &outputVars = parameters_.outputVariables.value();
  for (size_t i = 0; i < numNeighbors; ++i) {
    // Convert distances from metres to the requested units before writing
    std::vector<float> convertedDistances = distances[i];
    for (float &distance : convertedDistances) {
      if (distance != util::missingValue<float>()) {
        distance = this->convertDistanceFromMeters(distance, units);
      }
    }
    obsdb_.put_db(distanceOutputVars[i].group(),
                  distanceOutputVars[i].variable(), convertedDistances);
    obsdb_.put_db(outputVars[i].group(), outputVars[i].variable(),
                  neighborValues[i]);
  }
  oops::Log::trace() << "FindNearestNeighbors process complete" << std::endl;
}

// -----------------------------------------------------------------------------

void FindNearestNeighbors::validateParameterValues() const {
  const size_t numNeighbors = parameters_.outputVariables.value().size();
  if (numNeighbors == 0) {
    throw eckit::UserError("No output variables specified.", Here());
  }

  if (parameters_.distanceOutputVariables.value().size() != numNeighbors) {
    throw eckit::UserError(
        "Number of distance output variables must match number of output "
        "variables.",
        Here());
  }

  const std::string units = parameters_.distanceUnits.value();
  if (!(units == "m" || units == "metres" || units == "meters" ||
        units == "km" || units == "kilometres" || units == "kilometers" ||
        units == "mi" || units == "miles" || units == "nmi" ||
        units == "nautical miles")) {
    throw eckit::UserError("Unsupported distance units: " + units, Here());
  }

  if (parameters_.distanceMethod.value() != DistanceMethod::HAVERSINE) {
    throw eckit::NotImplemented(
        "Only haversine distance method is implemented.", Here());
  }

  const std::vector<Variable> &outputVars = parameters_.outputVariables.value();
  const std::vector<Variable> &distanceOutputVars =
      parameters_.distanceOutputVariables.value();
  for (size_t i = 0; i < numNeighbors; ++i) {
    if (outputVars[i].group() == "ObsValue" ||
        outputVars[i].group() == "DerivedObsValue") {
      throw eckit::BadParameter(
          "Output variables cannot be in group ObsValue or DerivedObsValue, "
          "since qc flag information cannot be updated. "
          "Variable " +
              outputVars[i].variable() + " is in group " +
              outputVars[i].group(),
          Here());
    }
    if (distanceOutputVars[i].group() == "ObsValue" ||
        distanceOutputVars[i].group() == "DerivedObsValue") {
      throw eckit::BadParameter(
          "Distance output variables cannot be in group ObsValue or "
          "DerivedObsValue, "
          "since qc flag information cannot be updated. "
          "Variable " +
              distanceOutputVars[i].variable() + " is in group " +
              distanceOutputVars[i].group(),
          Here());
    }
  }

  if (parameters_.algorithm.value() != "brute force") {
    throw eckit::NotImplemented("Only brute force algorithm is implemented.",
                                Here());
  }

  if (!obsdb_.has(parameters_.outputAssignment.value().group(),
                  parameters_.outputAssignment.value().variable())) {
    throw eckit::UserError("Output assignment variable " +
                               parameters_.outputAssignment.value().group() +
                               "/" +
                               parameters_.outputAssignment.value().variable() +
                               " not found in ObsSpace.",
                           Here());
  }

  // Check that the query and reference variables exist and are floats
  const std::string queryPointVarName =
      parameters_.queryPointVar.value().variable();
  const std::string queryPointVarGroup =
      parameters_.queryPointVar.value().group();
  const std::string referencePointVarGroup =
      parameters_.referencePointVar.value().group();
  const std::string referencePointVarName =
      parameters_.referencePointVar.value().variable();
  if (!obsdb_.has(queryPointVarGroup, queryPointVarName)) {
    throw eckit::UserError("Query point variable " + queryPointVarGroup + "/" +
                               queryPointVarName + " not found in ObsSpace.",
                           Here());
  }
  if (obsdb_.dtype(queryPointVarGroup, queryPointVarName) !=
      ioda::ObsDtype::Float) {
    throw eckit::UserError("Query point variable " + queryPointVarGroup + "/" +
                               queryPointVarName + " is not a float.",
                           Here());
  }
  if (!obsdb_.has(referencePointVarGroup, referencePointVarName)) {
    throw eckit::UserError(
        "Reference point variable " + referencePointVarGroup + "/" +
            referencePointVarName + " not found in ObsSpace.",
        Here());
  }
  if (obsdb_.dtype(referencePointVarGroup, referencePointVarName) !=
      ioda::ObsDtype::Float) {
    throw eckit::UserError("Reference point variable " +
                               referencePointVarGroup + "/" +
                               referencePointVarName + " is not a float.",
                           Here());
  }
}

// -----------------------------------------------------------------------------

template <typename T>
std::vector<std::vector<FindNearestNeighbors::Neighbor<T>>>
FindNearestNeighbors::computeNearestNeighbors(
    const std::vector<std::vector<float>> &queryCoords,
    const std::vector<std::vector<float>> &referenceCoords,
    const std::vector<T> &referenceValues, size_t numNeighbors,
    const DistanceMethod distanceMethod) const {
  const size_t nQuery = queryCoords[0].size();
  const size_t nReference = referenceCoords[0].size();

  // Result: A vector for each query point (hence this being a vector of
  // vectors), where the inner vector is in distance order from closest to
  // furthest
  std::vector<std::vector<Neighbor<T>>> allQueryPointNeighbors(nQuery);

  // Query points are lat/lon pairs
  using QueryPointType = std::pair<float, float>;

  // A map of query lat/lon pairs to each point's vector of nearest neighbors is
  // used to cache results for identical query points and avoid redundant
  // distance calculations and sorting.
  std::unordered_map<QueryPointType, std::vector<Neighbor<T>>, PairHash>
      queryPointToNeighborMap;

  for (size_t i = 0; i < nQuery; ++i) {
    const QueryPointType queryPoint{queryCoords[0][i], queryCoords[1][i]};

    // Check if we've already computed nearest neighbors for this query point
    const auto cachedNeighbor = queryPointToNeighborMap.find(queryPoint);
    if (cachedNeighbor != queryPointToNeighborMap.end()) {
      // If so, use the cached result instead of recomputing
      allQueryPointNeighbors[i] = cachedNeighbor->second;
      continue;
    }

    // Calculate distances for this query point
    std::vector<Neighbor<T>> queryPointNeighbors(nReference);
    for (size_t j = 0; j < nReference; ++j) {
      const float distance = calculateDistanceInMetres(
          queryCoords[0][i], queryCoords[1][i], referenceCoords[0][j],
          referenceCoords[1][j], distanceMethod);
      queryPointNeighbors[j] = Neighbor<T>{distance, referenceValues[j]};
    }

    // Sort the distances vector by the `distance` member of the struct
    std::sort(
        queryPointNeighbors.begin(), queryPointNeighbors.end(),
        [](const auto &a, const auto &b) { return a.distance < b.distance; });

    // Resize the distances vector to keep only the nearest numNeighbors
    if (queryPointNeighbors.size() > numNeighbors) {
      queryPointNeighbors.resize(numNeighbors);
    }

    // Store the computed result so identical query points can skip the work.
    queryPointToNeighborMap.emplace(queryPoint, queryPointNeighbors);

    // Assign the nearest neighbors to the result for this query point
    allQueryPointNeighbors[i] = std::move(queryPointNeighbors);
  }

  return allQueryPointNeighbors;
}

// -----------------------------------------------------------------------------

float FindNearestNeighbors::calculateDistanceInMetres(
    const float lat1, const float lon1, const float lat2, const float lon2,
    const DistanceMethod method) const {
  if (method == DistanceMethod::HAVERSINE) {
    return static_cast<float>(ufo::haversine(lat1, lon1, lat2, lon2));
  } else {
    throw eckit::UserError("Unsupported distance method enum value.", Here());
  }
}

float FindNearestNeighbors::convertDistanceFromMeters(
    const float distance, const std::string &units) const {
  if (units == "m" || units == "metres" || units == "meters") {
    return distance;
  } else if (units == "km" || units == "kilometres" || units == "kilometers") {
    return distance / 1000.0f;
  } else if (units == "mi" || units == "miles") {
    return distance / ufo::Constants::m_per_mile;
  } else if (units == "nmi" || units == "nautical miles") {
    return distance / ufo::Constants::m_per_nautical_mile;
  } else {
    throw eckit::UserError("Unsupported distance units: " + units, Here());
  }
}

// -----------------------------------------------------------------------------

void FindNearestNeighbors::print(std::ostream &os) const {
  os << "FindNearestNeighbors: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>, std::vector<size_t>>
FindNearestNeighbors::extractValidQueryPoints(
    const std::vector<bool> &apply, const size_t nLocsThisRank,
    const std::vector<float> &lats, const std::vector<float> &lons) const {
  assert(lats.size() == nLocsThisRank);
  assert(lons.size() == nLocsThisRank);
  const float missingFloat = util::missingValue<float>();

  std::vector<float> queryPointVariableVals;
  obsdb_.get_db(parameters_.queryPointVar.value().group(),
                parameters_.queryPointVar.value().variable(),
                queryPointVariableVals);
  assert(queryPointVariableVals.size() == nLocsThisRank);

  std::vector<float> validQueryLats;
  std::vector<float> validQueryLons;
  std::vector<size_t> validQueryIndices;

  for (size_t i = 0; i < nLocsThisRank; ++i) {
    if (apply[i] && queryPointVariableVals[i] != missingFloat &&
        lats[i] != missingFloat && lons[i] != missingFloat) {
      validQueryLats.push_back(lats[i]);
      validQueryLons.push_back(lons[i]);
      validQueryIndices.push_back(i);
    }
  }

  return {validQueryLats, validQueryLons, validQueryIndices};
}

template <typename T>
std::tuple<std::vector<T>, std::vector<float>, std::vector<float>>
FindNearestNeighbors::extractValidReferencePoints(
    const std::vector<bool> &apply, const size_t nLocsThisRank,
    const std::vector<bool> &isOwnedByThisRank, const std::vector<float> &lats,
    const std::vector<float> &lons,
    const std::vector<T> &outputAssignmentData) const {
  assert(lats.size() == nLocsThisRank);
  assert(lons.size() == nLocsThisRank);
  assert(isOwnedByThisRank.size() == nLocsThisRank);
  assert(outputAssignmentData.size() == nLocsThisRank);
  const float missingFloat = util::missingValue<float>();

  std::vector<float> referencePointVariableVals;
  obsdb_.get_db(parameters_.referencePointVar.value().group(),
                parameters_.referencePointVar.value().variable(),
                referencePointVariableVals);
  assert(referencePointVariableVals.size() == nLocsThisRank);

  std::vector<T> validReferencePointsOutputAssignmentData;
  std::vector<float> validReferenceLats;
  std::vector<float> validReferenceLons;
  for (size_t i = 0; i < nLocsThisRank; ++i) {
    if (apply[i] && isOwnedByThisRank[i]) {
      if (referencePointVariableVals[i] != missingFloat &&
          lats[i] != missingFloat && lons[i] != missingFloat) {
        validReferencePointsOutputAssignmentData.push_back(
            outputAssignmentData[i]);
        validReferenceLats.push_back(lats[i]);
        validReferenceLons.push_back(lons[i]);
      }
    }
  }
  return {validReferencePointsOutputAssignmentData, validReferenceLats,
          validReferenceLons};
}

template <typename T>
std::tuple<std::vector<T>, std::vector<float>, std::vector<float>>
FindNearestNeighbors::gatherGlobalReferencePoints(
    const std::vector<T> &localData, const std::vector<float> &localLats,
    const std::vector<float> &localLons) const {
  assert(localData.size() == localLats.size());
  assert(localData.size() == localLons.size());
  assert(localLats.size() == localLons.size());
  // Global vectors are initialized with local data and modified in place by
  // allGatherv
  std::vector<T> globalData = localData;
  std::vector<float> globalLats = localLats;
  std::vector<float> globalLons = localLons;
  // special handling for bool vectors since vector<bool> is not a standard
  // container and does not have a proper data() method which allGatherv needs
  // to access the underlying data.
  if constexpr (std::is_same<T, bool>::value) {
    std::vector<char> boolAsBytes(globalData.begin(), globalData.end());
    oops::mpi::allGatherv(obsdb_.comm(), boolAsBytes);
    globalData.assign(boolAsBytes.begin(), boolAsBytes.end());
  } else {
    oops::mpi::allGatherv(obsdb_.comm(), globalData);
  }
  oops::mpi::allGatherv(obsdb_.comm(), globalLats);
  oops::mpi::allGatherv(obsdb_.comm(), globalLons);

  assert(globalData.size() == globalLats.size());
  assert(globalData.size() == globalLons.size());
  assert(globalLats.size() == globalLons.size());

  return {globalData, globalLats, globalLons};
}

std::vector<size_t> FindNearestNeighbors::uniqueLatLonIndices(
    const std::vector<float> &lats, const std::vector<float> &lons) const {
  std::unordered_set<std::pair<float, float>, PairHash> uniqueLatLonSet;
  std::vector<size_t> uniqueIndices;

  assert(lats.size() == lons.size());
  for (size_t i = 0; i < lats.size(); ++i) {
    const auto latLon = std::make_pair(lats[i], lons[i]);
    if (uniqueLatLonSet.insert(latLon).second) {
      // If insertion succeeds, the pair is unique
      uniqueIndices.push_back(i);
    }
  }
  return uniqueIndices;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<float>, std::vector<float>>
FindNearestNeighbors::removeDuplicateReferencePoints(
    const std::vector<T> &data, const std::vector<float> &lats,
    const std::vector<float> &lons) const {
  assert(data.size() == lats.size());
  assert(data.size() == lons.size());
  assert(lats.size() == lons.size());

  std::vector<size_t> uniqueIndices = uniqueLatLonIndices(lats, lons);
  std::vector<float> uniqueLats;
  std::vector<float> uniqueLons;
  std::vector<T> uniqueData;
  for (size_t idx : uniqueIndices) {
    uniqueLats.push_back(lats[idx]);
    uniqueLons.push_back(lons[idx]);
    uniqueData.push_back(data[idx]);
  }
  return {uniqueData, uniqueLats, uniqueLons};
}

template <typename T>
std::pair<std::vector<std::vector<float>>, std::vector<std::vector<T>>>
FindNearestNeighbors::reconstructNeighborValuesAndDistances(
    const std::vector<std::vector<Neighbor<T>>> &nearestNeighbors,
    const std::vector<size_t> &validQueryIndices, const size_t nLocsThisRank,
    const size_t numNeighbors) const {
  // Initialize with missing values
  std::vector<std::vector<float>> distances(
      numNeighbors,
      std::vector<float>(nLocsThisRank, util::missingValue<float>()));
  std::vector<std::vector<T>> neighborValues(
      numNeighbors, std::vector<T>(nLocsThisRank, util::missingValue<T>()));

  // Fill in valid values
  for (size_t validIdx = 0; validIdx < validQueryIndices.size(); ++validIdx) {
    const size_t originalIdx = validQueryIndices[validIdx];
    const size_t numAvailableNeighbors = nearestNeighbors[validIdx].size();
    const size_t numNeighborsToAssign =
        std::min(numNeighbors, numAvailableNeighbors);
    for (size_t j = 0; j < numNeighborsToAssign; ++j) {
      distances[j][originalIdx] = nearestNeighbors[validIdx][j].distance;
      neighborValues[j][originalIdx] =
          nearestNeighbors[validIdx][j].referenceValue;
    }
  }
  return {distances, neighborValues};
}

}  // namespace ufo
