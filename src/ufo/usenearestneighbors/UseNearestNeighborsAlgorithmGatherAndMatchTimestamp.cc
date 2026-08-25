/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmGatherAndMatchTimestamp.h"
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::
    UseNearestNeighborsAlgorithmGatherAndMatchTimestamp(
        const UseNearestNeighborsGatherAndMatchTimestampParameters &params,
        const ObsFilterData &data, const std::vector<bool> &apply,
        const Variables &filtervars, const ioda::ObsDataVector<int> &flags,
        std::vector<std::vector<bool>> &flagged)
    : UseNearestNeighborsAlgorithmBase(params, data, apply, filtervars, flags,
                                       flagged) {}

void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::execute(
    const UseNearestNeighborsAlgorithmParametersBase &algParams,
    const UseNearestNeighborsParameters &options) const {
  const auto &params = dynamic_cast<
      const UseNearestNeighborsGatherAndMatchTimestampParameters &>(algParams);

  // Get and verify parameters:
  const Variable idVariable = options.idVar.value();
  const std::vector<Variable> &nearestNeighborIDVars =
      options.nearestNeighborIDVars.value();
  const std::vector<Variable> &outVars = params.outputVariables.value();
  const Variable gatherVariable = params.gatherVariable.value();

  if (nearestNeighborIDVars.size() != outVars.size()) {
    throw eckit::BadParameter(
        "Number of nearest neighbor identifier variables must match number of "
        "output variables.",
        Here());
  }

  // Determine types
  auto gatherDtype =
      data_.obsspace().dtype(gatherVariable.group(), gatherVariable.variable());
  auto idDtype =
      data_.obsspace().dtype(idVariable.group(), idVariable.variable());

  // Verify all nearestNeighborIDVars have the same type as idVariable
  verifyNearestNeighborIDTypesMatch(idVariable, nearestNeighborIDVars);

  // Nested dispatch based on both types
  if (gatherDtype == ioda::ObsDtype::Float) {
    dispatchOnIDType<float>(params, options, idDtype);
  } else if (gatherDtype == ioda::ObsDtype::Integer) {
    dispatchOnIDType<int>(params, options, idDtype);
  } else if (gatherDtype == ioda::ObsDtype::Integer_64) {
    dispatchOnIDType<int64_t>(params, options, idDtype);
  } else if (gatherDtype == ioda::ObsDtype::String) {
    dispatchOnIDType<std::string>(params, options, idDtype);
  } else if (gatherDtype == ioda::ObsDtype::DateTime) {
    dispatchOnIDType<util::DateTime>(params, options, idDtype);
  } else {
    throw eckit::UserError("Unsupported type for gather variable.", Here());
  }
}

// Helper dispatcher for second type
template <typename GatherType>
void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::dispatchOnIDType(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &params,
    const UseNearestNeighborsParameters &options,
    ioda::ObsDtype idDtype) const {
  if (idDtype == ioda::ObsDtype::Integer) {
    executeImpl<GatherType, int>(params, options);
  } else if (idDtype == ioda::ObsDtype::Integer_64) {
    executeImpl<GatherType, int64_t>(params, options);
  } else if (idDtype == ioda::ObsDtype::String) {
    executeImpl<GatherType, std::string>(params, options);
  } else if (idDtype == ioda::ObsDtype::DateTime) {
    executeImpl<GatherType, util::DateTime>(params, options);
  } else {
    throw eckit::UserError("Unsupported type for identifier variable.", Here());
  }
}

// Actual implementation with both types templated
template <typename GatherType, typename IDType>
void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &params,
    const UseNearestNeighborsParameters &options) const {
  oops::Log::info()
      << "Executing use nearest neighbors gather and match timestamp algorithm"
      << std::endl;

  // Get and verify parameters:
  const Variable idVariable = options.idVar.value();
  const std::vector<Variable> &nearestNeighborIDVars =
      options.nearestNeighborIDVars.value();
  const std::vector<Variable> &outVars = params.outputVariables.value();
  const Variable gatherVariable = params.gatherVariable.value();
  const Variable timestampMatchVariable = params.timestampMatchVariable.value();
  for (const Variable &var : outVars) {
    if (var.group() == "ObsValue" || var.group() == "DerivedObsValue") {
      throw eckit::BadParameter(
          "Output variables cannot be in group ObsValue or DerivedObsValue, "
          "since qc flag information cannot be updated. "
          "Variable " +
              var.variable() + " is in group " + var.group(),
          Here());
    }
  }
  const auto timestampMatchDtype = data_.obsspace().dtype(
      timestampMatchVariable.group(), timestampMatchVariable.variable());
  if (timestampMatchDtype != ioda::ObsDtype::DateTime) {
    throw eckit::BadParameter(
        "Timestamp match variable must be of DateTime type", Here());
  }

  // get variables (rank local)
  std::vector<IDType> idValues;
  data_.obsspace().get_db(idVariable.group(), idVariable.variable(), idValues);
  std::vector<GatherType> gatherValues;
  data_.obsspace().get_db(gatherVariable.group(), gatherVariable.variable(),
                          gatherValues);
  std::vector<util::DateTime> timestampMatchValues;
  data_.obsspace().get_db(timestampMatchVariable.group(),
                          timestampMatchVariable.variable(),
                          timestampMatchValues);
  std::vector<util::DateTime> timestampValues;
  data_.obsspace().get_db("MetaData", "dateTime", timestampValues);

  // Algorithm itself

  // For each of "gatherVariable", "timestampMatchVariable", "idVariable"
  // and any "nearestNeighborIDVars" find the indices of values which are
  // (a) owned by this rank, (b) are not missing and (c) are not masked out
  // by the apply vector. "timestampMatchVariable" additionally needs to be
  // compared to the rank-local timestamp variable - where these do not
  // match, indices should be excluded. The intersection of these indices gives
  // the rank-local indices of the "gatherVariable", "timestampMatchVariable"
  // and "idVariable" that we want to gather and match against. These are
  // gathered from all ranks to give global vectors of possible matches.
  //
  // Next, for each of the "nearestNeighborIDVars", the indices of values which
  // are (a) not missing and (b) not masked out by the apply vector are found.
  // For each of these indices, the nearest neighbor ID value and rank-local
  // timestamp variable are compared against the global vectors of possible
  // "timestampMatchVariable" and "idVariable" values. Where both match, the
  // associated "gatherVariable" value is assigned to the output variable.

  const size_t nLocsThisRank = data_.obsspace().nlocs();

  std::vector<bool> isOwnedByThisRank(nLocsThisRank);
  obsdb_.distribution()->patchObs(isOwnedByThisRank);

  const std::vector<size_t> validIDIndices =
      orderedValidIndices(idValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validGatherIndices =
      orderedValidIndices(gatherValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validTimestampMatchIndices =
      orderedValidIndices(timestampMatchValues, isOwnedByThisRank, apply_);

  // additional filtering of timestampMatchIndices to only those that match
  // the rank-local timestamp variable
  std::vector<size_t> matchedValidTimestampMatchIndices;
  for (size_t i : validTimestampMatchIndices) {
    if (timestampMatchValues[i] == timestampValues[i]) {
      matchedValidTimestampMatchIndices.push_back(i);
    }
  }

  // The intersection of indices gives indices of 3 smaller vectors which we
  // will gather and match against
  std::vector<size_t> matchedValidIndices = findSortedIntersection(
      validIDIndices, validGatherIndices, matchedValidTimestampMatchIndices);

  // Gather matched valid global vectors of idValues, gatherValues and
  // timestampMatchValues
  std::vector<IDType> matchedIDValues;
  std::vector<GatherType> matchedGatherValues;
  std::vector<util::DateTime> matchedTimestampMatchValues;
  matchedIDValues.reserve(matchedValidIndices.size());
  matchedGatherValues.reserve(matchedValidIndices.size());
  matchedTimestampMatchValues.reserve(matchedValidIndices.size());
  for (size_t idx : matchedValidIndices) {
    matchedIDValues.push_back(idValues[idx]);
    matchedGatherValues.push_back(gatherValues[idx]);
    matchedTimestampMatchValues.push_back(timestampMatchValues[idx]);
  }
  // Allgatherv gets global vectors - will modify in place.
  allGatherv(matchedIDValues);
  allGatherv(matchedGatherValues);
  allGatherv(matchedTimestampMatchValues);

  // Build lookup map: (ID, timestamp) -> gather value
  std::unordered_map<std::pair<IDType, util::DateTime>, GatherType,
                     PairHash<IDType, util::DateTime>>
      lookupMap;
  for (size_t j = 0; j < matchedIDValues.size(); ++j) {
    const auto key =
        std::make_pair(matchedIDValues[j], matchedTimestampMatchValues[j]);
    const auto[it, inserted] =
        lookupMap.try_emplace(key, matchedGatherValues[j]);
    if (!inserted) {
      if (it->second != matchedGatherValues[j]) {
        oops::Log::warning()
            << "Duplicate (identifier, timestamp) key found with conflicting "
            << "gather values while building lookup map for gather and match "
            << "timestamp algorithm. key=(" << key.first << ", " << key.second
            << "), existing=" << it->second
            << ", incoming=" << matchedGatherValues[j]
            << ". Keeping incoming value." << std::endl;
      }
      it->second = matchedGatherValues[j];
    }
  }

  // For each nearest neighbor ID variable, find matches and assign values
  for (size_t ivar = 0; ivar < outVars.size(); ++ivar) {
    const Variable &nearestNeighborIDVar = nearestNeighborIDVars[ivar];
    const Variable &outVar = outVars[ivar];

    oops::Log::info() << "Processing output variable: " << outVar.fullName()
                      << " using nearest neighbor ID variable: "
                      << nearestNeighborIDVar.fullName() << std::endl;

    std::vector<IDType> nearestNeighborIDValues;
    data_.obsspace().get_db(nearestNeighborIDVar.group(),
                            nearestNeighborIDVar.variable(),
                            nearestNeighborIDValues);
    // Mask out missing/not-used indices, but do NOT restrict to
    // owned-by-this-rank only, since we need to populate the output everywhere
    const std::vector<bool> allOnRankMask(nLocsThisRank, true);
    const std::vector<size_t> validNearestNeighborIDIndices =
        orderedValidIndices(nearestNeighborIDValues, allOnRankMask, apply_);

    // For each of these indices, find the matching "gather values" where (A)
    // the nearest neighbor ID matches the ID value, and (B) the timestampMatch
    // matches the timestamp variable.
    std::vector<GatherType> outValues(nLocsThisRank,
                                      util::missingValue<GatherType>());
    for (size_t iloc : validNearestNeighborIDIndices) {
      auto key = std::make_pair(nearestNeighborIDValues[iloc],
                                timestampMatchValues[iloc]);
      auto it = lookupMap.find(key);
      if (it != lookupMap.end()) {
        outValues[iloc] = it->second;
      }
    }
    // Put the output variable back into the ObsSpace
    data_.obsspace().put_db(outVar.group(), outVar.variable(), outValues);
  }
  oops::Log::info()
      << "Use nearest neighbors gather and match timestamp algorithm complete"
      << std::endl;
}

// Explicit template instantiations for commonly used type combinations

// GatherType = float
// (Most common: GatherType = float, IDType = std::string)
template void
UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::dispatchOnIDType<float>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    float, std::string>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    float, int>(const UseNearestNeighborsGatherAndMatchTimestampParameters &,
                const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    float, int64_t>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    float, util::DateTime>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;

// GatherType = int
template void
UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::dispatchOnIDType<int>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int, int>(const UseNearestNeighborsGatherAndMatchTimestampParameters &,
              const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int, int64_t>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int, std::string>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int, util::DateTime>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;

// GatherType = std::string
template void
UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::dispatchOnIDType<
    std::string>(const UseNearestNeighborsGatherAndMatchTimestampParameters &,
                 const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    std::string, int>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    std::string, int64_t>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    std::string, std::string>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    std::string, util::DateTime>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;

// GatherType = int64_t
template void
UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::dispatchOnIDType<int64_t>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int64_t, int>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int64_t, int64_t>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int64_t, std::string>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    int64_t, util::DateTime>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;

// GatherType = util::DateTime
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::
    dispatchOnIDType<util::DateTime>(
        const UseNearestNeighborsGatherAndMatchTimestampParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    util::DateTime, int>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    util::DateTime, int64_t>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    util::DateTime, std::string>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmGatherAndMatchTimestamp::executeImpl<
    util::DateTime, util::DateTime>(
    const UseNearestNeighborsGatherAndMatchTimestampParameters &,
    const UseNearestNeighborsParameters &) const;

static UseNearestNeighborsAlgorithmMaker<
    UseNearestNeighborsAlgorithmGatherAndMatchTimestamp>
    makerGatherAndMatchTimestamp("gather and match timestamp");

}  // namespace ufo
