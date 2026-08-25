/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmReferencePointVariablesMean.h"
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    UseNearestNeighborsAlgorithmReferencePointVariablesMean(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &params,
        const ObsFilterData &data, const std::vector<bool> &apply,
        const Variables &filtervars, const ioda::ObsDataVector<int> &flags,
        std::vector<std::vector<bool>> &flagged)
    : UseNearestNeighborsAlgorithmBase(params, data, apply, filtervars, flags,
                                       flagged) {}

void UseNearestNeighborsAlgorithmReferencePointVariablesMean::execute(
    const UseNearestNeighborsAlgorithmParametersBase &algParams,
    const UseNearestNeighborsParameters &options) const {
  const auto &params = dynamic_cast<
      const UseNearestNeighborsReferencePointVariablesMeanParameters &>(
      algParams);

  // Get and verify parameters:
  const Variable idVariable = options.idVar.value();
  const std::vector<Variable> &nearestNeighborIDVars =
      options.nearestNeighborIDVars.value();
  const Variable gatherVariable = params.gatherVariable.value();
  const Variable meanAboutVariable = params.meanAboutVariable.value();

  // Determine types
  const auto gatherDtype =
      data_.obsspace().dtype(gatherVariable.group(), gatherVariable.variable());
  const auto meanAboutDtype = data_.obsspace().dtype(
      meanAboutVariable.group(), meanAboutVariable.variable());
  const auto idDtype =
      data_.obsspace().dtype(idVariable.group(), idVariable.variable());

  // Verify mean about variable is an integer
  if (meanAboutDtype != ioda::ObsDtype::Integer &&
      meanAboutDtype != ioda::ObsDtype::Integer_64) {
    throw eckit::BadParameter("Mean about variable must be of integer type",
                              Here());
  }

  // Verify all nearestNeighborIDVars have the same type as idVariable
  verifyNearestNeighborIDTypesMatch(idVariable, nearestNeighborIDVars);

  // Dispatch based on gather type and ID type
  if (gatherDtype == ioda::ObsDtype::Float) {
    dispatchOnIDType<float>(params, options, idDtype, meanAboutDtype);
  } else if (gatherDtype == ioda::ObsDtype::Integer) {
    dispatchOnIDType<int>(params, options, idDtype, meanAboutDtype);
  } else if (gatherDtype == ioda::ObsDtype::Integer_64) {
    dispatchOnIDType<int64_t>(params, options, idDtype, meanAboutDtype);
  } else if (gatherDtype == ioda::ObsDtype::DateTime) {
    dispatchOnIDType<util::DateTime>(params, options, idDtype, meanAboutDtype);
  } else {
    throw eckit::BadParameter(
        "Unsupported gather variable data type for reference point variables "
        "mean",
        Here());
  }
}

// Helper dispatcher for ID type
template <typename GatherType>
void UseNearestNeighborsAlgorithmReferencePointVariablesMean::dispatchOnIDType(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &params,
    const UseNearestNeighborsParameters &options, ioda::ObsDtype idDtype,
    ioda::ObsDtype meanAboutDtype) const {
  if (idDtype == ioda::ObsDtype::Integer) {
    dispatchOnMeanAboutType<GatherType, int>(params, options, meanAboutDtype);
  } else if (idDtype == ioda::ObsDtype::Integer_64) {
    dispatchOnMeanAboutType<GatherType, int64_t>(params, options,
                                                 meanAboutDtype);
  } else if (idDtype == ioda::ObsDtype::String) {
    dispatchOnMeanAboutType<GatherType, std::string>(params, options,
                                                     meanAboutDtype);
  } else if (idDtype == ioda::ObsDtype::DateTime) {
    dispatchOnMeanAboutType<GatherType, util::DateTime>(params, options,
                                                        meanAboutDtype);
  } else {
    throw eckit::BadParameter(
        "Unsupported ID variable data type for reference point variables mean",
        Here());
  }
}

// Helper dispatcher for mean-about type
template <typename GatherType, typename IDType>
void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &params,
        const UseNearestNeighborsParameters &options,
        ioda::ObsDtype meanAboutDtype) const {
  if (meanAboutDtype == ioda::ObsDtype::Integer) {
    executeImpl<GatherType, IDType, int>(params, options);
  } else if (meanAboutDtype == ioda::ObsDtype::Integer_64) {
    executeImpl<GatherType, IDType, int64_t>(params, options);
  } else {
    throw eckit::BadParameter(
        "Unsupported mean about variable data type for reference point "
        "variables mean",
        Here());
  }
}

// Actual implementation with GatherType and IDType templated
template <typename GatherType, typename IDType, typename MeanAboutType>
void UseNearestNeighborsAlgorithmReferencePointVariablesMean::executeImpl(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &params,
    const UseNearestNeighborsParameters &options) const {
  oops::Log::info()
      << "Executing use nearest neighbors reference point variables mean "
         "algorithm"
      << std::endl;

  // Algorithm itself

  // First we get rid of locations which are never going to be included in the
  // calculation. These are locations in the local obsvector which (a) are
  // masked out by `apply`, (b) as missing and (c) are not owned by this rank.
  // We need to do this for the gather variable, the mean-about variable,
  // and the identifier variable.
  //
  // The intersection of the 3 sets of obsvector indices gives the rank-local
  // indices of the observations which could be included in the calculation.
  // These are gathered from all ranks a map of
  // (identifier variable value, mean-about value) -> (gather variable value)
  // can be used for the calculation of the mean. Note that there is expected
  // redundancy in the map - once a particular (ID, mean-about) pair is found,
  // all corresponding gather variable values are assumed to be identical and
  // are not included in the map.
  //
  // Next we apply filtering to each of the nearest neighbor ID
  // variables to remove missing values and those masked out by `apply`. These
  // are used to look up matching (gather variable value, mean-about value)
  // pairs in the map built previously.
  //
  // The aim is to compute, for however many nearest-neighbor ID variables
  // are specified, the mean of the gather variable values corresponding to
  // the mean-about variable values. The first output variable corresponds to
  // the value where the mean-about variable is 0, the second where it is 1, and
  // so on.

  // Get parameters
  const Variable idVariable = options.idVar.value();
  const std::vector<Variable> &nearestNeighborIDVars =
      options.nearestNeighborIDVars.value();
  const std::vector<Variable> &outVars = params.outputVariables.value();
  const Variable gatherVariable = params.gatherVariable.value();
  const Variable meanAboutVariable = params.meanAboutVariable.value();

  for (const Variable &var : params.outputVariables.value()) {
    if (var.group() == "ObsValue" || var.group() == "DerivedObsValue") {
      throw eckit::BadParameter(
          "Output variables cannot be in group ObsValue or DerivedObsValue, "
          "since qc flag information cannot be updated. "
          "Variable " +
              var.variable() + " is in group " + var.group(),
          Here());
    }
  }

  // Get variables (rank local)
  std::vector<IDType> idValues;
  data_.obsspace().get_db(idVariable.group(), idVariable.variable(), idValues);
  std::vector<GatherType> gatherValues;
  data_.obsspace().get_db(gatherVariable.group(), gatherVariable.variable(),
                          gatherValues);
  std::vector<MeanAboutType> meanAboutValues;
  data_.obsspace().get_db(meanAboutVariable.group(),
                          meanAboutVariable.variable(), meanAboutValues);

  const size_t nLocsThisRank = data_.obsspace().nlocs();

  std::vector<bool> isOwnedByThisRank(nLocsThisRank);
  obsdb_.distribution()->patchObs(isOwnedByThisRank);

  // Get valid indices where ID, gather, and meanAbout are not missing
  // For building the lookup map, only use observations owned by this rank
  const std::vector<size_t> validIDIndices =
      orderedValidIndices(idValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validGatherIndices =
      orderedValidIndices(gatherValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validMeanAboutIndices =
      orderedValidIndices(meanAboutValues, isOwnedByThisRank, apply_);

  // Find intersection of all three valid index sets
  const std::vector<size_t> matchedSortedIndices = findSortedIntersection(
      validIDIndices, validGatherIndices, validMeanAboutIndices);

  // Extract valid data only
  std::vector<IDType> matchedIDValues;
  std::vector<GatherType> matchedGatherValues;
  std::vector<MeanAboutType> matchedMeanAboutValues;
  matchedIDValues.reserve(matchedSortedIndices.size());
  matchedGatherValues.reserve(matchedSortedIndices.size());
  matchedMeanAboutValues.reserve(matchedSortedIndices.size());

  for (size_t idx : matchedSortedIndices) {
    matchedIDValues.push_back(idValues[idx]);
    matchedGatherValues.push_back(gatherValues[idx]);
    matchedMeanAboutValues.push_back(meanAboutValues[idx]);
  }

  // Gather data from all ranks
  allGatherv(matchedIDValues);
  allGatherv(matchedGatherValues);
  allGatherv(matchedMeanAboutValues);

  // Build lookup map: (ID, mean-about) -> gather value
  const auto lookupMap = buildLookupMap<GatherType, IDType, MeanAboutType>(
      matchedIDValues, matchedGatherValues, matchedMeanAboutValues);

  // Get all nearest neighbor ID variables
  std::vector<std::vector<IDType>> allNearestNeighborIDValues(
      nearestNeighborIDVars.size());
  for (size_t jnn = 0; jnn < nearestNeighborIDVars.size(); ++jnn) {
    data_.obsspace().get_db(nearestNeighborIDVars[jnn].group(),
                            nearestNeighborIDVars[jnn].variable(),
                            allNearestNeighborIDValues[jnn]);
  }

  // Initialize output arrays for all variables
  std::vector<std::vector<GatherType>> allOutValues(outVars.size());
  for (size_t iOutVar = 0; iOutVar < outVars.size(); ++iOutVar) {
    allOutValues[iOutVar].resize(nLocsThisRank,
                                 util::missingValue<GatherType>());
  }

  // For each location in the obsvector
  for (size_t iloc = 0; iloc < nLocsThisRank; ++iloc) {
    // Skip if masked out by apply
    if (!apply_[iloc]) {
      continue;
    }

    // For each output variable (one per mean-about value: 0, 1, 2, ...)
    for (size_t iOutVar = 0; iOutVar < outVars.size(); ++iOutVar) {
      const MeanAboutType targetMeanAboutValue =
          static_cast<MeanAboutType>(iOutVar);

      // Collect values from all nearest neighbors for averaging at this
      // mean-about value
      const std::vector<GatherType> valuesToAverage =
          collectValuesForAveraging<GatherType, IDType, MeanAboutType>(
              iloc, targetMeanAboutValue, allNearestNeighborIDValues,
              lookupMap);

      // Compute mean if we have any values to average
      if (!valuesToAverage.empty()) {
        allOutValues[iOutVar][iloc] = computeMean<GatherType>(valuesToAverage);
      }
    }
  }

  // Save all output values to the ObsSpace
  for (size_t iOutVar = 0; iOutVar < outVars.size(); ++iOutVar) {
    const Variable &outVar = outVars[iOutVar];
    const std::string outputGroupName =
        outVar.group().empty() ? "DerivedMetaData" : outVar.group();
    const std::string outputVariableName = outVar.variable();
    data_.obsspace().put_db(outputGroupName, outputVariableName,
                            allOutValues[iOutVar]);
  }
}

// Helper method implementations

template <typename GatherType, typename IDType, typename MeanAboutType>
std::map<std::pair<IDType, MeanAboutType>, GatherType>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap(
    const std::vector<IDType> &matchedIDValues,
    const std::vector<GatherType> &matchedGatherValues,
    const std::vector<MeanAboutType> &matchedMeanAboutValues) const {
  // Build lookup map: (ID, mean-about) -> gather value
  // Note: Due to expected redundancy, we only store the first occurrence of
  // each (ID, mean-about) pair
  std::map<std::pair<IDType, MeanAboutType>, GatherType> lookupMap;
  for (size_t j = 0; j < matchedIDValues.size(); ++j) {
    // Only insert if not already present (keep first value, ignore redundant
    // ones)
    auto key = std::make_pair(matchedIDValues[j], matchedMeanAboutValues[j]);
    if (lookupMap.find(key) == lookupMap.end()) {
      lookupMap[key] = matchedGatherValues[j];
    }
  }
  return lookupMap;
}

template <typename GatherType>
GatherType UseNearestNeighborsAlgorithmReferencePointVariablesMean::computeMean(
    const std::vector<GatherType> &valuesToAverage) const {
  GatherType sum = std::accumulate(valuesToAverage.begin(),
                                   valuesToAverage.end(), GatherType{});
  return sum / static_cast<GatherType>(valuesToAverage.size());
}

// Specialization required for util::DateTime
template <>
util::DateTime
UseNearestNeighborsAlgorithmReferencePointVariablesMean::computeMean<
    util::DateTime>(const std::vector<util::DateTime> &valuesToAverage) const {
  // Use first datetime as reference epoch
  const util::DateTime &reference = valuesToAverage[0];
  int64_t totalOffsetSeconds = 0;
  for (size_t i = 1; i < valuesToAverage.size(); ++i) {
    totalOffsetSeconds += (valuesToAverage[i] - reference).toSeconds();
  }
  int64_t meanOffsetSeconds =
      totalOffsetSeconds / static_cast<int64_t>(valuesToAverage.size());
  return reference + util::Duration(meanOffsetSeconds);
}

template <typename GatherType, typename IDType, typename MeanAboutType>
std::vector<GatherType>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging(
        size_t iloc, MeanAboutType targetMeanAboutValue,
        const std::vector<std::vector<IDType>> &allNearestNeighborIDValues,
        const std::map<std::pair<IDType, MeanAboutType>, GatherType> &lookupMap)
        const {
  std::vector<GatherType> valuesToAverage;
  valuesToAverage.reserve(allNearestNeighborIDValues.size());

  // For each nearest neighbor ID variable, look up the value
  // for this location's nearest neighbor station at the target mean-about value
  for (size_t jnn = 0; jnn < allNearestNeighborIDValues.size(); ++jnn) {
    const IDType &nearestNeighborID = allNearestNeighborIDValues[jnn][iloc];

    // Skip if this nearest neighbor ID is missing
    if (nearestNeighborID == util::missingValue<IDType>()) {
      continue;
    }

    // Look up (nearestNeighborID, targetMeanAboutValue) in the map
    auto key = std::make_pair(nearestNeighborID, targetMeanAboutValue);
    auto it = lookupMap.find(key);

    if (it != lookupMap.end()) {
      // Found matching (ID, mean-about) pair - collect value for averaging
      valuesToAverage.push_back(it->second);
    }
  }

  return valuesToAverage;
}

// Explicit template instantiations for commonly used type combinations

// GatherType = float
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::dispatchOnIDType<
    float>(const UseNearestNeighborsReferencePointVariablesMeanParameters &,
           const UseNearestNeighborsParameters &, ioda::ObsDtype,
           ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<float, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<float, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<float, std::string>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<float, util::DateTime>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::executeImpl<float, int,
                                                                     int>(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &,
    const UseNearestNeighborsParameters &) const;
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::executeImpl<float, int,
                                                                     int64_t>(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<float, int64_t, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<float, int64_t, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<float, std::string, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<float, std::string, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<float, util::DateTime, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<float, util::DateTime, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template std::map<std::pair<int, int>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, int, int>(const std::vector<int> &, const std::vector<float> &,
                     const std::vector<int> &) const;
template std::map<std::pair<int, int64_t>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, int, int64_t>(const std::vector<int> &, const std::vector<float> &,
                         const std::vector<int64_t> &) const;
template std::map<std::pair<int64_t, int>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, int64_t, int>(const std::vector<int64_t> &,
                         const std::vector<float> &,
                         const std::vector<int> &) const;
template std::map<std::pair<int64_t, int64_t>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, int64_t, int64_t>(const std::vector<int64_t> &,
                             const std::vector<float> &,
                             const std::vector<int64_t> &) const;
template std::map<std::pair<std::string, int>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, std::string, int>(const std::vector<std::string> &,
                             const std::vector<float> &,
                             const std::vector<int> &) const;
template std::map<std::pair<std::string, int64_t>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, std::string, int64_t>(const std::vector<std::string> &,
                                 const std::vector<float> &,
                                 const std::vector<int64_t> &) const;
template std::map<std::pair<util::DateTime, int>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, util::DateTime, int>(const std::vector<util::DateTime> &,
                                const std::vector<float> &,
                                const std::vector<int> &) const;
template std::map<std::pair<util::DateTime, int64_t>, float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    float, util::DateTime, int64_t>(const std::vector<util::DateTime> &,
                                    const std::vector<float> &,
                                    const std::vector<int64_t> &) const;
template float
UseNearestNeighborsAlgorithmReferencePointVariablesMean::computeMean<float>(
    const std::vector<float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, int, int>(
        size_t, int, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, int, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int64_t>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, int64_t, int>(
        size_t, int, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, int64_t, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int64_t>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, std::string, int>(
        size_t, int, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, std::string, int64_t>(
        size_t, int64_t, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int64_t>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, util::DateTime, int>(
        size_t, int, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int>, float> &) const;
template std::vector<float>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<float, util::DateTime, int64_t>(
        size_t, int64_t, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int64_t>, float> &) const;

// GatherType = int
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::dispatchOnIDType<int>(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype,
    ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int, std::string>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int, util::DateTime>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::executeImpl<int, int,
                                                                     int>(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &,
    const UseNearestNeighborsParameters &) const;
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::executeImpl<int, int,
                                                                     int64_t>(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int, int64_t, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int, int64_t, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int, std::string, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int, std::string, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int, util::DateTime, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int, util::DateTime, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template std::map<std::pair<int, int>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, int, int>(const std::vector<int> &, const std::vector<int> &,
                   const std::vector<int> &) const;
template std::map<std::pair<int, int64_t>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, int, int64_t>(const std::vector<int> &, const std::vector<int> &,
                       const std::vector<int64_t> &) const;
template std::map<std::pair<int64_t, int>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, int64_t, int>(const std::vector<int64_t> &, const std::vector<int> &,
                       const std::vector<int> &) const;
template std::map<std::pair<int64_t, int64_t>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, int64_t, int64_t>(const std::vector<int64_t> &,
                           const std::vector<int> &,
                           const std::vector<int64_t> &) const;
template std::map<std::pair<std::string, int>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, std::string, int>(const std::vector<std::string> &,
                           const std::vector<int> &,
                           const std::vector<int> &) const;
template std::map<std::pair<std::string, int64_t>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, std::string, int64_t>(const std::vector<std::string> &,
                               const std::vector<int> &,
                               const std::vector<int64_t> &) const;
template std::map<std::pair<util::DateTime, int>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, util::DateTime, int>(const std::vector<util::DateTime> &,
                              const std::vector<int> &,
                              const std::vector<int> &) const;
template std::map<std::pair<util::DateTime, int64_t>, int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int, util::DateTime, int64_t>(const std::vector<util::DateTime> &,
                                  const std::vector<int> &,
                                  const std::vector<int64_t> &) const;
template int
UseNearestNeighborsAlgorithmReferencePointVariablesMean::computeMean<int>(
    const std::vector<int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, int, int>(
        size_t, int, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, int, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int64_t>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, int64_t, int>(
        size_t, int, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, int64_t, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int64_t>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, std::string, int>(
        size_t, int, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, std::string, int64_t>(
        size_t, int64_t, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int64_t>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, util::DateTime, int>(
        size_t, int, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int>, int> &) const;
template std::vector<int>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int, util::DateTime, int64_t>(
        size_t, int64_t, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int64_t>, int> &) const;

// GatherType = int64_t
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::dispatchOnIDType<
    int64_t>(const UseNearestNeighborsReferencePointVariablesMeanParameters &,
             const UseNearestNeighborsParameters &, ioda::ObsDtype,
             ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int64_t, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int64_t, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int64_t, std::string>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<int64_t, util::DateTime>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void
UseNearestNeighborsAlgorithmReferencePointVariablesMean::executeImpl<int64_t,
                                                                     int, int>(
    const UseNearestNeighborsReferencePointVariablesMeanParameters &,
    const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, int, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, int64_t, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, int64_t, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, std::string, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, std::string, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, util::DateTime, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<int64_t, util::DateTime, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template std::map<std::pair<int, int>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, int, int>(const std::vector<int> &, const std::vector<int64_t> &,
                       const std::vector<int> &) const;
template std::map<std::pair<int, int64_t>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, int, int64_t>(const std::vector<int> &,
                           const std::vector<int64_t> &,
                           const std::vector<int64_t> &) const;
template std::map<std::pair<int64_t, int>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, int64_t, int>(const std::vector<int64_t> &,
                           const std::vector<int64_t> &,
                           const std::vector<int> &) const;
template std::map<std::pair<int64_t, int64_t>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, int64_t, int64_t>(const std::vector<int64_t> &,
                               const std::vector<int64_t> &,
                               const std::vector<int64_t> &) const;
template std::map<std::pair<std::string, int>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, std::string, int>(const std::vector<std::string> &,
                               const std::vector<int64_t> &,
                               const std::vector<int> &) const;
template std::map<std::pair<std::string, int64_t>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, std::string, int64_t>(const std::vector<std::string> &,
                                   const std::vector<int64_t> &,
                                   const std::vector<int64_t> &) const;
template std::map<std::pair<util::DateTime, int>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, util::DateTime, int>(const std::vector<util::DateTime> &,
                                  const std::vector<int64_t> &,
                                  const std::vector<int> &) const;
template std::map<std::pair<util::DateTime, int64_t>, int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    int64_t, util::DateTime, int64_t>(const std::vector<util::DateTime> &,
                                      const std::vector<int64_t> &,
                                      const std::vector<int64_t> &) const;
template int64_t
UseNearestNeighborsAlgorithmReferencePointVariablesMean::computeMean<int64_t>(
    const std::vector<int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, int, int>(
        size_t, int, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, int, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int64_t>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, int64_t, int>(
        size_t, int, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, int64_t, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int64_t>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, std::string, int>(
        size_t, int, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, std::string, int64_t>(
        size_t, int64_t, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int64_t>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, util::DateTime, int>(
        size_t, int, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int>, int64_t> &) const;
template std::vector<int64_t>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<int64_t, util::DateTime, int64_t>(
        size_t, int64_t, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int64_t>, int64_t> &) const;

// GatherType = util::DateTime
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnIDType<util::DateTime>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype,
        ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<util::DateTime, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<util::DateTime, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<util::DateTime, std::string>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    dispatchOnMeanAboutType<util::DateTime, util::DateTime>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, int, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, int, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, int64_t, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, int64_t, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, std::string, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, std::string, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, util::DateTime, int>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template void UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    executeImpl<util::DateTime, util::DateTime, int64_t>(
        const UseNearestNeighborsReferencePointVariablesMeanParameters &,
        const UseNearestNeighborsParameters &) const;
template std::map<std::pair<int, int>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, int, int>(const std::vector<int> &,
                              const std::vector<util::DateTime> &,
                              const std::vector<int> &) const;
template std::map<std::pair<int, int64_t>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, int, int64_t>(const std::vector<int> &,
                                  const std::vector<util::DateTime> &,
                                  const std::vector<int64_t> &) const;
template std::map<std::pair<int64_t, int>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, int64_t, int>(const std::vector<int64_t> &,
                                  const std::vector<util::DateTime> &,
                                  const std::vector<int> &) const;
template std::map<std::pair<int64_t, int64_t>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, int64_t, int64_t>(const std::vector<int64_t> &,
                                      const std::vector<util::DateTime> &,
                                      const std::vector<int64_t> &) const;
template std::map<std::pair<std::string, int>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, std::string, int>(const std::vector<std::string> &,
                                      const std::vector<util::DateTime> &,
                                      const std::vector<int> &) const;
template std::map<std::pair<std::string, int64_t>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, std::string, int64_t>(const std::vector<std::string> &,
                                          const std::vector<util::DateTime> &,
                                          const std::vector<int64_t> &) const;
template std::map<std::pair<util::DateTime, int>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, util::DateTime, int>(const std::vector<util::DateTime> &,
                                         const std::vector<util::DateTime> &,
                                         const std::vector<int> &) const;
template std::map<std::pair<util::DateTime, int64_t>, util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::buildLookupMap<
    util::DateTime, util::DateTime, int64_t>(
    const std::vector<util::DateTime> &, const std::vector<util::DateTime> &,
    const std::vector<int64_t> &) const;
// Note: computeMean<util::DateTime> is a template specialization, not an
// instantiation
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, int, int>(
        size_t, int, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int>, util::DateTime> &) const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, int, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int>> &,
        const std::map<std::pair<int, int64_t>, util::DateTime> &) const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, int64_t, int>(
        size_t, int, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int>, util::DateTime> &) const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, int64_t, int64_t>(
        size_t, int64_t, const std::vector<std::vector<int64_t>> &,
        const std::map<std::pair<int64_t, int64_t>, util::DateTime> &) const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, std::string, int>(
        size_t, int, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int>, util::DateTime> &) const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, std::string, int64_t>(
        size_t, int64_t, const std::vector<std::vector<std::string>> &,
        const std::map<std::pair<std::string, int64_t>, util::DateTime> &)
        const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, util::DateTime, int>(
        size_t, int, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int>, util::DateTime> &) const;
template std::vector<util::DateTime>
UseNearestNeighborsAlgorithmReferencePointVariablesMean::
    collectValuesForAveraging<util::DateTime, util::DateTime, int64_t>(
        size_t, int64_t, const std::vector<std::vector<util::DateTime>> &,
        const std::map<std::pair<util::DateTime, int64_t>, util::DateTime> &)
        const;

static UseNearestNeighborsAlgorithmMaker<
    UseNearestNeighborsAlgorithmReferencePointVariablesMean>
    makerReferencePointVariablesMean("reference point variables mean");

}  // namespace ufo
