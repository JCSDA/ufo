/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmLocalPlaneFit.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmBase.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// Type aliases for brevity
using LocalPlaneFit = UseNearestNeighborsAlgorithmLocalPlaneFit;

template <typename GatherType>
using NeighborData = typename LocalPlaneFit::template NeighborData<GatherType>;

template <typename GatherType>
using PlaneFitResult =
    typename LocalPlaneFit::template PlaneFitResult<GatherType>;

// Helper: compute inverse distance weighted value from pre-computed
// weights
template <typename T>
static T computeInverseDistanceWeighting(
    const std::vector<T> &values, const std::vector<float> &normalizedWeights) {
  const size_t k = values.size();
  double result =
      0.0;  // accumulate in double regardless of T to minimize rounding errors
  for (size_t i = 0; i < k; ++i) {
    result += static_cast<double>(values[i]) *
              static_cast<double>(normalizedWeights[i]);
  }
  return static_cast<T>(result);
}

LocalPlaneFit::UseNearestNeighborsAlgorithmLocalPlaneFit(
    const UseNearestNeighborsLocalPlaneFitParameters &params,
    const ObsFilterData &data, const std::vector<bool> &apply,
    const Variables &filtervars, const ioda::ObsDataVector<int> &flags,
    std::vector<std::vector<bool>> &flagged)
    : UseNearestNeighborsAlgorithmBase(params, data, apply, filtervars, flags,
                                       flagged) {}

void LocalPlaneFit::execute(
    const UseNearestNeighborsAlgorithmParametersBase &algParams,
    const UseNearestNeighborsParameters &options) const {
  const auto &params =
      dynamic_cast<const UseNearestNeighborsLocalPlaneFitParameters &>(
          algParams);

  // Get and verify parameters
  const Variable idVariable = options.idVar.value();
  const std::vector<Variable> &nearestNeighborIDVars =
      options.nearestNeighborIDVars.value();
  const std::vector<Variable> &distanceVars = params.distanceVariables.value();

  // Validate that number of distance variables matches number of nearest
  // neighbor ID variables
  if (distanceVars.size() != nearestNeighborIDVars.size()) {
    throw eckit::BadParameter(
        "Number of distance variables (" + std::to_string(distanceVars.size()) +
            ") must match number of nearest neighbor identifier variables (" +
            std::to_string(nearestNeighborIDVars.size()) + ").",
        Here());
  }

  // Verify all nearestNeighborIDVars have the same type as idVariable
  verifyNearestNeighborIDTypesMatch(idVariable, nearestNeighborIDVars);

  // Determine gather type (reference/query point variable type)
  const Variable referencePointVariable = params.referencePointVariable.value();
  auto gatherDtype = data_.obsspace().dtype(referencePointVariable.group(),
                                            referencePointVariable.variable());
  auto idDtype =
      data_.obsspace().dtype(idVariable.group(), idVariable.variable());

  // Dispatch based on gather type
  if (gatherDtype == ioda::ObsDtype::Float) {
    dispatchOnIDType<float>(params, options, idDtype);
  } else if (gatherDtype == ioda::ObsDtype::Integer) {
    dispatchOnIDType<int>(params, options, idDtype);
  } else if (gatherDtype == ioda::ObsDtype::Integer_64) {
    dispatchOnIDType<int64_t>(params, options, idDtype);
  } else {
    throw eckit::UserError(
        "Local plane fit algorithm requires numeric (float, int, or int64) "
        "reference point variable type.",
        Here());
  }
}

// Helper: dispatcher for ID type
template <typename GatherType>
void LocalPlaneFit::dispatchOnIDType(
    const UseNearestNeighborsLocalPlaneFitParameters &params,
    const UseNearestNeighborsParameters &options,
    const ioda::ObsDtype idDtype) const {
  const ioda::ObsDtype matchDtype =
      data_.obsspace().dtype(params.matchVariable.value().group(),
                             params.matchVariable.value().variable());
  if (idDtype == ioda::ObsDtype::Integer) {
    dispatchOnMatchType<GatherType, int>(params, options, matchDtype);
  } else if (idDtype == ioda::ObsDtype::Integer_64) {
    dispatchOnMatchType<GatherType, int64_t>(params, options, matchDtype);
  } else if (idDtype == ioda::ObsDtype::String) {
    dispatchOnMatchType<GatherType, std::string>(params, options, matchDtype);
  } else if (idDtype == ioda::ObsDtype::DateTime) {
    dispatchOnMatchType<GatherType, util::DateTime>(params, options,
                                                    matchDtype);
  } else {
    throw eckit::UserError("Unsupported type for identifier variable.", Here());
  }
}

// Helper: dispatcher for match type
template <typename GatherType, typename IDType>
void LocalPlaneFit::dispatchOnMatchType(
    const UseNearestNeighborsLocalPlaneFitParameters &params,
    const UseNearestNeighborsParameters &options,
    const ioda::ObsDtype matchDtype) const {
  if (matchDtype == ioda::ObsDtype::Integer) {
    executeImpl<GatherType, IDType, int>(params, options);
  } else if (matchDtype == ioda::ObsDtype::Integer_64) {
    executeImpl<GatherType, IDType, int64_t>(params, options);
  } else if (matchDtype == ioda::ObsDtype::String) {
    executeImpl<GatherType, IDType, std::string>(params, options);
  } else if (matchDtype == ioda::ObsDtype::DateTime) {
    executeImpl<GatherType, IDType, util::DateTime>(params, options);
  } else {
    throw eckit::UserError("Unsupported type for match variable.", Here());
  }
}

// Helper: Gather reference data and build lookup map
template <typename GatherType, typename IDType, typename MatchType>
typename LocalPlaneFit::template LookupMap<GatherType, IDType, MatchType>
LocalPlaneFit::gatherReferenceData(
    const Variable &idVariable, const Variable &referencePointVariable,
    const std::vector<bool> &isOwnedByThisRank,
    const std::vector<IDType> &idValues,
    const std::vector<GatherType> &referencePointValues,
    const std::vector<float> &latValues, const std::vector<float> &lonValues,
    const std::vector<MatchType> &matchValues) const {
  // Find valid indices for reference data (owned, not missing, not masked)
  const std::vector<size_t> validIDIndices =
      orderedValidIndices(idValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validReferenceIndices =
      orderedValidIndices(referencePointValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validLatIndices =
      orderedValidIndices(latValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validLonIndices =
      orderedValidIndices(lonValues, isOwnedByThisRank, apply_);
  const std::vector<size_t> validMatchIndices =
      orderedValidIndices(matchValues, isOwnedByThisRank, apply_);

  // Find intersection of all valid indices
  std::vector<size_t> matchedValidIndices = findSortedIntersection(
      validIDIndices, validReferenceIndices, validLatIndices, validLonIndices,
      validMatchIndices);

  // Gather matched valid global vectors
  std::vector<IDType> matchedIDValues(matchedValidIndices.size());
  std::vector<GatherType> matchedReferenceValues(matchedValidIndices.size());
  std::vector<float> matchedLatValues(matchedValidIndices.size());
  std::vector<float> matchedLonValues(matchedValidIndices.size());
  std::vector<MatchType> matchedMatchValues(matchedValidIndices.size());

  for (size_t j = 0; j < matchedValidIndices.size(); ++j) {
    const size_t idx = matchedValidIndices[j];
    matchedIDValues[j] = idValues[idx];
    matchedReferenceValues[j] = referencePointValues[idx];
    matchedLatValues[j] = latValues[idx];
    matchedLonValues[j] = lonValues[idx];
    matchedMatchValues[j] = matchValues[idx];
  }

  // Gather from all ranks
  allGatherv(matchedIDValues);
  allGatherv(matchedReferenceValues);
  allGatherv(matchedLatValues);
  allGatherv(matchedLonValues);
  allGatherv(matchedMatchValues);

  // Build lookup map keyed on (stationID, matchValue)
  typename LocalPlaneFit::template LookupMap<GatherType, IDType, MatchType>
      lookupMap;
  for (size_t j = 0; j < matchedIDValues.size(); ++j) {
    auto key = std::make_pair(matchedIDValues[j], matchedMatchValues[j]);
    if (lookupMap.count(key) > 0) {
      throw eckit::BadValue(
          "The match variable is not unique for the reference data at the "
          "reference identifier.",
          Here());
    }
    lookupMap[key] = std::make_tuple(matchedReferenceValues[j],
                                     matchedLatValues[j], matchedLonValues[j]);
  }

  return lookupMap;
}

// Helper: Collect neighbor data for a single location
template <typename GatherType, typename IDType, typename MatchType>
NeighborData<GatherType> LocalPlaneFit::collectNeighborData(
    size_t iloc, size_t numNeighbors,
    const std::vector<std::vector<IDType>> &allNearestNeighborIDs,
    const std::vector<std::vector<float>> &allDistances,
    const typename LocalPlaneFit::template LookupMap<GatherType, IDType,
                                                     MatchType> &lookupMap,
    MatchType matchValue) const {
  NeighborData<GatherType> result;
  result.isValid = false;
  result.values.reserve(numNeighbors);
  result.latitudes.reserve(numNeighbors);
  result.longitudes.reserve(numNeighbors);
  result.distances.reserve(numNeighbors);

  for (size_t i = 0; i < numNeighbors; ++i) {
    // Check if neighbor ID, distance or match value is missing
    if (allNearestNeighborIDs[i][iloc] == util::missingValue<IDType>() ||
        allDistances[i][iloc] == util::missingValue<float>() ||
        matchValue == util::missingValue<MatchType>()) {
      return result;
    }

    auto it = lookupMap.find({allNearestNeighborIDs[i][iloc], matchValue});
    if (it == lookupMap.end()) {
      return result;
    }

    const auto &[refValue, refLat, refLon] = it->second;

    // Check if reference value is missing
    if (refValue == util::missingValue<GatherType>()) {
      return result;
    }

    result.values.push_back(refValue);
    result.latitudes.push_back(refLat);
    result.longitudes.push_back(refLon);
    result.distances.push_back(allDistances[i][iloc]);
  }

  result.isValid = true;
  return result;
}

// Helper: Compute plane fit, indicating where fallback to IDW is needed
template <typename GatherType>
PlaneFitResult<GatherType> LocalPlaneFit::computePlaneFit(
    size_t iloc, float lat_q, float lon_q, float relativeErrorThreshold,
    const std::vector<GatherType> &neighborValues,
    const std::vector<float> &neighborLats,
    const std::vector<float> &neighborLons,
    const std::vector<float> &weights) const {
  PlaneFitResult<GatherType> result;
  result.success = false;
  result.value = util::missingValue<GatherType>();
  result.reason = PlaneFitResult<GatherType>::FailureReason::None;

  const size_t k = neighborValues.size();
  const float lat_q_rad = lat_q * M_PI / 180.0f;
  const float cos_lat_q = std::cos(lat_q_rad);

  // Earth radius in km (for equirectangular approximation)
  const float R = ufo::Constants::mean_earth_rad_m / 1000.0f;

  // Compute local equirectangular coordinates centered at query point
  Eigen::VectorXf dx(k);
  Eigen::VectorXf dy(k);
  Eigen::VectorXf z(k);

  for (size_t i = 0; i < k; ++i) {
    float dlon = (neighborLons[i] - lon_q) * M_PI / 180.0f;
    float dlat = (neighborLats[i] - lat_q) * M_PI / 180.0f;
    dx(i) = dlon * cos_lat_q * R;
    dy(i) = dlat * R;
    z(i) = neighborValues[i];
  }

  // Build design matrix A = [dx, dy, 1]
  Eigen::MatrixXf A(k, 3);
  A.col(0) = dx;
  A.col(1) = dy;
  A.col(2) = Eigen::VectorXf::Ones(k);

  // Build diagonal weight matrix
  Eigen::VectorXf W(k);
  for (size_t i = 0; i < k; ++i) {
    W(i) = weights[i];
  }

  try {
    // Solve weighted least squares: (A^T W A) coeffs = A^T W z
    Eigen::MatrixXf AtW = A.transpose() * W.asDiagonal();
    Eigen::MatrixXf AtWA = AtW * A;
    Eigen::VectorXf AtWz = AtW * z;

    // Use LDLT decomposition (more robust than Cholesky for near-singular
    // matrices, e.g. collinear points where k=3)
    Eigen::LDLT<Eigen::MatrixXf> ldlt(AtWA);

    if (ldlt.info() != Eigen::Success || !ldlt.isPositive()) {
      // LDLT failed or matrix is not positive semidefinite, fall back to IDW
      oops::Log::debug()
          << "Local plane fit: LDLT decomposition failed at location " << iloc
          << " (k=" << k << "), falling back to IDW" << std::endl;
      result.reason = PlaneFitResult<GatherType>::FailureReason::LDLT;
      return result;
    }

    Eigen::VectorXf coeffs = ldlt.solve(AtWz);

    // Check relative error to ensure plane fit is reasonable
    Eigen::VectorXf z_predicted = A * coeffs;
    Eigen::VectorXf errors = z - z_predicted;

    // Compute weighted RMS error consistent with the least-squares system:
    // sqrt(sum_i w_i * e_i^2)
    const float weightedRmsError =
        std::sqrt((errors.array().square() * W.array()).sum());

    // Compute weighted RMS of neighbor values as scale:
    // sqrt(sum_i w_i * z_i^2)
    const float weightedRmsNeighbors =
        std::sqrt((z.array().square() * W.array()).sum());

    const float relativeError =
        weightedRmsError / (weightedRmsNeighbors + 1e-10f);

    // If relative error exceeds threshold, fall back to IDW
    if (relativeError > relativeErrorThreshold) {
      oops::Log::debug() << "Local plane fit: relative error " << relativeError
                         << " exceeds threshold " << relativeErrorThreshold
                         << " at location " << iloc << ", falling back to IDW"
                         << std::endl;
      result.reason = PlaneFitResult<GatherType>::FailureReason::RelativeError;
      return result;
    }

    // Query point is at origin of local coordinate system, so intercept is the
    // result
    result.success = true;
    result.value = static_cast<GatherType>(coeffs(2));
    return result;
  } catch (const std::exception &e) {
    // If any exception occurs, fall back to IDW
    oops::Log::debug() << "Local plane fit: exception caught at location "
                       << iloc << " (" << e.what() << "), falling back to IDW"
                       << std::endl;
    result.reason = PlaneFitResult<GatherType>::FailureReason::Exception;
    return result;
  }
}

// Actual implementation with GatherType, IDType, and MatchType templated
template <typename GatherType, typename IDType, typename MatchType>
void LocalPlaneFit::executeImpl(
    const UseNearestNeighborsLocalPlaneFitParameters &params,
    const UseNearestNeighborsParameters &options) const {
  oops::Log::info()
      << "Executing use nearest neighbors local plane fit algorithm"
      << std::endl;

  // Get parameters
  const Variable idVariable = options.idVar.value();
  const std::vector<Variable> &nearestNeighborIDVars =
      options.nearestNeighborIDVars.value();
  const Variable queryPointVariable = params.queryPointVariable.value();
  const Variable referencePointVariable = params.referencePointVariable.value();
  const std::vector<Variable> &distanceVars = params.distanceVariables.value();
  const Variable outputVariable = params.outputVariable.value();
  const float power = params.inverseDistanceWeightingPower.value();
  const float relativeErrorThreshold =
      params.relativeErrorThreshold.value().value_or(0.25f);

  if (outputVariable.group() == "ObsValue" ||
      outputVariable.group() == "DerivedObsValue") {
    throw eckit::BadParameter(
        "Output variables cannot be in group ObsValue or DerivedObsValue, "
        "since qc flag information cannot be updated. "
        "Variable " +
            outputVariable.variable() + " is in group " +
            outputVariable.group(),
        Here());
  }

  const size_t numNeighbors = nearestNeighborIDVars.size();

  // Get rank-local data
  const size_t nLocsThisRank = data_.obsspace().nlocs();

  std::vector<bool> isOwnedByThisRank(nLocsThisRank);
  obsdb_.distribution()->patchObs(isOwnedByThisRank);

  // Read ID variable
  std::vector<IDType> idValues;
  data_.obsspace().get_db(idVariable.group(), idVariable.variable(), idValues);

  // Read query point variable
  std::vector<GatherType> queryPointValues;
  data_.obsspace().get_db(queryPointVariable.group(),
                          queryPointVariable.variable(), queryPointValues);

  // Read reference point variable
  std::vector<GatherType> referencePointValues;
  data_.obsspace().get_db(referencePointVariable.group(),
                          referencePointVariable.variable(),
                          referencePointValues);

  // Read latitude and longitude (for equirectangular coordinates)
  std::vector<float> latValues;
  std::vector<float> lonValues;
  data_.obsspace().get_db("MetaData", "latitude", latValues);
  data_.obsspace().get_db("MetaData", "longitude", lonValues);

  // Read match variable
  std::vector<MatchType> matchValues(nLocsThisRank);
  const Variable &matchVar = params.matchVariable.value();
  data_.obsspace().get_db(matchVar.group(), matchVar.variable(), matchValues);

  // Gather reference data keyed on (stationID, matchValue)
  auto lookupMap = gatherReferenceData<GatherType, IDType, MatchType>(
      idVariable, referencePointVariable, isOwnedByThisRank, idValues,
      referencePointValues, latValues, lonValues, matchValues);

  // Read all nearest neighbor IDs and distances
  std::vector<std::vector<IDType>> allNearestNeighborIDs(numNeighbors);
  std::vector<std::vector<float>> allDistances(numNeighbors);
  for (size_t i = 0; i < numNeighbors; ++i) {
    data_.obsspace().get_db(nearestNeighborIDVars[i].group(),
                            nearestNeighborIDVars[i].variable(),
                            allNearestNeighborIDs[i]);
    data_.obsspace().get_db(distanceVars[i].group(), distanceVars[i].variable(),
                            allDistances[i]);
  }

  // Initialize output variable with missing values
  std::vector<GatherType> outValues(nLocsThisRank,
                                    util::missingValue<GatherType>());

  // Counters for fallback to IDW cases
  size_t fallbackLDLT = 0;
  size_t fallbackRelativeError = 0;
  size_t fallbackException = 0;

  // Counters for successful cases
  size_t colocatedPoints = 0;
  size_t planeFitSuccess = 0;
  size_t idwKLessThan3 = 0;

  // Process all locations (not just owned ones, for output everywhere)
  for (size_t iloc = 0; iloc < nLocsThisRank; ++iloc) {
    if (!apply_[iloc]) continue;

    // Skip if query point value is missing
    if (queryPointValues[iloc] == util::missingValue<GatherType>()) continue;

    // Skip if query location has missing lat/lon
    if (latValues[iloc] == util::missingValue<float>() ||
        lonValues[iloc] == util::missingValue<float>()) {
      continue;
    }

    const float lat_q = latValues[iloc];
    const float lon_q = lonValues[iloc];

    // Collect nearest neighbor data
    auto neighborData = collectNeighborData<GatherType, IDType, MatchType>(
        iloc, numNeighbors, allNearestNeighborIDs, allDistances, lookupMap,
        matchValues[iloc]);

    if (!neighborData.isValid) {
      continue;  // Skip if any data is invalid
    }

    const size_t k = neighborData.values.size();

    // If query point coincides with a neighbor (distance < 1e-10 km = 0.01 m)
    float minDist = *std::min_element(neighborData.distances.begin(),
                                      neighborData.distances.end());
    if (minDist < 1e-10f) {
      size_t minIdx = std::min_element(neighborData.distances.begin(),
                                       neighborData.distances.end()) -
                      neighborData.distances.begin();
      outValues[iloc] = neighborData.values[minIdx];
      ++colocatedPoints;
      continue;
    }

    // Compute inverse distance weights
    std::vector<float> weights(k);
    float weightSum = 0.0f;
    for (size_t i = 0; i < k; ++i) {
      weights[i] = 1.0f / (std::pow(neighborData.distances[i], power) + 1e-10f);
      weightSum += weights[i];
    }
    // Normalize weights
    for (size_t i = 0; i < k; ++i) {
      weights[i] /= weightSum;
    }

    // For k < 3, use inverse distance weighting
    if (k < 3) {
      outValues[iloc] =
          computeInverseDistanceWeighting(neighborData.values, weights);
      ++idwKLessThan3;
      continue;
    }

    // For k >= 3, attempt weighted least squares plane fit
    auto fitResult = computePlaneFit<GatherType>(
        iloc, lat_q, lon_q, relativeErrorThreshold, neighborData.values,
        neighborData.latitudes, neighborData.longitudes, weights);

    if (fitResult.success) {
      outValues[iloc] = fitResult.value;
      ++planeFitSuccess;
    } else {
      // Plane fit failed, fall back to IDW
      outValues[iloc] =
          computeInverseDistanceWeighting(neighborData.values, weights);

      // Update fallback counters based on failure reason
      switch (fitResult.reason) {
        case PlaneFitResult<GatherType>::FailureReason::LDLT:
          ++fallbackLDLT;
          break;
        case PlaneFitResult<GatherType>::FailureReason::RelativeError:
          ++fallbackRelativeError;
          break;
        case PlaneFitResult<GatherType>::FailureReason::Exception:
          ++fallbackException;
          break;
        default:
          break;
      }
    }
  }

  // Write output to ObsSpace
  data_.obsspace().put_db(outputVariable.group(), outputVariable.variable(),
                          outValues);

  size_t totalFallback =
      fallbackLDLT + fallbackRelativeError + fallbackException;
  oops::Log::info()
      << "Use nearest neighbors local plane fit algorithm complete. "
      << "Processed " << nLocsThisRank << " locations: "
      << "colocated=" << colocatedPoints << ", "
      << "plane fit=" << planeFitSuccess << ", "
      << "IDW (k<3)=" << idwKLessThan3 << ", "
      << "IDW fallback=" << totalFallback << " (LDLT failure: " << fallbackLDLT
      << ", rel. error exceeded: " << fallbackRelativeError
      << ", other exceptions: " << fallbackException << ")" << std::endl;
}

// Explicit template instantiations for supported type combinations

// GatherType = float
template void LocalPlaneFit::dispatchOnIDType<float>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;

template void LocalPlaneFit::dispatchOnMatchType<float, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<float, int64_t>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<float, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<float, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;

template void LocalPlaneFit::executeImpl<float, int, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, int, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, int, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, std::string, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, std::string, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, std::string, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, util::DateTime, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, util::DateTime, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<float, util::DateTime, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;

// GatherType = int
template void LocalPlaneFit::dispatchOnIDType<int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;

template void LocalPlaneFit::dispatchOnMatchType<int, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<int, int64_t>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<int, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<int, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;

template void LocalPlaneFit::executeImpl<int, int, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, int, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, int, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, std::string, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, std::string, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, std::string, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, util::DateTime, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, util::DateTime, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int, util::DateTime, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;

// GatherType = int64_t
template void LocalPlaneFit::dispatchOnIDType<int64_t>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;

template void LocalPlaneFit::dispatchOnMatchType<int64_t, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<int64_t, int64_t>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<int64_t, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;
template void LocalPlaneFit::dispatchOnMatchType<int64_t, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &, ioda::ObsDtype) const;

template void LocalPlaneFit::executeImpl<int64_t, int, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, int, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, int, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, std::string, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, std::string, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, std::string, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, util::DateTime, int>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void LocalPlaneFit::executeImpl<int64_t, util::DateTime, std::string>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;
template void
LocalPlaneFit::executeImpl<int64_t, util::DateTime, util::DateTime>(
    const UseNearestNeighborsLocalPlaneFitParameters &,
    const UseNearestNeighborsParameters &) const;

static UseNearestNeighborsAlgorithmMaker<LocalPlaneFit> makerLocalPlaneFit(
    "local plane fit");

}  // namespace ufo
