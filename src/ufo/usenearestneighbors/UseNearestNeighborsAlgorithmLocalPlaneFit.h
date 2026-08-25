/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMLOCALPLANEFIT_H_
#define UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMLOCALPLANEFIT_H_

#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

/// \brief Parameters for the "local plane fit" algorithm.
///
/// For each query point, this algorithm interpolates (or extrapolates) a scalar
/// value from its nearest reference observations using weighted least-squares
/// plane fitting in a local equirectangular coordinate system centered on the
/// query point. If the plane fit is not well-determined (LDLT failure, fewer
/// than 3 neighbors, or large relative residuals), the algorithm falls back to
/// inverse distance weighting (IDW).
///
/// This algorithm uses the query and reference point variables in two ways:
/// (1) to identify which observations are query/reference points (non-missing
/// values), and (2) to indicate the locations whose latitude/longitude coordinates
/// are used for the plane fit. This is similar to the Find Nearest Neighbors filter,
/// where these variables identify both the observation type and the coordinate source.
///
/// The \c match_variable ensures that query and reference values are compared
/// at equivalent conditions (e.g. the same datetime). Only reference
/// observations whose \c match_variable value equals that of the query point
/// are used.
///
/// Example YAML:
/// \code{.yaml}
/// algorithm:
///   name: local plane fit
///   query point variable: ObsValue/airTemperatureObservation
///   reference point variable: ObsValue/airTemperatureReference
///   match variable: MetaData/matchDateTime
///   inverse distance weighting power: 2.0
///   relative error threshold: 0.25
///   distance variables:
///     - name: DerivedMetaData/firstNearestReferenceDistance
///     - name: DerivedMetaData/secondNearestReferenceDistance
///     - name: DerivedMetaData/thirdNearestReferenceDistance
///     - name: DerivedMetaData/fourthNearestReferenceDistance
///   output variable: DerivedMetaData/interpolatedTemperature
/// \endcode
///
/// For each ObsSpace row \c i:
/// - The reference location coordinates and values are looked up by matching the
///   nearest-neighbor identifier variable values at row \c i with reference data
///   whose \c match_variable equals the value at row \c i. For example, if
///   \c firstNearestReferenceStationID[i] = "A123" and \c matchIndex[i] = 1, the
///   algorithm searches for all reference rows with
///   \c referenceStationIdentification = "A123" and \c matchIndex = 1.
/// - Reference row coordinates (latitude, longitude) and the
///   \c reference_point_variable value are extracted; distances from (lat, lon)
///   at row \c i are read from the \c distance variables.
/// - If fewer than 3 matched neighbors, IDW of the neighbor values is written to
///   the output variable at row \c i.
/// - If 3 or more matched neighbors, a weighted least-squares plane fit is attempted.
///   If the fit succeeds (positive-definite LDLT and relative error below threshold),
///   the intercept (fitted value at the query location) is written. Otherwise,
///   IDW is written. If any required data is missing at row \c i, the output
///   remains missing at row \c i.
class UseNearestNeighborsLocalPlaneFitParameters
    : public UseNearestNeighborsAlgorithmParametersBase {
  OOPS_CONCRETE_PARAMETERS(UseNearestNeighborsLocalPlaneFitParameters,
                           UseNearestNeighborsAlgorithmParametersBase)
 public:
  /// \brief Variable at query-point locations. Locations where this variable
  /// is missing are skipped and are not used for interpolation. The latitude and
  /// longitude coordinates at these non-missing locations are used to define the
  /// query points for the plane fit. Supported types: float, integer.
  oops::RequiredParameter<Variable> queryPointVariable{"query point variable",
                                                       this};

  /// \brief Variable at reference-point locations to be
  /// interpolated/extrapolated. Locations where this variable is missing are
  /// excluded from the reference set. The latitude and longitude coordinates
  /// at these non-missing locations are used as reference points for the plane fit.
  /// If this variable is missing at any required reference point identified by
  /// the nearest-neighbor matching, the output remains missing at the query point.
  oops::RequiredParameter<Variable> referencePointVariable{
      "reference point variable", this};

  /// \brief Variables holding the distances (in km) from each query point to
  /// each of its nearest reference observations, one per nearest-neighbor
  /// identifier variable. Must be the same length as \c nearest neighbor
  /// identifier variables.
  oops::RequiredParameter<std::vector<Variable>> distanceVariables{
      "distance variables", this};

  /// \brief Power parameter \f$p\f$ for inverse distance weighting. Weights
  /// are proportional to \f$1 / d^p\f$.
  oops::RequiredParameter<float> inverseDistanceWeightingPower{
      "inverse distance weighting power", this};

  /// \brief Variable where the interpolated/extrapolated values are written.
  /// Must not be in the ObsValue or DerivedObsValue group.
  oops::RequiredParameter<Variable> outputVariable{"output variable", this};

  /// \brief Relative weighted-2-norm error threshold for the plane fit quality
  /// check. If the ratio \f$\sqrt{\sum_k w_k e_k^2} / (\sqrt{\sum_k w_k z_k^2} +
  /// 10^{-10})\f$ (where \f$e = z - A\beta\f$) exceeds this value, the
    /// algorithm falls back to IDW. Default is 0.25. Setting this to 0.0
    /// accepts only exact zero-relative-error fits; all other fits fall back to
    /// IDW.
  oops::OptionalParameter<float> relativeErrorThreshold{
      "relative error threshold", this};

  /// \brief Variable used to match query points with reference points. Only
  /// reference observations whose \c match variable value equals that of the
  /// query point are used. Supported types: integer, 64-bit integer, string,
  /// datetime.
  oops::RequiredParameter<Variable> matchVariable{"match variable", this};
};

/// \brief Implementation of the "local plane fit" nearest-neighbor algorithm.
/// See UseNearestNeighborsLocalPlaneFitParameters for the full description and
/// YAML example.
class UseNearestNeighborsAlgorithmLocalPlaneFit
    : public UseNearestNeighborsAlgorithmBase {
 public:
  using Parameters_ = UseNearestNeighborsLocalPlaneFitParameters;

  /// \brief Lookup map type: (ID, match value) → (reference value, latitude,
  /// longitude).
  template <typename GatherType, typename IDType, typename MatchType>
  using LookupMap = std::unordered_map<std::pair<IDType, MatchType>,
                                       std::tuple<GatherType, float, float>,
                                       PairHash<IDType, MatchType>>;

  UseNearestNeighborsAlgorithmLocalPlaneFit(
      const UseNearestNeighborsLocalPlaneFitParameters& params,
      const ObsFilterData& data, const std::vector<bool>& apply,
      const Variables& filtervars, const ioda::ObsDataVector<int>& flags,
      std::vector<std::vector<bool>>& flagged);

 private:
  // Helper structs

  /// \brief Data gathered from the nearest neighbors of a single query point.
  template <typename GatherType>
  struct NeighborData {
    bool isValid;                       ///< False if any required data is missing.
    std::vector<GatherType> values;     ///< Reference-point variable values.
    std::vector<float> latitudes;       ///< Reference-point latitudes (degrees).
    std::vector<float> longitudes;      ///< Reference-point longitudes (degrees).
    std::vector<float> distances;       ///< Distances from query to each neighbor (km).
  };

  /// \brief Result of a local plane fit at a single query point.
  template <typename GatherType>
  struct PlaneFitResult {
    bool success;          ///< True if the plane fit was accepted.
    GatherType value;      ///< Fitted value at the query point (intercept).
    /// \brief Reason the plane fit was not accepted (if \c success is false).
    enum class FailureReason { None, LDLT, RelativeError, Exception };
    FailureReason reason;
  };

  void execute(const UseNearestNeighborsAlgorithmParametersBase& algParams,
               const UseNearestNeighborsParameters& options) const override;

  /// \brief Dispatch on the ID variable type after the gather type is known.
  template <typename GatherType>
  void dispatchOnIDType(
      const UseNearestNeighborsLocalPlaneFitParameters& params,
      const UseNearestNeighborsParameters& options,
      const ioda::ObsDtype idDtype) const;

  /// \brief Dispatch on the match variable type after gather and ID types are
  /// known.
  template <typename GatherType, typename IDType>
  void dispatchOnMatchType(
      const UseNearestNeighborsLocalPlaneFitParameters& params,
      const UseNearestNeighborsParameters& options,
      const ioda::ObsDtype matchDtype) const;

  /// \brief Fully typed implementation of the local plane fit algorithm.
  template <typename GatherType, typename IDType, typename MatchType>
  void executeImpl(const UseNearestNeighborsLocalPlaneFitParameters& params,
                   const UseNearestNeighborsParameters& options) const;

  /// \brief Gather reference data from all MPI ranks and build a lookup map
  /// keyed on (station ID, match value) → (reference value, latitude,
  /// longitude). Throws if the same key appears more than once (the match
  /// variable must be unique per reference identifier).
  template <typename GatherType, typename IDType, typename MatchType>
  LookupMap<GatherType, IDType, MatchType> gatherReferenceData(
      const Variable& idVariable, const Variable& referencePointVariable,
      const std::vector<bool>& isOwnedByThisRank,
      const std::vector<IDType>& idValues,
      const std::vector<GatherType>& referencePointValues,
      const std::vector<float>& latValues, const std::vector<float>& lonValues,
      const std::vector<MatchType>& matchValues) const;

  /// \brief Collect neighbor data for a single query-point location \p iloc.
  /// Returns an invalid NeighborData if any required neighbor field is missing.
  template <typename GatherType, typename IDType, typename MatchType>
  NeighborData<GatherType> collectNeighborData(
      size_t iloc, size_t numNeighbors,
      const std::vector<std::vector<IDType>>& allNearestNeighborIDs,
      const std::vector<std::vector<float>>& allDistances,
      const LookupMap<GatherType, IDType, MatchType>& lookupMap,
      MatchType matchValue) const;

  /// \brief Attempt a weighted least-squares local plane fit for a single
  /// query point. Returns a PlaneFitResult indicating success or the reason
  /// for fallback to IDW.
  template <typename GatherType>
  PlaneFitResult<GatherType> computePlaneFit(
      size_t iloc, float lat_q, float lon_q, float relativeErrorThreshold,
      const std::vector<GatherType>& neighborValues,
      const std::vector<float>& neighborLats,
      const std::vector<float>& neighborLons,
      const std::vector<float>& weights) const;
};

}  // namespace ufo

#endif  // UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMLOCALPLANEFIT_H_
