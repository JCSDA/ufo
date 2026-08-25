/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMGATHERANDMATCHTIMESTAMP_H_
#define UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMGATHERANDMATCHTIMESTAMP_H_

#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

/// \brief Parameters for the "gather and match timestamp" algorithm.
///
/// For each query point, this algorithm gathers a value from each nearest
/// reference observation using an exact match on identifier and timestamp. One
/// output variable is written per nearest-neighbor identifier variable (i.e.
/// output variable 1 corresponds to the first nearest neighbor, output variable 2
/// to the second nearest neighbor, and so on).
///
/// The algorithm builds a global lookup map over reference candidates satisfying
/// \c timestamp_match_variable == MetaData/dateTime, then for each location and
/// neighbor slot looks up the pair (neighbor ID, timestamp match variable).
///
/// Example YAML:
/// \code{.yaml}
/// algorithm:
///   name: gather and match timestamp
///   gather variable: ObsValue/airTemperatureReference
///   timestamp match variable: DerivedMetaData/binnedDateTime
///   output variables:
///     - name: DerivedMetaData/matchedFirstNearestTemperature
///     - name: DerivedMetaData/matchedSecondNearestTemperature
///     - name: DerivedMetaData/matchedThirdNearestTemperature
/// \endcode
///
/// For each ObsSpace row \c i:
/// - \c DerivedMetaData/matchedFirstNearestTemperature[i] is looked up using
///   the key ( \c DerivedMetaData/firstNearestReferenceStationID[i],
///   \c DerivedMetaData/binnedDateTime[i] ), matched against rows where
///   \c MetaData/referenceStationIdentification and
///   \c DerivedMetaData/binnedDateTime == \c MetaData/dateTime.
/// - \c DerivedMetaData/matchedSecondNearestTemperature[i] and
///   \c DerivedMetaData/matchedThirdNearestTemperature[i] are looked up
///   analogously using the second and third nearest-neighbor station IDs.
/// - If any key is absent from the lookup map, the corresponding output is
///   written as missing at row \c i .
class UseNearestNeighborsGatherAndMatchTimestampParameters
    : public UseNearestNeighborsAlgorithmParametersBase {
  OOPS_CONCRETE_PARAMETERS(UseNearestNeighborsGatherAndMatchTimestampParameters,
                           UseNearestNeighborsAlgorithmParametersBase)
 public:
  /// \brief The variable to gather from the reference observations. Supported
  /// types: float, integer, string, datetime.
  oops::RequiredParameter<Variable> gatherVariable{"gather variable", this};

  /// \brief Variable holding the timestamp used to select which value to
  /// retrieve from each nearest reference observation. The value at the query
  /// point location in the ObsSpace is matched against values at the reference
  /// observation locations (where the identifier variable matches) and only
  /// reference observations for which
  /// \c timestamp_match_variable == \c MetaData/dateTime are considered.
  /// Supported type: datetime.
  oops::RequiredParameter<Variable> timestampMatchVariable{
      "timestamp match variable", this};

  /// \brief Output variables where matched values will be written, one per
  /// entry in \c nearest_neighbor_identifier_variables. Must be the same
  /// length as that list. Variables in the ObsValue or DerivedObsValue groups
  /// are not permitted.
  oops::RequiredParameter<std::vector<Variable>> outputVariables{
      "output variables", this};
};

/// \brief Implementation of the "gather and match timestamp" nearest-neighbor
/// algorithm. See UseNearestNeighborsGatherAndMatchTimestampParameters for the
/// full description and YAML example.
class UseNearestNeighborsAlgorithmGatherAndMatchTimestamp
    : public UseNearestNeighborsAlgorithmBase {
 public:
  using Parameters_ = UseNearestNeighborsGatherAndMatchTimestampParameters;

  UseNearestNeighborsAlgorithmGatherAndMatchTimestamp(
      const UseNearestNeighborsGatherAndMatchTimestampParameters& params,
      const ObsFilterData& data, const std::vector<bool>& apply,
      const Variables& filtervars, const ioda::ObsDataVector<int>& flags,
      std::vector<std::vector<bool>>& flagged);

 private:
  void execute(const UseNearestNeighborsAlgorithmParametersBase& algParams,
               const UseNearestNeighborsParameters& options) const override;

  template <typename GatherType>
  void dispatchOnIDType(
      const UseNearestNeighborsGatherAndMatchTimestampParameters& params,
      const UseNearestNeighborsParameters& options,
      ioda::ObsDtype idDtype) const;

  template <typename GatherType, typename IDType>
  void executeImpl(
      const UseNearestNeighborsGatherAndMatchTimestampParameters& params,
      const UseNearestNeighborsParameters& options) const;
};

}  // namespace ufo

#endif  // UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMGATHERANDMATCHTIMESTAMP_H_
