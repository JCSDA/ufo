/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMREFERENCEPOINTVARIABLESMEAN_H_
#define UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMREFERENCEPOINTVARIABLESMEAN_H_

#include <map>
#include <utility>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/Variable.h"
#include "UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

/// \brief Parameters for the "reference point variables mean" algorithm.
///
/// For each query point, this algorithm gathers a variable from its nearest
/// reference observations and computes means grouped by an integer
/// \c mean_about_variable. The first output variable receives the mean for bin
/// 0, the second for bin 1, etc. (determined by list position, not by variable
/// name).
///
/// It is assumed that, for a given combination of identifier and
/// \c mean_about_variable integer, there is at most one unique
/// \c gather_variable value.
///
/// Example YAML:
/// \code{.yaml}
/// algorithm:
///   name: reference point variables mean
///   gather variable: ObsValue/airTemperatureReference
///   mean about: DerivedMetaData/dateTimeAveragingBin # 0, 1, etc. (only 0 and
///   1 used here) output variables:
///     - name: DerivedMetaData/averageThreeClosestTemperaturesMostRecent  # bin
///     0
///     - name: DerivedMetaData/averageThreeClosestTemperaturesPrevious    # bin
///     1
/// \endcode
///
/// For each ObsSpace row \c i:
/// - The nearest reference observations are identified by the nearest-neighbor
///   identifier variable values at row \c i (e.g.
///   \c DerivedMetaData/firstNearestReferenceStationID[i] ,
///   \c DerivedMetaData/secondNearestReferenceStationID[i] ,
///   \c DerivedMetaData/thirdNearestReferenceStationID[i]) .
/// - \c ObsValue/airTemperatureReference values are gathered from the
///   corresponding reference rows and grouped by their
///   \c DerivedMetaData/dateTimeAveragingBin value.
/// - \c DerivedMetaData/averageThreeClosestTemperaturesMostRecent[i] is the
///   mean of the gathered values for which the bin index is 0.
/// - \c DerivedMetaData/averageThreeClosestTemperaturesPrevious[i] is the
///   mean for bin index 1.
/// - If no reference observations fall into a bin, the corresponding output
///   is written as missing at row \c i .
class UseNearestNeighborsReferencePointVariablesMeanParameters
    : public UseNearestNeighborsAlgorithmParametersBase {
  OOPS_CONCRETE_PARAMETERS(
      UseNearestNeighborsReferencePointVariablesMeanParameters,
      UseNearestNeighborsAlgorithmParametersBase)
 public:
  /// \brief Variable to gather from the reference observations and average.
  /// Supported types: float, integer, datetime.
  oops::RequiredParameter<Variable> gatherVariable{"gather variable", this};

  /// \brief Integer variable defining the bin index used to group values
  /// before averaging.
  oops::RequiredParameter<Variable> meanAboutVariable{"mean about", this};

  /// \brief Output variables where per-bin means are written, one per bin.
  /// The i-th output variable (0-indexed) receives the mean for bin i. The
  /// variable names are labels only; the mapping is determined by list
  /// position. Variables in the ObsValue or DerivedObsValue groups are not
  /// permitted.
  oops::RequiredParameter<std::vector<Variable>> outputVariables{
      "output variables", this};
};

/// \brief Implementation of the "reference point variables mean"
/// nearest-neighbor algorithm. See
/// UseNearestNeighborsReferencePointVariablesMeanParameters for the full
/// description and YAML example.
class UseNearestNeighborsAlgorithmReferencePointVariablesMean
    : public UseNearestNeighborsAlgorithmBase {
 public:
  using Parameters_ = UseNearestNeighborsReferencePointVariablesMeanParameters;

  UseNearestNeighborsAlgorithmReferencePointVariablesMean(
      const UseNearestNeighborsReferencePointVariablesMeanParameters& params,
      const ObsFilterData& data, const std::vector<bool>& apply,
      const Variables& filtervars, const ioda::ObsDataVector<int>& flags,
      std::vector<std::vector<bool>>& flagged);

 private:
  void execute(const UseNearestNeighborsAlgorithmParametersBase& algParams,
               const UseNearestNeighborsParameters& options) const override;

  template <typename GatherType>
  void dispatchOnIDType(
      const UseNearestNeighborsReferencePointVariablesMeanParameters& params,
      const UseNearestNeighborsParameters& options, ioda::ObsDtype idDtype,
      ioda::ObsDtype meanAboutDtype) const;

  template <typename GatherType, typename IDType>
  void dispatchOnMeanAboutType(
      const UseNearestNeighborsReferencePointVariablesMeanParameters& params,
      const UseNearestNeighborsParameters& options,
      ioda::ObsDtype meanAboutDtype) const;

  template <typename GatherType, typename IDType, typename MeanAboutType>
  void executeImpl(
      const UseNearestNeighborsReferencePointVariablesMeanParameters& params,
      const UseNearestNeighborsParameters& options) const;

  template <typename GatherType, typename IDType, typename MeanAboutType>
  std::map<std::pair<IDType, MeanAboutType>, GatherType> buildLookupMap(
      const std::vector<IDType>& matchedIDValues,
      const std::vector<GatherType>& matchedGatherValues,
      const std::vector<MeanAboutType>& matchedMeanAboutValues) const;

  template <typename GatherType>
  GatherType computeMean(const std::vector<GatherType>& valuesToAverage) const;

  template <typename GatherType, typename IDType, typename MeanAboutType>
  std::vector<GatherType> collectValuesForAveraging(
      size_t iloc, MeanAboutType targetMeanAboutValue,
      const std::vector<std::vector<IDType>>& allNearestNeighborIDValues,
      const std::map<std::pair<IDType, MeanAboutType>, GatherType>& lookupMap)
      const;
};

}  // namespace ufo

#endif  // UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMREFERENCEPOINTVARIABLESMEAN_H_
