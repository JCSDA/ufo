/*
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_FINDNEARESTNEIGHBORS_H_
#define UFO_FILTERS_FINDNEARESTNEIGHBORS_H_

#include <variant>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "oops/util/NamedEnumerator.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
}  // namespace ioda

namespace ufo {

enum class DistanceMethod { HAVERSINE };

struct DistanceMethodParameterTraitsHelper {
  typedef DistanceMethod EnumType;
  static constexpr char enumTypeName[] = "DistanceMethod";
  static constexpr util::NamedEnumerator<DistanceMethod> namedValues[] = {
      {DistanceMethod::HAVERSINE, "haversine"}};
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::DistanceMethod>
    : public EnumParameterTraits<ufo::DistanceMethodParameterTraitsHelper> {};

}  // namespace oops

namespace ufo {

/// Parameters controlling the operation of the FindNearestNeighbors filter.
class FindNearestNeighborsParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(FindNearestNeighborsParameters, FilterParametersBase)

 public:
  /// Query-point selector: locations with non-missing values mark query points
  /// for which nearest neighbors are found. Coordinates are read from
  /// MetaData/latitude and MetaData/longitude. Must be of float type.
  oops::RequiredParameter<Variable> queryPointVar{"query point variable", this};

  /// Reference-point selector: locations with non-missing values mark candidate
  /// nearest neighbors. Coordinates are read from MetaData/latitude and
  /// MetaData/longitude. Must be of float type.
  oops::RequiredParameter<Variable> referencePointVar{
      "reference point variable", this};

  /// Output variable names for distances to nearest neighbors. Must have the
  /// same length as outputVariables.
  oops::RequiredParameter<std::vector<Variable>> distanceOutputVariables{
      "distance output variables", this};

  /// Values from this variable at reference points are moved from the reference
  /// point locations in the obs space to the query point locations in the obs
  /// space according to the nearest neighbor search results.
  oops::RequiredParameter<Variable> outputAssignment{"output assignment", this};

  /// Output variable names receiving nearest-neighbor values from
  /// outputAssignment. Must have the same length as distanceOutputVariables.
  oops::RequiredParameter<std::vector<Variable>> outputVariables{
      "output variables", this};

  /// Nearest-neighbor search algorithm. Currently only "brute force" is
  /// implemented.
  oops::Parameter<std::string> algorithm{"algorithm", "brute force", this};

  /// Distance method used to rank nearest neighbors. Currently only
  /// "haversine" is implemented.
  oops::Parameter<DistanceMethod> distanceMethod{
      "distance method", DistanceMethod::HAVERSINE, this};

  /// Units for distance outputs: m/metres/meters, km/kilometres/kilometers,
  /// mi/miles, nmi/nautical miles.
  oops::Parameter<std::string> distanceUnits{"distance units", "m", this};
};

/// \brief FindNearestNeighbors filter: finds nearest neighbors for query
/// points.
///
/// For each query point with valid coordinates (latitude/longitude), this
/// filter finds the nearest reference point(s) and outputs their distances and
/// associated values. The search uses a configurable distance-method parameter
/// (currently only haversine is implemented) and various distance units
/// (meters, kilometers, miles, nautical miles).
///
/// The filter performs a rearrangement of observation data or metadata in the
/// obs space: it transfers values from reference point locations to query point
/// locations in the obs space based on nearest neighbor relationships.
///
/// The filter requires specification of:
/// - Query point variable: a variable with non-missing values marking query
/// point locations
/// - Reference point variable: a variable with non-missing values marking
/// reference point
///   locations
/// - Output assignment variable: values at reference point locations to be
/// transferred to
///   query point locations
/// - Output variables receiving the transferred values and their distances
///
/// Coordinates for distance calculations are taken from MetaData/latitude and
/// MetaData/longitude.
///
/// Example YAML configuration:
///
/// \code{.yaml}
/// - filter: Find Nearest Neighbors
///   query point variable: MetaData/station_location
///   reference point variable: MetaData/observation_location
///   output assignment: ObsValue/temperature
///   output variables:
///   - name: DerivedMetaData/nearestTemperature
///   - name: DerivedMetaData/secondNearestTemperature
///   distance output variables:
///   - name: DerivedMetaData/distanceToNearest
///   - name: DerivedMetaData/distanceToSecondNearest
///   distance units: km
/// \endcode
///
/// In this example, for each station (query point), the filter finds the two
/// nearest observation locations (reference points) and transfers their
/// temperature values into the obs space at the station locations, along with
/// the corresponding distances.
class FindNearestNeighbors : public FilterBase,
                             private util::ObjectCounter<FindNearestNeighbors> {
 public:
  using Parameters_ = FindNearestNeighborsParameters;

  static const std::string classname() { return "ufo::FindNearestNeighbors"; }

  FindNearestNeighbors(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                       ioda::ObsDataVector<int> &flags,
                       ioda::ObsDataVector<float> &obserr);
  ~FindNearestNeighbors();

 private:
  // Standard filter methods
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &apply, const Variables &filtervars,
                   std::vector<std::vector<bool>> &flagged) const override;
  int qcFlag() const override { return QCflags::pass; }
  Parameters_ parameters_;

  // Specific methods and types for this filter
  template <typename T>
  void process(const std::vector<T> &outputAssignmentData,
               const std::vector<bool> &apply) const;

  template <typename T>
  struct Neighbor {
    float distance;
    T referenceValue;
  };

  template <typename T>
  std::vector<std::vector<Neighbor<T>>> computeNearestNeighbors(
      const std::vector<std::vector<float>> &queryCoords,
      const std::vector<std::vector<float>> &referenceCoords,
      const std::vector<T> &referenceValues, size_t numNeighbors,
      const DistanceMethod distanceMethod) const;

  float calculateDistanceInMetres(const float lat1, const float lon1,
                                  const float lat2, const float lon2,
                                  const DistanceMethod method) const;

  float convertDistanceFromMeters(const float distance,
                                  const std::string &units) const;

  void validateParameterValues() const;

  std::tuple<std::vector<float>, std::vector<float>, std::vector<size_t>>
  extractValidQueryPoints(const std::vector<bool> &apply,
                          const size_t nLocsThisRank,
                          const std::vector<float> &lats,
                          const std::vector<float> &lons) const;

  template <typename T>
  std::tuple<std::vector<T>, std::vector<float>, std::vector<float>>
  extractValidReferencePoints(const std::vector<bool> &apply,
                              const size_t nLocsThisRank,
                              const std::vector<bool> &isOwnedByThisRank,
                              const std::vector<float> &lats,
                              const std::vector<float> &lons,
                              const std::vector<T> &outputAssignmentData) const;

  template <typename T>
  std::tuple<std::vector<T>, std::vector<float>, std::vector<float>>
  gatherGlobalReferencePoints(const std::vector<T> &localData,
                              const std::vector<float> &localLats,
                              const std::vector<float> &localLons) const;

  std::vector<size_t> uniqueLatLonIndices(const std::vector<float> &lats,
                                          const std::vector<float> &lons) const;

  template <typename T>
  std::tuple<std::vector<T>, std::vector<float>, std::vector<float>>
  removeDuplicateReferencePoints(const std::vector<T> &data,
                                 const std::vector<float> &lats,
                                 const std::vector<float> &lons) const;

  template <typename T>
  std::pair<std::vector<std::vector<float>>, std::vector<std::vector<T>>>
  reconstructNeighborValuesAndDistances(
      const std::vector<std::vector<Neighbor<T>>> &nearestNeighbors,
      const std::vector<size_t> &validQueryIndices, const size_t nLocsThisRank,
      const size_t numNeighbors) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FINDNEARESTNEIGHBORS_H_
