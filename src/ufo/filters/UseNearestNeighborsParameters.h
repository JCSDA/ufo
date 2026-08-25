/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_USENEARESTNEIGHBORSPARAMETERS_H_
#define UFO_FILTERS_USENEARESTNEIGHBORSPARAMETERS_H_

#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterParametersBase.h"
#include "ufo/filters/Variable.h"

namespace ufo {

// Forward declarations (full definitions in UseNearestNeighborsAlgorithmBase.h)
class UseNearestNeighborsAlgorithmParametersBase;
class UseNearestNeighborsAlgorithmFactory;

/// \brief Wrapper holding the polymorphic algorithm parameters for the Use
/// Nearest Neighbors filter. The \c name field selects which algorithm to use
/// (e.g. "gather and match timestamp", "reference point variables mean",
/// "local plane fit") and reads the corresponding concrete parameter class.
class UseNearestNeighborsAlgorithmParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(UseNearestNeighborsAlgorithmParametersWrapper,
                           oops::Parameters)
 public:
  oops::RequiredPolymorphicParameter<UseNearestNeighborsAlgorithmParametersBase,
                                     UseNearestNeighborsAlgorithmFactory>
      name{"name", this};
};

/// \brief Top-level filter parameters for the Use Nearest Neighbors filter.
///
/// Example YAML:
/// \code{.yaml}
/// filter: Use Nearest Neighbors
/// identifier variable: MetaData/referenceStationIdentification
/// nearest neighbor identifier variables:
///   - name: DerivedMetaData/firstNearestReferenceStationID
///   - name: DerivedMetaData/secondNearestReferenceStationID
/// algorithm:
///   {...} # See individual algorithm parameter classes for details
/// \endcode
class UseNearestNeighborsParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(UseNearestNeighborsParameters, FilterParametersBase)

 public:
  /// \brief Algorithm wrapper; selects which algorithm to run and holds its
  /// parameters.
  oops::RequiredParameter<UseNearestNeighborsAlgorithmParametersWrapper>
      algorithmParameters{"algorithm",
                          "Algorithm wrapper; selects which algorithm to run "
                          "and holds its parameters.",
                          this};
  /// \brief Variable whose value uniquely identifies each reference
  /// observation (i.e., each ObsSpace row used as a reference point).
  /// This is often a station ID shared by multiple rows (for example,
  /// different times at the same physical location). Values of this variable
  /// on ObsSpace rows not used as reference points may be missing or arbitrary,
  /// as they are not used by the algorithm. Supported types: integer, string,
  /// datetime.
  oops::RequiredParameter<Variable> idVar{
      "identifier variable",
      "Variable whose value uniquely identifies each reference observation.",
      this};
  /// \brief List of variables holding, at each query-point location, the
  /// identifier of the first, second, third, etc. nearest reference
  /// observation row in the ObsSpace (nearest in space, typically produced by
  /// the Find Nearest Neighbors filter). At query row \c i, each value should match
  /// the \c identifier_variable value found on the corresponding reference
  /// row. Must be the same type as \c identifier_variable. The length of this
  /// list determines how many nearest neighbors are considered.
  oops::RequiredParameter<std::vector<Variable>> nearestNeighborIDVars{
      "nearest neighbor identifier variables",
      "Ordered list of nearest-neighbor identifier variables at query points.",
      this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_USENEARESTNEIGHBORSPARAMETERS_H_
