/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMBASE_H_
#define UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMBASE_H_

#include <algorithm>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterParametersBase.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variables.h"

#include "ufo/filters/UseNearestNeighborsParameters.h"

namespace ufo {

/// \brief Hash functor for types used as unordered_map keys, including
/// util::DateTime (which has no std::hash specialization).
template <typename IDType>
struct IDHash {
  std::size_t operator()(const IDType &id) const {
    if constexpr (std::is_same_v<IDType, util::DateTime>) {
      // DateTime needs toString() since it doesn't have std::hash
      return std::hash<std::string>{}(id.toString());
    } else {
      // Other types are assumed to have std::hash
      return std::hash<IDType>{}(id);
    }
  }
};

/// \brief Hash functor for std::pair<FirstType, SecondType>. Combines the
/// hashes from IDHash<FirstType> and IDHash<SecondType>. Needed whenever a
/// pair is used as an unordered_map key.
template <typename FirstType, typename SecondType>
struct PairHash {
  std::size_t operator()(const std::pair<FirstType, SecondType> &p) const {
    // Combine two hashes with a non-commutative mix so (a, b) and (b, a)
    // naturally hash to different values when used as unordered_map keys.
    // This follows boost::hash_combine style mixing (see
    // https://www.boost.org/doc/libs/1_34_1/doc/html/boost/hash_combine.html)
    const auto h1 = IDHash<FirstType>{}(p.first);
    const auto h2 = IDHash<SecondType>{}(p.second);
    auto seed = h1;
    seed ^= h2 + 0x9e3779b9UL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

/// \brief Abstract base class for algorithm parameters for the Use Nearest
/// Neighbors filter. Concrete algorithm parameter classes inherit from this
/// and add their own algorithm-specific parameters. The \c name (YAML key)
/// parameter selects which registered algorithm to create.
class UseNearestNeighborsAlgorithmParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(UseNearestNeighborsAlgorithmParametersBase,
                           oops::Parameters)
 public:
  oops::RequiredParameter<std::string> algorithmName{"name", this};
};

// Forward declaration of top-level filter parameters
class UseNearestNeighborsParameters;

/// \brief Abstract base class for Use Nearest Neighbors algorithms.
///
/// Subclasses implement the \c execute() method to perform the specific
/// nearest-neighbor operation (e.g. gathering a value, computing a mean, or
/// fitting a local plane). The public interface is \c runAlgorithm(), which
/// extracts the concrete algorithm parameter class (e.g.
/// UseNearestNeighborsGatherAndMatchTimestampParameters) from the abstract
/// base type and delegates to \c execute().
///
/// The base class provides helper methods shared by all algorithms:
/// - \c allGatherv(): gathers a rank-local vector to global across all MPI
///   ranks, including special handling for boolean vectors.
/// - \c orderedValidIndices(): finds indices that are owned, non-missing, and
///   not masked out by the \c where clause.
/// - \c findSortedIntersection(): finds the sorted intersection of two or more
///   sorted index vectors.
/// - \c verifyNearestNeighborIDTypesMatch(): checks that all nearest-neighbor
///   identifier variables have the same ObsDtype as the identifier variable.
class UseNearestNeighborsAlgorithmBase {
 public:
  explicit UseNearestNeighborsAlgorithmBase(
      const UseNearestNeighborsAlgorithmParametersBase &, const ObsFilterData &,
      const std::vector<bool> &, const Variables &,
      const ioda::ObsDataVector<int> &, std::vector<std::vector<bool>> &);
  virtual ~UseNearestNeighborsAlgorithmBase() {}

  /// \brief Run the algorithm, dispatching to the concrete \c execute()
  /// implementation with the appropriate typed parameters extracted from
  /// \p options.
  void runAlgorithm(const UseNearestNeighborsParameters &options) const;

  /// \brief Gather a rank-local vector to a global vector across all MPI
  /// ranks. Modifies \p data in place.
  template <typename T>
  void allGatherv(std::vector<T> &data) const;

 protected:
  /// \brief Execute the algorithm. Subclasses must override this method to
  /// implement the specific nearest-neighbor operation. \p algParams is
  /// received as the abstract base type, but the concrete implementation
  /// casts it via dynamic_cast to its specific parameter class for type safety.
  /// \p options holds the top-level filter parameters (common to all algorithms).
  virtual void execute(const UseNearestNeighborsAlgorithmParametersBase &algParams,
                       const UseNearestNeighborsParameters &options) const = 0;

  /// \brief Return a sorted vector of location indices for which \p variableData
  /// is non-missing, the location is owned by this MPI rank, and the location
  /// is selected by the \c where clause (i.e. \p apply is true).
  template <typename T>
  std::vector<size_t> orderedValidIndices(
      const std::vector<T> &variableData,
      const std::vector<bool> &isOwnedByThisRank,
      const std::vector<bool> &apply) const;

  /// \brief Return the sorted intersection of two or more sorted index
  /// vectors. At least two vectors must be provided.
  template <typename... Args>
  std::vector<size_t> findSortedIntersection(const std::vector<size_t> &first,
                                             const std::vector<size_t> &second,
                                             const Args &...rest) const;

  /// \brief Verify that all variables in \p nearestNeighborIDVars have the
  /// same ObsDtype as \p idVariable. Throws if any mismatch is found.
  void verifyNearestNeighborIDTypesMatch(
      const Variable &idVariable,
      const std::vector<Variable> &nearestNeighborIDVars) const;

  const ObsFilterData &data_;
  ioda::ObsSpace &obsdb_;
  const std::vector<bool> apply_;

 private:
  const Variables &filtervars_;
  const ioda::ObsDataVector<int> &flags_;
  std::vector<std::vector<bool>> &flagged_;
};

/// \brief Factory for creating UseNearestNeighbors algorithm objects.
///
/// New algorithms are registered by creating a static instance of
/// UseNearestNeighborsAlgorithmMaker<MyAlgorithm> with the algorithm name
/// string. The factory looks up this name to instantiate the correct concrete
/// class.
class UseNearestNeighborsAlgorithmFactory {
 public:
  static std::unique_ptr<UseNearestNeighborsAlgorithmBase> create(
      const UseNearestNeighborsAlgorithmParametersBase &, const ObsFilterData &,
      const std::vector<bool> &, const Variables &,
      const ioda::ObsDataVector<int> &, std::vector<std::vector<bool>> &);
  static std::unique_ptr<UseNearestNeighborsAlgorithmParametersBase>
  createParameters(const std::string &);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

 protected:
  explicit UseNearestNeighborsAlgorithmFactory(const std::string &);
  virtual ~UseNearestNeighborsAlgorithmFactory() = default;

 private:
  virtual std::unique_ptr<UseNearestNeighborsAlgorithmBase> make(
      const UseNearestNeighborsAlgorithmParametersBase &, const ObsFilterData &,
      const std::vector<bool> &, const Variables &,
      const ioda::ObsDataVector<int> &, std::vector<std::vector<bool>> &) = 0;
  virtual std::unique_ptr<UseNearestNeighborsAlgorithmParametersBase>
  makeParameters() const = 0;

  static std::map<std::string, UseNearestNeighborsAlgorithmFactory *> &
  getMakers() {
    static std::map<std::string, UseNearestNeighborsAlgorithmFactory *> makers_;
    return makers_;
  }
};

/// \brief Template-based concrete maker for UseNearestNeighbors algorithm
/// subclasses.
///
/// Instantiate a static member of this class (e.g. in the algorithm's .cc
/// file) to register the algorithm under its name:
/// \code
///   static UseNearestNeighborsAlgorithmMaker<MyAlgorithm> maker("my name");
/// \endcode
template <class T>
class UseNearestNeighborsAlgorithmMaker
    : public UseNearestNeighborsAlgorithmFactory {
  using Parameters_ = oops::TParameters_IfAvailableElseFallbackType_t<
      T, UseNearestNeighborsAlgorithmParametersBase>;

  std::unique_ptr<UseNearestNeighborsAlgorithmBase> make(
      const UseNearestNeighborsAlgorithmParametersBase &params,
      const ObsFilterData &data, const std::vector<bool> &apply,
      const Variables &filtervars, const ioda::ObsDataVector<int> &flags,
      std::vector<std::vector<bool>> &flagged) override {
    const auto &stronglyTyped = dynamic_cast<const Parameters_ &>(params);
    return std::unique_ptr<UseNearestNeighborsAlgorithmBase>(
        new T(stronglyTyped, data, apply, filtervars, flags, flagged));
  }

  std::unique_ptr<UseNearestNeighborsAlgorithmParametersBase> makeParameters()
      const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit UseNearestNeighborsAlgorithmMaker(const std::string &name)
      : UseNearestNeighborsAlgorithmFactory(name) {}
};

template <typename T>
void UseNearestNeighborsAlgorithmBase::allGatherv(std::vector<T> &data) const {
  // special handling for bool vectors since vector<bool> is not a standard
  // container and does not have a proper data() method which allGatherv needs
  // to access the underlying data.
  if constexpr (std::is_same<T, bool>::value) {
    std::vector<char> boolAsBytes(data.begin(), data.end());
    oops::mpi::allGatherv(obsdb_.comm(), boolAsBytes);
    data.assign(boolAsBytes.begin(), boolAsBytes.end());
  } else {
    oops::mpi::allGatherv(obsdb_.comm(), data);
  }
}

template <typename T>
std::vector<size_t> UseNearestNeighborsAlgorithmBase::orderedValidIndices(
    const std::vector<T> &variableData,
    const std::vector<bool> &isOwnedByThisRank,
    const std::vector<bool> &apply) const {
  std::vector<size_t> validIndices;
  for (size_t i = 0; i < variableData.size(); ++i) {
    if (isOwnedByThisRank[i] && apply[i] &&
        variableData[i] != util::missingValue<T>()) {
      validIndices.push_back(i);
    }
  }
  return validIndices;
}

// Variadic template implementation for findSortedIntersection
// Recursive case: intersect first two vectors, then recurse
template <typename... Args>
std::vector<size_t> UseNearestNeighborsAlgorithmBase::findSortedIntersection(
    const std::vector<size_t> &first, const std::vector<size_t> &second,
    const Args &...rest) const {
  // Intersect first two vectors
  std::vector<size_t> intersection;
  std::set_intersection(first.begin(), first.end(), second.begin(),
                        second.end(), std::back_inserter(intersection));

  // If no more arguments, return the intersection
  if constexpr (sizeof...(rest) == 0) {
    return intersection;
  } else {
    // Otherwise, recurse with the intersection and remaining arguments
    return findSortedIntersection(intersection, rest...);
  }
}

}  // namespace ufo

#endif  // UFO_USENEARESTNEIGHBORS_USENEARESTNEIGHBORSALGORITHMBASE_H_
