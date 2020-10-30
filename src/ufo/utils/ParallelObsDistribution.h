/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_PARALLELOBSDISTRIBUTION_H_
#define UFO_UTILS_PARALLELOBSDISTRIBUTION_H_

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"

namespace ufo {

/// \brief Describes how observations in an ObsSpace are distributed across MPI processes.
class ParallelObsDistribution {
 public:
  /// \brief Construct an object describing the distribution of observations in \p obsspace across
  /// MPI processes.
  explicit ParallelObsDistribution(const ioda::ObsSpace &obsspace);

  /// \brief Return the total number of observations held by all MPI processes.
  size_t globalObsCount() const { return globalObsCount_; }

  /// \brief Return a vector whose ith element is the number of observations held by the MPI process
  /// with rank i.
  ///
  /// The returned vector can be passed to the \c recvcounts parameter of
  /// eckit::mpi::Comm::allGatherv().
  const std::vector<int> &localObsCounts() const { return localObsCounts_; }

  /// \brief Return a vector whose ith element is the total number of observations held by the MPI
  /// processes with ranks less than i.
  ///
  /// The returned vector can be passed to the \c recvcounts parameter of
  /// eckit::mpi::Comm::displs().
  const std::vector<int> &localObsIdDisplacements() const { return localObsIdDisplacements_; }

 private:
  size_t globalObsCount_;
  std::vector<int> localObsCounts_;
  std::vector<int> localObsIdDisplacements_;
};

/// \brief Return a vector containing the values of variable \c variable@group for all observations
/// held by process 0, then all observations held by process 1 etc.
///
/// \tparam T
///   Type of the variable values. Must be int, float, double or util::DateTime,
///   otherwise a linking error will occur.
///
/// \related ParallelObsDistribution
template <typename T>
std::vector<T> getGlobalVariableValues(const ioda::ObsSpace &obsspace,
                                       const ParallelObsDistribution &obsDistribution,
                                       const std::string &variable,
                                       const std::string &group);

}  // namespace ufo

#endif  // UFO_UTILS_PARALLELOBSDISTRIBUTION_H_
