/*
 * (C) Crown copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_STATISTIC_H_
#define UFO_FILTERS_OBSFUNCTIONS_STATISTIC_H_

#include <functional>
#include <string>
#include <vector>

#include "oops/mpi/mpi.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Enumeration of available statistic methods
enum class StatisticMethod {
  ARITHMETIC_MEAN,
  HARMONIC_MEAN,
  MEDIAN,
  MODE,
  WEIGHTED_MEAN,
  STANDARD_DEVIATION,
  VARIANCE
};

/// \brief Parameter traits helper for StatisticMethod enum
struct StatisticMethodParameterTraitsHelper {
  typedef StatisticMethod EnumType;
  static constexpr char enumTypeName[] = "StatisticMethod";
  static constexpr util::NamedEnumerator<StatisticMethod> namedValues[] = {
      {StatisticMethod::ARITHMETIC_MEAN, "arithmetic mean"},
      {StatisticMethod::HARMONIC_MEAN, "harmonic mean"},
      {StatisticMethod::MEDIAN, "median"},
      {StatisticMethod::MODE, "mode"},
      {StatisticMethod::WEIGHTED_MEAN, "weighted mean"},
      {StatisticMethod::STANDARD_DEVIATION, "standard deviation"},
      {StatisticMethod::VARIANCE, "variance"}};
};

}  // namespace ufo

namespace oops {

/// \brief Register StatisticMethod enum with oops parameter system
template <>
struct ParameterTraits<ufo::StatisticMethod>
    : public EnumParameterTraits<ufo::StatisticMethodParameterTraitsHelper> {};

}  // namespace oops

namespace ufo {

/// \brief Parameters for the Statistic ObsFunction
class StatisticParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StatisticParameters, Parameters)

 public:
  /// Input variable for which to compute the statistic
  oops::RequiredParameter<Variable> variable{"variable", this};

  /// Type of statistic to compute
  oops::RequiredParameter<StatisticMethod> statistic{"statistic", this};

  /// Weight variable (required only for weighted mean)
  oops::OptionalParameter<Variable> weightVariable{"weight variable", this};

  /// Delta degrees of freedom for standard deviation and variance (default 0
  /// for population)
  oops::Parameter<int> deltaDegreesOfFreedom{"delta degrees of freedom", 0,
                                             this};

  /// For harmonic mean: whether to throw an exception on invalid values
  /// (zero or negative). If false, returns missing and logs warning instead.
  /// Default is true (population behavior).
  oops::Parameter<bool> abortIfInvalidOperation{"abort if invalid operation",
                                                true, this};
};

/// \brief Compute global statistics across all MPI ranks
///
/// This ObsFunction computes various statistics (arithmetic mean, harmonic
/// mean, median, mode, weighted mean, standard deviation, variance) by
/// gathering non-missing values from all MPI ranks. The same global statistic
/// value is assigned to all locations; the Variable Assignment filter's where
/// clause controls which locations receive the assignment.
///
/// Note, at present this ObsFunction performs an MPI AllGather of all
/// non-missing values to compute the statistic, which may result in high memory
/// usage for large datasets.
///
/// For weighted mean, both the variable and weight must be non-missing for a
/// value to contribute to the statistic.
///
/// Mode behavior:
/// - If no modal value exists (all values unique), returns missing
/// - If multiple modes exist (multimodal), returns the smallest modal value
///
/// Harmonic mean behavior:
/// - If any zero or negative values are encountered, logs a warning and returns
/// missing
class Statistic : public ObsFunctionBase<float> {
 public:
  /// Constructor takes configuration and validates parameters
  explicit Statistic(const eckit::LocalConfiguration &);

  /// Destructor
  ~Statistic();

  /// Compute the statistic across all MPI ranks
  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;

  /// Return the required variables for this ObsFunction
  const ufo::Variables &requiredVariables() const;

 private:
  /// Compute arithmetic mean across all ranks
  void computeArithmeticMean(const ObsFilterData &,
                             ioda::ObsDataVector<float> &) const;

  /// Compute harmonic mean across all ranks
  void computeHarmonicMean(const ObsFilterData &,
                           ioda::ObsDataVector<float> &) const;

  /// Compute median across all ranks
  void computeMedian(const ObsFilterData &, ioda::ObsDataVector<float> &) const;

  /// Compute mode across all ranks
  void computeMode(const ObsFilterData &, ioda::ObsDataVector<float> &) const;

  /// Compute weighted mean across all ranks
  void computeWeightedMean(const ObsFilterData &,
                           ioda::ObsDataVector<float> &) const;

  /// Compute variance across all ranks
  void computeVariance(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const;

  /// Compute standard deviation across all ranks (calls computeVariance)
  void computeStdDev(const ObsFilterData &, ioda::ObsDataVector<float> &) const;

  /// Helper method to gather global non-missing owned values for a given
  /// channel
  std::vector<float> gatherGlobalNonMissingValues(
      const ObsFilterData &in, const ioda::ObsDataVector<float> &data,
      const std::vector<bool> &isOwnedByThisRank, size_t channelIndex) const;

  /// Template method to reduce boilerplate in statistic computation methods
  /// Handles getting input data, ownership, looping over variables, and
  /// assigning results
  void applyStatisticFunction(
      const ObsFilterData &in, ioda::ObsDataVector<float> &out,
      std::function<float(const std::vector<float> &)> computeFunc) const;

  /// Parameters for this ObsFunction
  StatisticParameters options_;

  /// Required variables
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_STATISTIC_H_
