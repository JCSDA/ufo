/*
 * (C) Crown copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/Statistic.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <vector>

#include "ioda/distribution/Distribution.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// Define static members for enum registration
constexpr char StatisticMethodParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<StatisticMethod>
    StatisticMethodParameterTraitsHelper::namedValues[];

// Register this ObsFunction with the factory
static ObsFunctionMaker<Statistic> makerStatistic_("Statistic");

// -----------------------------------------------------------------------------

Statistic::Statistic(const eckit::LocalConfiguration& conf) : invars_() {
  // Validate and deserialize parameters
  options_.validateAndDeserialize(conf);

  // Validate that weight variable is present for weighted mean
  if (options_.statistic.value() == StatisticMethod::WEIGHTED_MEAN &&
      options_.weightVariable.value() == boost::none) {
    throw eckit::BadParameter(
        "Weight variable must be specified for weighted mean statistic",
        Here());
  }

  // Add main variable to required variables
  invars_ += options_.variable.value();

  // Add weight variable to required variables if needed
  if (options_.weightVariable.value() != boost::none) {
    invars_ += *options_.weightVariable.value();
  }

  oops::Log::debug() << "Statistic: configured" << std::endl;
}

// -----------------------------------------------------------------------------

Statistic::~Statistic() {}

// -----------------------------------------------------------------------------

void Statistic::compute(const ObsFilterData& in,
                        ioda::ObsDataVector<float>& out) const {
  // Check that input variable has float type
  const Variable& var = options_.variable.value();
  if (in.dtype(var) != ioda::ObsDtype::Float) {
    // Probably would be fine for integer types too but not tested. DateTime
    // should be easy (with DateTime hash function) but not implemented yet.
    // Strings could work for mode but not tested either.
    throw eckit::NotImplemented(
        "Statistic ObsFunction only supports float input variables at present",
        Here());
  }

  // Dispatch to appropriate computation method based on statistic type
  switch (options_.statistic.value()) {
    case StatisticMethod::ARITHMETIC_MEAN:
      computeArithmeticMean(in, out);
      break;

    case StatisticMethod::HARMONIC_MEAN:
      computeHarmonicMean(in, out);
      break;

    case StatisticMethod::MEDIAN:
      computeMedian(in, out);
      break;

    case StatisticMethod::MODE:
      computeMode(in, out);
      break;

    case StatisticMethod::WEIGHTED_MEAN:
      computeWeightedMean(in, out);
      break;

    case StatisticMethod::STANDARD_DEVIATION:
      computeStdDev(in, out);
      break;

    case StatisticMethod::VARIANCE:
      computeVariance(in, out);
      break;

    default:
      throw eckit::BadValue("Unsupported statistic method", Here());
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables& Statistic::requiredVariables() const { return invars_; }

// -----------------------------------------------------------------------------

std::vector<float> Statistic::gatherGlobalNonMissingValues(
    const ObsFilterData& in, const ioda::ObsDataVector<float>& data,
    const std::vector<bool>& isOwnedByThisRank, size_t channelIndex) const {
  const float missing = util::missingValue<float>();
  const size_t nlocs = in.nlocs();

  // Collect local non-missing owned values
  std::vector<float> localValues;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (isOwnedByThisRank[iloc] && data[channelIndex][iloc] != missing) {
      localValues.push_back(data[channelIndex][iloc]);
    }
  }

  // Gather values from all ranks
  std::vector<float> globalValues(localValues);
  oops::mpi::allGatherv(in.obsspace().comm(), globalValues);

  return globalValues;
}

// -----------------------------------------------------------------------------

void Statistic::applyStatisticFunction(
    const ObsFilterData& in, ioda::ObsDataVector<float>& out,
    std::function<float(const std::vector<float>&)> computeFunc) const {
  const size_t nlocs = in.nlocs();

  // Get input variable data
  const Variable& var = options_.variable.value();
  ioda::ObsDataVector<float> inputData(in.obsspace(), var.toOopsObsVariables());
  in.get(var, inputData);

  // Get ownership information
  std::vector<bool> isOwnedByThisRank(nlocs);
  in.obsspace().distribution()->patchObs(isOwnedByThisRank);

  // Process each variable
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
    const std::vector<float> globalNonMissingValues =
        gatherGlobalNonMissingValues(in, inputData, isOwnedByThisRank, ivar);

    const float result = computeFunc(globalNonMissingValues);

    // Assign to all locations
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ivar][iloc] = result;
    }
  }
}

// -----------------------------------------------------------------------------

void Statistic::computeArithmeticMean(const ObsFilterData& in,
                                      ioda::ObsDataVector<float>& out) const {
  applyStatisticFunction(in, out, [](const std::vector<float>& globalValues) {
    if (globalValues.empty()) return util::missingValue<float>();

    float sum = 0.0f;
    for (float val : globalValues) {
      sum += val;
    }
    return sum / static_cast<float>(globalValues.size());
  });
}

// -----------------------------------------------------------------------------

void Statistic::computeHarmonicMean(const ObsFilterData& in,
                                    ioda::ObsDataVector<float>& out) const {
  const bool abortOnInvalid = options_.abortIfInvalidOperation.value();
  applyStatisticFunction(in, out, [abortOnInvalid](const std::vector<float>& globalValues) {
    const float missing = util::missingValue<float>();
    if (globalValues.empty()) return missing;

    // Check for zeros or negative values
    for (float val : globalValues) {
      if (val <= 0.0f) {
        if (abortOnInvalid) {
          throw eckit::BadValue(
              "Statistic: harmonic mean encountered zero or negative "
              "value",
              Here());
        } else {
          oops::Log::warning()
              << "Statistic: harmonic mean encountered zero or negative "
              << "values, returning missing value" << std::endl;
          return missing;
        }
      }
    }

    float sumReciprocals = 0.0f;
    for (float val : globalValues) {
      sumReciprocals += 1.0f / val;
    }
    return static_cast<float>(globalValues.size()) / sumReciprocals;
  });
}

// -----------------------------------------------------------------------------

void Statistic::computeMedian(const ObsFilterData& in,
                              ioda::ObsDataVector<float>& out) const {
  applyStatisticFunction(in, out, [](const std::vector<float>& globalValues) {
    if (globalValues.empty()) return util::missingValue<float>();

    // Need to sort, so make a copy
    std::vector<float> sortedValues(globalValues);
    std::sort(sortedValues.begin(), sortedValues.end());

    const size_t n = sortedValues.size();
    if (n % 2 == 1) {
      // Odd number of elements: take middle element
      return sortedValues[n / 2];
    } else {
      // Even number of elements: take average of two middle elements
      return (sortedValues[n / 2 - 1] + sortedValues[n / 2]) / 2.0f;
    }
  });
}

// -----------------------------------------------------------------------------

void Statistic::computeMode(const ObsFilterData& in,
                            ioda::ObsDataVector<float>& out) const {
  applyStatisticFunction(in, out, [](const std::vector<float>& globalValues) {
    const float missing = util::missingValue<float>();
    if (globalValues.empty()) return missing;

    // Count occurrences of each value
    std::unordered_map<float, size_t> valueCounts;
    for (float val : globalValues) {
      valueCounts[val]++;
    }

    // Find maximum count
    size_t maxCount = 0;
    float modeValue = missing;
    for (const auto& pair : valueCounts) {
      if (pair.second > maxCount) {  // new mode candidate, new max count
        maxCount = pair.second;
        modeValue = pair.first;
      } else if (pair.second == maxCount && pair.first < modeValue) {
        modeValue = pair.first;  // tie for mode candidate: keep smallest
      }
    }

    if (maxCount <= 1) return missing;  // No mode exists (all values unique)

    return modeValue;
  });
}

// -----------------------------------------------------------------------------

void Statistic::computeWeightedMean(const ObsFilterData& in,
                                    ioda::ObsDataVector<float>& out) const {
  // NOTE: Since weighted mean requires both value and weight to be non-missing,
  // we cannot use the generic applyStatisticFunction method here.
  const float missing = util::missingValue<float>();
  const size_t nlocs = in.nlocs();

  // Get input variable data
  const Variable& var = options_.variable.value();
  ioda::ObsDataVector<float> inputData(in.obsspace(), var.toOopsObsVariables());
  in.get(var, inputData);

  // Get weight variable data
  const Variable& weightVar = *options_.weightVariable.value();
  ioda::ObsDataVector<float> weightData(in.obsspace(),
                                        weightVar.toOopsObsVariables());
  in.get(weightVar, weightData);

  // Get ownership information
  std::vector<bool> isOwnedByThisRank(nlocs);
  in.obsspace().distribution()->patchObs(isOwnedByThisRank);

  // Process each variable
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
    // Collect local non-missing owned values where BOTH value and weight are
    // non-missing
    std::vector<float> localValues;
    std::vector<float> localWeights;
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (isOwnedByThisRank[iloc] && inputData[ivar][iloc] != missing &&
          weightData[ivar][iloc] != missing) {
        localValues.push_back(inputData[ivar][iloc]);
        localWeights.push_back(weightData[ivar][iloc]);
      }
    }

    // Gather values and weights from all ranks
    std::vector<float> globalValues(localValues);
    std::vector<float> globalWeights(localWeights);
    oops::mpi::allGatherv(in.obsspace().comm(), globalValues);
    oops::mpi::allGatherv(in.obsspace().comm(), globalWeights);

    // Compute weighted mean
    float weightedMean = missing;
    if (!globalValues.empty()) {
      float sumWeightedValues = 0.0f;
      float sumWeights = 0.0f;
      for (size_t i = 0; i < globalValues.size(); ++i) {
        sumWeightedValues += globalValues[i] * globalWeights[i];
        sumWeights += globalWeights[i];
      }
      if (sumWeights > 0.0f) {
        weightedMean = sumWeightedValues / sumWeights;
      }
    }

    // Assign to all locations
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ivar][iloc] = weightedMean;
    }
  }
}

// -----------------------------------------------------------------------------

void Statistic::computeVariance(const ObsFilterData& in,
                                ioda::ObsDataVector<float>& out) const {
  const int ddof = options_.deltaDegreesOfFreedom.value();
  if (ddof < 0) {
    throw eckit::BadParameter(
        "delta degrees of freedom must be non-negative", Here());
  }

  applyStatisticFunction(
      in, out, [ddof](const std::vector<float>& globalValues) {
        const size_t n = globalValues.size();
        if (n <= static_cast<size_t>(ddof)) return util::missingValue<float>();

        // Compute mean
        float sum = 0.0f;
        for (float val : globalValues) {
          sum += val;
        }
        const float mean = sum / static_cast<float>(n);

        // Compute sum of squared differences
        float sumSquaredDiff = 0.0f;
        for (float val : globalValues) {
          const float diff = val - mean;
          sumSquaredDiff += diff * diff;
        }

        // Compute variance with delta degrees of freedom
        return sumSquaredDiff / static_cast<float>(n - ddof);
      });
}

// -----------------------------------------------------------------------------

void Statistic::computeStdDev(const ObsFilterData& in,
                              ioda::ObsDataVector<float>& out) const {
  // Compute variance first (MPI gathering is handled inside computeVariance),
  // then take the element-wise square root below.
  computeVariance(in, out);

  const float missing = util::missingValue<float>();
  const size_t nlocs = in.nlocs();

  // Take square root of variance to get standard deviation
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (out[ivar][iloc] != missing) {
        out[ivar][iloc] = std::sqrt(out[ivar][iloc]);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
