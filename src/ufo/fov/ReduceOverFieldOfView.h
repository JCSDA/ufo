/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cstddef>
#include <vector>

namespace util {
template <typename T> class Range;
}

namespace ufo {
namespace fov {

// -----------------------------------------------------------------------------

/// Compute average of field samples w.r.t. weights.
///
/// Acts on multiple obs, reducing each obs from many samples to one average.
///
/// Expect following vector sizes:
/// - average : resized to nlocs on output; one per obs
/// - sample_ranges : size nlocs; one per obs
/// - sample_values : size sum(sample_nlocs); one per FOV sample over all obs
/// - sample_weights : size sum(sample_nlocs)
void average(std::vector<double> & average,
             const std::vector<util::Range<size_t>> & sample_ranges,
             const std::vector<double> & sample_values,
             const std::vector<double> & sample_weights);

/// Compute average of field samples w.r.t. weights, accounting for masking
///
/// Acts on multiple obs, reducing each obs from many samples to one average.
///
/// Expect following vector sizes:
/// - average : resized to nlocs on output; one per obs
/// - sample_ranges : size nlocs; one per obs
/// - sample_values : size sum(sample_nlocs); one per FOV sample over all obs
/// - sample_weights : size sum(sample_nlocs)
/// - mask : size nlocs; where this is 0, average will be set to valueOutsideMask
/// - sample_mask : size sum(sample_nlocs)
void averageWithMask(std::vector<double> & average,
                     const std::vector<util::Range<size_t>> & sample_ranges,
                     const std::vector<double> & sample_values,
                     const std::vector<double> & sample_weights,
                     const std::vector<double> & mask,
                     const std::vector<double> & sample_mask,
                     double valueOutsideMask);

/// Compute dominant integer value in field samples w.r.t. weights, accounting for masking
///
/// Acts on multiple obs, reducing each obs from many samples to one average.
///
/// Expect following vector sizes:
/// - dominant_int : resized to nlocs on output; one per obs; should hold int values
/// - sample_ranges : size nlocs; one per obs
/// - sample_int_values : size sum(sample_nlocs); one per FOV sample over all obs; int values
/// - sample_weights : size sum(sample_nlocs)
/// - mask : size nlocs; where this is 0, average will be set to valueOutsideMask
/// - sample_mask : size sum(sample_nlocs)
void dominantIntegerValue(std::vector<double> & dominant_int,
                          const std::vector<util::Range<size_t>> & sample_ranges,
                          const std::vector<double> & sample_int_values,
                          const std::vector<double> & sample_weights,
                          const std::vector<double> & mask,
                          const std::vector<double> & sample_mask,
                          int valueOutsideMask);

// -----------------------------------------------------------------------------

}  // namespace fov
}  // namespace ufo
