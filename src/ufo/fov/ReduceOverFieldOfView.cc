/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/fov/ReduceOverFieldOfView.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "oops/util/Range.h"

#include "ufo/fov/SampleFieldOfView.h"

namespace ufo {
namespace fov {

// -----------------------------------------------------------------------------

void average(std::vector<double> & average,
             const std::vector<util::Range<size_t>> & sample_ranges,
             const std::vector<double> & sample_values,
             const std::vector<double> & sample_weights) {
  const size_t nlocs = sample_ranges.size();
  // The statement below assumes the samples run over obs locs in increasing order
  const size_t nsamples = sample_ranges.back().end;
  ASSERT(sample_values.size() == nsamples);
  ASSERT(sample_weights.size() == nsamples);

  if (average.size() != nlocs) {
    average.resize(nlocs);
  }

  for (size_t ifov = 0; ifov < nlocs; ++ifov) {
    average[ifov] = 0.0;
    double norm = 0.0;
    for (size_t i = sample_ranges[ifov].begin; i < sample_ranges[ifov].end; ++i) {
      average[ifov] += sample_values[i] * sample_weights[i];
      norm += sample_weights[i];
    }
    ASSERT(norm > 0.0);
    average[ifov] /= norm;
  }
}

// -----------------------------------------------------------------------------

void averageWithMask(std::vector<double> & average,
                     const std::vector<util::Range<size_t>> & sample_ranges,
                     const std::vector<double> & sample_values,
                     const std::vector<double> & sample_weights,
                     const std::vector<double> & mask,
                     const std::vector<double> & sample_mask,
                     const double valueOutsideMask) {
  const size_t nlocs = sample_ranges.size();
  // The statement below assumes the samples run over obs locs in increasing order
  const size_t nsamples = sample_ranges.back().end;
  ASSERT(sample_values.size() == nsamples);
  ASSERT(sample_weights.size() == nsamples);
  ASSERT(mask.size() == nlocs);
  ASSERT(sample_mask.size() == nsamples);

  if (average.size() != nlocs) {
    average.resize(nlocs);
  }

  for (size_t ifov = 0; ifov < nlocs; ++ifov) {
    if (mask[ifov] > 0.0) {
      // compute masked average
      average[ifov] = 0.0;
      double norm = 0.0;
      for (size_t i = sample_ranges[ifov].begin; i < sample_ranges[ifov].end; ++i) {
        average[ifov] += sample_values[i] * sample_weights[i] * sample_mask[i];
        norm += sample_weights[i] * sample_mask[i];
      }
      ASSERT(norm > 0.0);
      average[ifov] /= norm;
    } else {
      // outside mask; fall back to alternative value
      average[ifov] = valueOutsideMask;
    }
  }
}

// -----------------------------------------------------------------------------

void dominantIntegerValue(std::vector<double> & dominant_int,
                     const std::vector<util::Range<size_t>> & sample_ranges,
                     const std::vector<double> & sample_int_values,
                     const std::vector<double> & sample_weights,
                     const std::vector<double> & mask,
                     const std::vector<double> & sample_mask,
                     const int valueOutsideMask) {
  const size_t nlocs = sample_ranges.size();
  // The statement below assumes the samples run over obs locs in increasing order
  const size_t nsamples = sample_ranges.back().end;
  ASSERT(sample_int_values.size() == nsamples);
  ASSERT(sample_weights.size() == nsamples);
  ASSERT(mask.size() == nlocs);
  ASSERT(sample_mask.size() == nsamples);

  if (dominant_int.size() != nlocs) {
    dominant_int.resize(nlocs);
  }

  int minval = std::numeric_limits<int>().max();
  int maxval = std::numeric_limits<int>().min();
  for (size_t i = 0; i < nsamples; ++i) {
    // Exclude missingValues that arise outside interpolation masks
    if (sample_mask[i] > 0.0) {
      minval = std::min(minval, static_cast<int>(std::round(sample_int_values[i])));
      maxval = std::max(maxval, static_cast<int>(std::round(sample_int_values[i])));
    }
  }

  std::vector<double> int_weights(maxval - minval + 1);

  for (size_t ifov = 0; ifov < nlocs; ++ifov) {
    if (mask[ifov] > 0.0) {
      std::fill(int_weights.begin(), int_weights.end(), 0.0);
      for (size_t i = sample_ranges[ifov].begin; i < sample_ranges[ifov].end; ++i) {
        if (sample_mask[i] > 0.0) {
          const int this_int = static_cast<int>(std::round(sample_int_values[i]));
          int_weights[this_int - minval] += sample_weights[i] * sample_mask[i];
        }
      }
      dominant_int[ifov] = minval + std::distance(int_weights.begin(),
                                    std::max_element(int_weights.begin(), int_weights.end()));
    } else {
      // outside mask; fall back to alternative value
      dominant_int[ifov] = valueOutsideMask;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace fov
}  // namespace ufo
