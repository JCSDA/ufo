/*
 * (C) Crown Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_STEPCHECKPARAMETERS_H_
#define UFO_FILTERS_STEPCHECKPARAMETERS_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/TrackCheckUtilsParameters.h"

namespace ufo {

/// \brief Options controlling the operation of step check filter.

/// Inherits from TrackCheckUtilsParameters to gain access to the
/// station_id_variable parameter, which can control the ObsAccessor grouping
/// behaviour if an obs grouping has not been specified in the ObsSpace.
class StepCheckParameters : public TrackCheckUtilsParameters {
  OOPS_CONCRETE_PARAMETERS(StepCheckParameters, TrackCheckUtilsParameters)

 public:
/// The threshold for step magnitude between observations to trigger flagging
/// at a station. Flags are applied to the observation where the step occurs (the
/// second observation in the pair) unless chunking is used (see below).
  oops::RequiredParameter<float> stepThreshold{"step threshold", this};

  /// Whether the step threshold is inclusive (>=) or exclusive (>).
  /// Default is true (inclusive).
  oops::Parameter<bool> inclusiveStepThreshold{"inclusive step threshold", true,
                                               this};

  /// The number of steps that are allowed before a timeseries is flagged.
  /// Not compatible with "percentage step tolerance" or "use average step".
  oops::OptionalParameter<size_t> numberStepTolerance{"number step tolerance",
                                                      this};

  /// The percentage of steps (out of total valid steps) that are allowed before
  /// flagging begins.
  /// Not compatible with "number step tolerance" or "use average step".
  oops::OptionalParameter<float> percentageStepTolerance{
      "percentage step tolerance", this};

  /// If specified, enables circular difference calculation for periodic data.
  /// The value represents the period (e.g., 360.0 for degrees).
  oops::OptionalParameter<float> circularPeriod{"circular period", this};

  /// If specified, observations are split into chunks of this size and the mean
  /// step within each chunk is calculated (using circular difference if
  /// specified). These per-chunk mean steps are then averaged across all
  /// non-stuck chunks, and that overall average is compared to the threshold.
  /// If exceeded, all observations in the record are flagged.
  /// Requires "use average step" to be true.
  oops::OptionalParameter<size_t> chunkSize{"chunk size", this};

  /// Whether to ignore the last chunk if it contains fewer observations than
  /// chunk size.
  /// Default is false.
  oops::Parameter<bool> ignoreLastChunkIfIncomplete{
      "ignore last chunk if incomplete", false, this};

  /// Whether to remove chunks where all values are the same (within "chunk
  /// stuck tolerance"). Stuck chunks are excluded before calculating the
  /// average of the per-chunk mean steps. This avoids sets of identical values
  /// artificially deflating the overall average step magnitude.
  /// Default is false.
  oops::Parameter<bool> removeStuckChunks{"remove stuck chunks", false, this};

  /// Tolerance for determining if values in a chunk are the same.
  /// Used only if remove stuck chunks is true.
  /// Default is 0.0.
  oops::Parameter<float> chunkStuckTolerance{"chunk stuck tolerance", 0.0,
                                             this};

  /// Whether to calculate the average of all steps and compare to threshold.
  /// If the average step exceeds the threshold, all observations in the record
  /// are flagged. Required when using "chunk size".
  /// Not compatible with "number step tolerance" or "percentage step
  /// tolerance". Default is false.
  oops::Parameter<bool> useAverageStep{"use average step", false, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_STEPCHECKPARAMETERS_H_
