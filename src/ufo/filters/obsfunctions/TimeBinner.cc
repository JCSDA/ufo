/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/TimeBinner.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionTimeBinner> makerObsFunctionTimeBinner_(
    "TimeBinner");

ObsFunctionTimeBinner::ObsFunctionTimeBinner(
    const eckit::LocalConfiguration &conf) {
  oops::Log::trace() << "ObsFunctionTimeBinner constructor" << std::endl;
  options_.validateAndDeserialize(conf);

  // Add the input timestamp variable to the required variables
  invars_ += options_.inputTimestampVariable.value();

  if (options_.binnedTimestampVariable.value() != boost::none) {
    const boost::optional<Variable> &binnedVar =
        options_.binnedTimestampVariable.value();
    invars_ += *binnedVar;
  }

  // Validate that the bin window does not overlap with adjacent bins
  const util::Duration binIntervalUnit =
      getBinIntervalUnit(options_.binIntervalUnit.value());
  const util::Duration binWindowLowerBound{
      options_.binWindowLowerBound.value()};
  const util::Duration binWindowUpperBound{
      options_.binWindowUpperBound.value()};

  if (binWindowLowerBound > binWindowUpperBound) {
    throw eckit::BadValue(
        "The bin window lower bound must be less than or equal to the bin "
        "window upper bound.",
        Here());
  }
  if ((binWindowUpperBound - binWindowLowerBound).toSeconds() >=
      binIntervalUnit.toSeconds()) {
    throw eckit::BadValue(
        "Bin window (lower bound to upper bound) must be less than the bin "
        "interval unit. "
        "This ensures no overlap with adjacent bins.",
        Here());
  }
  // neither the bin window upper bound or lower bound can have a magnitude of
  // greater than the bin interval unit, otherwise the flooring process to
  // derive the bin label would be incorrect.
  if (std::abs(binWindowLowerBound.toSeconds()) > binIntervalUnit.toSeconds() ||
      std::abs(binWindowUpperBound.toSeconds()) > binIntervalUnit.toSeconds()) {
    throw eckit::BadValue(
        "Bin window lower and upper bounds must not exceed the bin interval "
        "unit in magnitude. "
        "This ensures correct bin label derivation.",
        Here());
  }
}

ObsFunctionTimeBinner::~ObsFunctionTimeBinner() {
  oops::Log::trace() << "ObsFunctionTimeBinner destructor" << std::endl;
}

void ObsFunctionTimeBinner::compute(const ObsFilterData &in,
                                    ioda::ObsDataVector<int> &out) const {
  oops::Log::trace() << "ObsFunctionTimeBinner compute start" << std::endl;

  const int missingBinValue =
      util::missingValue<int>();  // todo change  to missingValue<T> when this
                                  // is allowed to be not just an IntObsFunction

  const size_t nlocs = in.nlocs();

  // Retrieve input timestamps
  std::vector<util::DateTime> timestamps(nlocs);
  in.get(options_.inputTimestampVariable.value(), timestamps);

  // Parse bin interval and window
  const util::Duration binIntervalUnit =
      getBinIntervalUnit(options_.binIntervalUnit.value());
  const int64_t binIntervalUnitSecs = binIntervalUnit.toSeconds();
  const util::Duration binWindowLowerBound{
      options_.binWindowLowerBound.value()};
  const util::Duration binWindowUpperBound{
      options_.binWindowUpperBound.value()};
  const bool reverseOrder = options_.reverseChronologicalOrder.value();

  // A fixed epoch is used for integer-second arithmetic
  const util::DateTime epoch{"1970-01-01T00:00:00Z"};

  // Bins are labelled with timestamps which are derived from a "floor" derived
  // from the bin interval unit. Consider the timestamp 2018-04-17 01:20:02 UTC:
  // if the bin interval unit is "hour" the floored time is 2018-04-17 01:00:00
  // UTC, if the bin interval unit is "minute" the floored time is 2018-04-17
  // 01:20:00 UTC, if the bin interval unit is "day" the floored time is
  // 2018-04-17 00:00:00 UTC, etc.
  //
  // Which "floor"ed timestamp label is applied to a given timestamp depends on
  // the bin window start and end. For example, with a bin interval unit of
  // "hour" and a bin window start of -25 minutes and end of +25 minutes, the
  // timestamp 2018-04-17 01:20:02 UTC would be labelled with the timestamp
  // 2018-04-17 01:00:00 UTC (since it falls within the window of that bin), as
  // would the timestamp 2018-04-17 00:50:00 UTC. All timestamps in that window
  // would be assigned to the same bin number.
  //
  // Note that reverseOrder only affects the numbering of the bins, not the
  // timestamps assigned to them.
  //
  // Bin windows do not have to cover the label timestamp. Bin window starts
  // can be positive, implying that the window starts after the bin label
  // timestamp. Bin window ends can also be negative, implying that the window
  // ends before the bin label timestamp. However, the bin window (end - start)
  // must not exceed the bin interval unit, to avoid overlap with adjacent bins
  // and the start must be less than or equal to the end.

  // Iteratively find a valid reference timestamp that falls within a bin window
  const util::DateTime firstBinLabelTimestamp = findFirstBinLabelTimestamp(
      in.obsspace(), timestamps, binWindowLowerBound, binWindowUpperBound,
      binIntervalUnit, epoch, reverseOrder);

  // Compute the bin numbers and optionally the binned timestamps
  std::vector<util::DateTime> binLabelTimestamps(
      nlocs, util::missingValue<util::DateTime>());
  if (firstBinLabelTimestamp == util::missingValue<util::DateTime>()) {
    // No valid reference timestamp found, so all bin outputs are missing
    std::fill(out[0].begin(), out[0].end(), missingBinValue);
  } else {
    // Assigns output bin index (respecting reverse ordering) and label
    // timestamp.
    const auto assignBin = [&](size_t iloc, int64_t bin,
                               const util::DateTime &binLabelTimestamp) {
      out[0][iloc] =
          reverseOrder ? static_cast<int>(-bin) : static_cast<int>(bin);
      binLabelTimestamps[iloc] = binLabelTimestamp;
    };
    // Returns true if the timestamp falls within the given bin label's window.
    const auto inWindow = [&](const util::DateTime &timestamp,
                              const util::DateTime &binLabelTimestamp) {
      return timestamp >= binLabelTimestamp + binWindowLowerBound &&
             timestamp <= binLabelTimestamp + binWindowUpperBound;
    };

    const int64_t firstBinLabelTimestampSecs =
        (firstBinLabelTimestamp - epoch).toSeconds();
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (timestamps[iloc] == util::missingValue<util::DateTime>()) {
        out[0][iloc] = missingBinValue;
        continue;
      }
      // Find the candidate bin label timestamp for the current timestamp by
      // flooring to the nearest bin interval unit. This uses modular arithmetic:
      // timestampSecs - (timestampSecs % binIntervalUnitSecs) floors to the
      // nearest multiple of binIntervalUnitSecs. Example: for a timestamp of
      // 2018-04-17 01:20:02 UTC with hourly bins (binIntervalUnitSecs=3600):
      // timestampSecs=1523928002, 1523928002 % 3600 = 1202, so the floor is
      // 1523926800 which corresponds to 2018-04-17 01:00:00 UTC.
      const int64_t timestampSecs = (timestamps[iloc] - epoch).toSeconds();
      const int64_t candidateBinLabelTimestampSecs =
          timestampSecs - (timestampSecs % binIntervalUnitSecs);
      const util::DateTime candidateBinLabelTimestamp =
          epoch + util::Duration(candidateBinLabelTimestampSecs);

      // The candidate bin label should be exactly divisible by the bin interval
      // unit, allowing the probable bin number to be deduced
      assert((candidateBinLabelTimestampSecs % binIntervalUnitSecs) == 0);
      int64_t candidateBinNumber =
          (candidateBinLabelTimestampSecs - firstBinLabelTimestampSecs) /
          binIntervalUnitSecs;

      if (inWindow(timestamps[iloc], candidateBinLabelTimestamp)) {
        assignBin(iloc, candidateBinNumber, candidateBinLabelTimestamp);
        continue;
      }
      // The timestamp doesn't fall within the candidate bin's window, but may
      // still fall within the window of a neighboring bin, so check those
      // before giving up and assigning missing. Note that only the
      // immediately adjacent bins need to be checked, due to the validation
      // that the bin window starts and ends cannot exceed the bin interval
      // unit. Note also we don't check the previous bin if the candidate bin
      // number is 0, since we don't allow negative bin numbers in the current
      // implementation (see the documentation).
      out[0][iloc] = missingBinValue;
      const util::DateTime nextBinLabelTimestamp =
          candidateBinLabelTimestamp + binIntervalUnit;
      if (inWindow(timestamps[iloc], nextBinLabelTimestamp)) {
        assignBin(iloc, candidateBinNumber + 1, nextBinLabelTimestamp);
      } else if (candidateBinNumber > 0) {
        const util::DateTime prevBinLabelTimestamp =
            candidateBinLabelTimestamp - binIntervalUnit;
        if (inWindow(timestamps[iloc], prevBinLabelTimestamp)) {
          assignBin(iloc, candidateBinNumber - 1, prevBinLabelTimestamp);
        }
      }
    }
  }
  // Save bin label timestamps if requested
  if (options_.binnedTimestampVariable.value() != boost::none) {
    const Variable &binnedVar = *options_.binnedTimestampVariable.value();
    in.obsspace().put_db(binnedVar.group(), binnedVar.variable(),
                         binLabelTimestamps);
  }

  oops::Log::trace() << "ObsFunctionTimeBinner compute complete" << std::endl;
}

const util::Duration ObsFunctionTimeBinner::getBinIntervalUnit(
    const std::string &unit) const {
  if (unit == "second" || unit == "seconds") {
    return util::Duration("PT1S");
  } else if (unit == "minute" || unit == "minutes") {
    return util::Duration("PT1M");
  } else if (unit == "hour" || unit == "hours") {
    return util::Duration("PT1H");
  } else if (unit == "day" || unit == "days") {
    return util::Duration("P1D");
  } else {
    throw eckit::BadValue("Unsupported bin interval unit: " + unit, Here());
  }
}

const ufo::Variables &ObsFunctionTimeBinner::requiredVariables() const {
  return invars_;
}

const util::DateTime ObsFunctionTimeBinner::findFirstBinLabelTimestamp(
    const ioda::ObsSpace &obsSpace,
    const std::vector<util::DateTime> &timestamps,
    const util::Duration &binWindowLowerBound,
    const util::Duration &binWindowUpperBound,
    const util::Duration &binIntervalUnit, const util::DateTime &epoch,
    bool reverseOrder) const {
  // Find a valid reference timestamp (min or max, depending on order) that
  // falls within a bin window. This determines the first bin label and allows
  // all other bin numbers to be computed consistently. Iteratively try
  // reference timestamps until one is found that fits, or return missing if
  // none fit.
  auto workingTimestamps = timestamps;
  const int64_t binWindowLowerBoundSecs = binWindowLowerBound.toSeconds();
  const int64_t binIntervalUnitSecs = binIntervalUnit.toSeconds();

  while (!workingTimestamps.empty()) {
    // Select the earliest (chronological) or latest (reverse) timestamp.
    const util::DateTime trialReferenceTimestamp =
        reverseOrder
            ? getGlobalMaxTimestamp(obsSpace, workingTimestamps, epoch)
            : getGlobalMinTimestamp(obsSpace, workingTimestamps, epoch);
    const int64_t trialReferenceTimestampSecs =
        (trialReferenceTimestamp - epoch).toSeconds();

    // Compute the candidate bin label for this reference timestamp. The label
    // is derived by flooring to the bin interval unit, but first we adjust the
    // timestamp to account for the window lower bound offset. This ensures that
    // timestamps near (but after) a potential bin label are mapped to the
    // correct label, not the previous one. For example: with a 1-hour interval
    // and a relative window of [-30 min, +30 min), a reference at 19:40
    // belongs to the 20:00 bin window [19:30, 20:30), not the 19:00 window.
    // The adjustment below handles this by shifting negative window lower
    // bounds forward before flooring.
    int64_t trialReferenceTimestampAdjustedForWindowSecs =
        trialReferenceTimestampSecs;
    if (binWindowLowerBoundSecs < 0) {
      trialReferenceTimestampAdjustedForWindowSecs -= binWindowLowerBoundSecs;
    }

    const int64_t candidateBinLabelTimestampSecs =
        trialReferenceTimestampAdjustedForWindowSecs -
        (trialReferenceTimestampAdjustedForWindowSecs % binIntervalUnitSecs);
    const util::DateTime candidateBinLabelTimestamp =
        epoch + util::Duration(candidateBinLabelTimestampSecs);

    // Check if the reference timestamp falls within this candidate bin's
    // window.
    const util::DateTime candidateBinWindowLowerBound =
        candidateBinLabelTimestamp + binWindowLowerBound;
    const util::DateTime candidateBinWindowUpperBound =
        candidateBinLabelTimestamp + binWindowUpperBound;

    if (trialReferenceTimestamp >= candidateBinWindowLowerBound &&
        trialReferenceTimestamp <= candidateBinWindowUpperBound) {
      // Found a valid first bin label.
      return candidateBinLabelTimestamp;
    }

    // This reference timestamp doesn't fit any bin window; remove it and
    // iterate again using the next min/max timestamp.
    workingTimestamps.erase(
        std::remove(workingTimestamps.begin(), workingTimestamps.end(),
                    trialReferenceTimestamp),
        workingTimestamps.end());
  }
  return util::missingValue<util::DateTime>();
}

const util::DateTime ObsFunctionTimeBinner::getGlobalMinTimestamp(
    const ioda::ObsSpace &obsSpace,
    const std::vector<util::DateTime> &timestamps,
    const util::DateTime &epoch) const {
  const util::DateTime minTimestamp = *std::min_element(
      timestamps.begin(), timestamps.end(),
      [](const util::DateTime &a, const util::DateTime &b) { return a < b; });
  int64_t minTimestampSecs = (minTimestamp - epoch).toSeconds();
  obsSpace.comm().allReduceInPlace(minTimestampSecs, eckit::mpi::min());
  return epoch + util::Duration(minTimestampSecs);
}

const util::DateTime ObsFunctionTimeBinner::getGlobalMaxTimestamp(
    const ioda::ObsSpace &obsSpace,
    const std::vector<util::DateTime> &timestamps,
    const util::DateTime &epoch) const {
  const util::DateTime maxTimestamp = *std::max_element(
      timestamps.begin(), timestamps.end(),
      [](const util::DateTime &a, const util::DateTime &b) { return a < b; });
  int64_t maxTimestampSecs = (maxTimestamp - epoch).toSeconds();
  obsSpace.comm().allReduceInPlace(maxTimestampSecs, eckit::mpi::max());
  return epoch + util::Duration(maxTimestampSecs);
}

}  // namespace ufo
