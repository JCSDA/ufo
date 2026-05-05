/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_TIMEBINNER_H_
#define UFO_FILTERS_OBSFUNCTIONS_TIMEBINNER_H_

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Parameters controlling IntObsFunction/TimeBinner.
///
/// TimeBinner assigns each observation timestamp to a bin label timestamp and
/// returns a bin number as an integer ObsFunction output.
///
/// Bin membership rule:
/// An observation timestamp t belongs to bin label T when
///   T + bin window lower bound <= t <= T + bin window upper bound
///
/// Key constraints (validated in the implementation):
/// - bin window lower bound <= bin window upper bound
/// - (bin window upper bound - bin window lower bound) < bin interval unit
/// - |bin window lower bound| and |bin window upper bound| are each <= bin
/// interval unit
///
/// If an observation falls outside all bin windows, the output bin number is
/// missing.
class TimeBinnerParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TimeBinnerParameters, Parameters)

 public:
  /// \brief Bin interval granularity.
  ///
  /// Accepted values: second, seconds, minute, minutes, hour, hours, day, days.
  oops::RequiredParameter<std::string> binIntervalUnit{"bin interval unit",
                                                       this};

  /// \brief Lower bound of the bin window relative to the bin label timestamp.
  ///
  /// Example values: -PT30M, PT0S, PT1S.
  oops::RequiredParameter<std::string> binWindowLowerBound{
      "bin window lower bound", this};

  /// \brief Upper bound of the bin window relative to the bin label timestamp.
  ///
  /// Example values: PT0M, PT30M, PT1H.
  oops::RequiredParameter<std::string> binWindowUpperBound{
      "bin window upper bound", this};

  /// \brief Input timestamp variable.
  ///
  /// Default: MetaData/dateTime.
  oops::Parameter<Variable> inputTimestampVariable{
      "input timestamp variable", Variable{"MetaData/dateTime"}, this};

  /// \brief Optional output variable holding bin label timestamps.
  ///
  /// For observations assigned to bins, the label timestamp is written.
  /// For observations outside all bin windows, missing DateTime is written.
  oops::OptionalParameter<Variable> binnedTimestampVariable{
      "binned timestamp variable", this};

  /// \brief If true, bin numbering is reverse chronological.
  ///
  /// - false: bin 0 is the earliest valid bin label.
  /// - true:  bin 0 is the latest valid bin label.
  oops::Parameter<bool> reverseChronologicalOrder{"reverse chronological order",
                                                  false, this};
};

// -----------------------------------------------------------------------------

/// \brief Bins observation timestamps into discrete time intervals.
///
/// This is an IntObsFunction registered as IntObsFunction/TimeBinner.
///
/// The output integer is the bin number. Optionally, a second variable can be
/// populated with the associated bin label timestamp.
///
/// Example A: window ending at the hour boundary
/// \code{.yaml}
/// - filter: Variable Assignment
///   assignments:
///     - name: MetaData/binNumbers
///       type: int
///       function:
///         name: IntObsFunction/TimeBinner
///         options:
///           bin interval unit: hour
///           bin window lower bound: -PT30M
///           bin window upper bound: PT0M
///           binned timestamp variable: MetaData/binTimestamps
/// \endcode
/// | input timestamp      | bin number | binned timestamp     |
/// |----------------------|------------|----------------------|
/// | 2018-04-16T23:50:00Z | 0          | 2018-04-17T00:00:00Z |
/// | 2018-04-17T00:15:00Z | missing    | missing              |
/// | 2018-04-17T00:35:00Z | 1          | 2018-04-17T01:00:00Z |
///
/// Example B: window spanning the hour boundary
/// \code{.yaml}
/// - filter: Variable Assignment
///   assignments:
///     - name: MetaData/binNumbers
///       type: int
///       function:
///         name: IntObsFunction/TimeBinner
///         options:
///           bin interval unit: hour
///           bin window lower bound: -PT30M
///           bin window upper bound: PT29M59S
///           binned timestamp variable: MetaData/binTimestamps
/// \endcode
/// | input timestamp      | bin number | binned timestamp     |
/// |----------------------|------------|----------------------|
/// | 2018-04-17T00:15:00Z | 0          | 2018-04-17T00:00:00Z |
/// | 2018-04-17T00:35:00Z | 1          | 2018-04-17T01:00:00Z |
///
/// Example C: reverse chronological order (same window as Example B)
/// \code{.yaml}
/// - filter: Variable Assignment
///   assignments:
///     - name: MetaData/binNumbers
///       type: int
///       function:
///         name: IntObsFunction/TimeBinner
///         options:
///           bin interval unit: hour
///           bin window lower bound: -PT30M
///           bin window upper bound: PT29M59S
///           binned timestamp variable: MetaData/binTimestamps
///           reverse chronological order: true
/// \endcode
/// | input timestamp      | bin number | binned timestamp     |
/// |----------------------|------------|----------------------|
/// | 2018-04-17T00:15:00Z | 1          | 2018-04-17T00:00:00Z |
/// | 2018-04-17T00:35:00Z | 0          | 2018-04-17T01:00:00Z |
class ObsFunctionTimeBinner : public ObsFunctionBase<int> {
 public:
  /// \brief Construct TimeBinner from YAML/local configuration.
  ///
  /// Validates TimeBinner options and declares required input variables.
  ///
  /// \param conf Local configuration containing TimeBinner options.
  explicit ObsFunctionTimeBinner(const eckit::LocalConfiguration &conf);

  /// \brief Destroy the TimeBinner instance.
  ~ObsFunctionTimeBinner();

  /// \brief Compute bin numbers for all observations.
  ///
  /// Writes integer bin numbers to \p out. If configured, also writes the
  /// corresponding bin label timestamps to ObsSpace. Observations with missing
  /// input timestamps receive missing (invalid) bin numbers.
  ///
  /// \param in Input data accessor providing observation values.
  /// \param out Output integer vector storing bin numbers.
  void compute(const ObsFilterData &in,
               ioda::ObsDataVector<int> &out) const override;

  /// \brief Return variables required by this ObsFunction.
  ///
  /// Includes the input timestamp variable and optionally the configured output
  /// binned timestamp variable.
  ///
  /// \return Set of variables required by TimeBinner.
  const ufo::Variables &requiredVariables() const override;

 private:
  TimeBinnerParameters options_;
  ufo::Variables invars_;

  /// \brief Find the first valid bin label timestamp used as bin-number origin.
  ///
  /// Iteratively selects a trial reference timestamp (min or max depending on
  /// \p reverseOrder), derives a candidate bin label, and returns the first
  /// candidate whose window contains that trial timestamp. If no timestamps fit
  /// into any bin window, returns missing DateTime.
  ///
  /// \param obsSpace ObsSpace used for MPI-aware min/max operations.
  /// \param timestamps Candidate timestamps (must not be empty).
  /// \param binWindowLowerBound Window lower bound offset from bin label.
  /// \param binWindowUpperBound Window upper bound offset from bin label.
  /// \param binIntervalUnit Bin interval unit duration (hour, minute, second,
  /// day).
  /// \param epoch Epoch used for integer-second arithmetic.
  /// \param reverseOrder If true, select global max first; otherwise min first.
  /// \return First valid bin label timestamp, or missing DateTime if none
  /// found.
  const util::DateTime findFirstBinLabelTimestamp(
      const ioda::ObsSpace &obsSpace,
      const std::vector<util::DateTime> &timestamps,
      const util::Duration &binWindowLowerBound,
      const util::Duration &binWindowUpperBound,
      const util::Duration &binIntervalUnit, const util::DateTime &epoch,
      bool reverseOrder) const;

  /// \brief Return global minimum timestamp across MPI ranks.
  ///
  /// \param obsSpace ObsSpace providing communicator.
  /// \param timestamps Local timestamps (must not be empty).
  /// \param epoch Epoch used for integer-second arithmetic.
  /// \return Global minimum timestamp.
  const util::DateTime getGlobalMinTimestamp(
      const ioda::ObsSpace &obsSpace,
      const std::vector<util::DateTime> &timestamps,
      const util::DateTime &epoch) const;

  /// \brief Return global maximum timestamp across MPI ranks.
  ///
  /// \param obsSpace ObsSpace providing communicator.
  /// \param timestamps Local timestamps (must not be empty).
  /// \param epoch Epoch used for integer-second arithmetic.
  /// \return Global maximum timestamp.
  const util::DateTime getGlobalMaxTimestamp(
      const ioda::ObsSpace &obsSpace,
      const std::vector<util::DateTime> &timestamps,
      const util::DateTime &epoch) const;

  /// \brief Convert bin interval unit string to util::Duration.
  ///
  /// Supported units are second(s), minute(s), hour(s), and day(s).
  ///
  /// \param unit Unit string from configuration.
  /// \return Duration corresponding to \p unit.
  /// \throw eckit::BadValue if \p unit is unsupported.
  const util::Duration getBinIntervalUnit(const std::string &unit) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_TIMEBINNER_H_
