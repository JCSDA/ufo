/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_RECORDHANDLER_H_
#define UFO_UTILS_RECORDHANDLER_H_

#include <vector>

namespace ioda {
  class ObsSpace;
}

namespace ufo {

  /// \brief Class which is used to ensure records are treated as single observations in the spatial
  /// and temporal thinning routines.
  /// \details If a data set has been divided into records, it is sometimes desirable to treat each
  /// record as a single observation in the Gaussian Thinning and Temporal Thinning filters.
  /// For example, this can be used to match the behaviour of OPS.
  /// The routines in the `RecordHandler` class can be used to modify the `apply` and `isThinned`
  /// vectors that are used in both of the aforementioned filters.
  /// The modifications are described in the documentation of each routine.
class RecordHandler
{
 public:
  explicit RecordHandler(const ioda::ObsSpace & obsdb);

  /// Modify the input `apply` vector if records are treated as single observations.
  /// This function first finds the location of the earliest non-missing datetime in each record,
  /// corresponding to the spatial and temporal location of the launch site.
  /// The value of `apply` at this location is preserved (as either true or false), and all other
  /// values of `apply` in the record are set to false.
  std::vector<bool> changeApplyIfRecordsAreSingleObs(const std::vector<bool> &) const;

  /// Modify the input `isThinned` vector if records are treated as single observations.
  /// This function first finds the location of the earliest non-missing datetime in each record,
  /// corresponding to the spatial and temporal location of the launch site.
  /// All values of `isThinned` in the record are set to the the value of `isThinned` at the
  /// location of the launch site.
  std::vector<bool> changeThinnedIfRecordsAreSingleObs(const std::vector<bool> &) const;

 private:
  /// ObsSpace.
  const ioda::ObsSpace & obsdb_;

  /// Obtain 'launch' position associated with each record, which is the location in the record
  /// with the earliest non-missing datetime.
  std::vector<std::size_t> getLaunchPositions() const;
};

}  // namespace ufo

#endif  // UFO_UTILS_RECORDHANDLER_H_
