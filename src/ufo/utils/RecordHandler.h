/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_RECORDHANDLER_H_
#define UFO_UTILS_RECORDHANDLER_H_

#include <vector>

#include "ioda/ObsDataVector.h"

#include "ufo/filters/Variables.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
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
  explicit RecordHandler(const ioda::ObsSpace & obsdb,
                         const Variables & filtervars,
                         const ioda::ObsDataVector<int> & flags,
                         const bool retainOnlyIfAllFilterVariablesAreValid);

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

  /// Check each record only contains one value of the category variable. This determines the type
  /// of the category variable and calls the function `checkRecordCategoriesImpl`.
  void checkRecordCategories(const Variable & categoryVariableName) const;

 private:
  /// ObsSpace.
  const ioda::ObsSpace & obsdb_;

  /// Filter variables.
  const Variables & filtervars_;

  /// Filter QC flags.
  const ioda::ObsDataVector<int> & flags_;

  /// Option to choose how to treat observations where there are multiple filter variables.
  /// This quantity is set in the thinning routine that instantiates this class.
  /// If true, an observation is valid only so long as all filter variables have passed QC.
  /// For invalid observation locations (selected by a where clause) any remaining unflagged filter
  /// variables are rejected.
  /// If false, an observation location is valid if any filter variables have not been rejected.
  /// The default value of this parameter is false.
  const bool retainOnlyIfAllFilterVariablesAreValid_;

  /// Obtain 'launch' position associated with each record, which is the location in the record
  /// with the earliest non-missing datetime.
  std::vector<std::size_t> getLaunchPositions() const;

  /// Type-dependent implementation of the check that each record only contains one value of the
  /// category variable.
  template <typename VariableType>
  void checkRecordCategoriesImpl(const Variable & categoryVariableName) const;
};

}  // namespace ufo

#endif  // UFO_UTILS_RECORDHANDLER_H_
