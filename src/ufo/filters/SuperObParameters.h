/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_SUPEROBPARAMETERS_H_
#define UFO_FILTERS_SUPEROBPARAMETERS_H_

#include <vector>
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "ufo/filters/FilterParametersBase.h"

namespace ufo {

// Forward declarations (full definitions elsewhere)
class SuperObParametersBase;
class SuperObFactory;

class SuperObParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SuperObParametersWrapper, oops::Parameters)
 public:
  /// Name of the superobbing algorithm.
  /// Valid names are specified using a `SuperObMaker` in subclasses of
  /// SuperObBase in the ufo/superob directory.
  oops::RequiredPolymorphicParameter<SuperObParametersBase, SuperObFactory>
      superObName{"name", this};
};

class SuperObParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SuperObParameters, FilterParametersBase)

 public:
  /// Parameter that contains details of the algorithm to use.
  oops::RequiredParameter<SuperObParametersWrapper> algorithmParameters{
      "algorithm", this};

  /// If `true` (default), when a where clause is used, locations selected by
  /// the where clause are set to missing before superobs are assigned.
  /// If `false`, existing values are preserved outside the where selection.
  oops::OptionalParameter<std::vector<bool>>
      setValuesOutsideWhereClauseToMissing{
          "set values outside where clause to missing", this};

  /// If `increment if non-missing` is true, a variable specified by
  /// `variables to increment` is incremented by the corresponding
  /// `increment values` entry for each superob computed. Variables to increment
  /// must be integer valued and already exist in the ObsSpace. If
  /// `increment whole record` is true, the increment is applied to all
  /// locations in the record selected by the where clause , otherwise only to
  /// the location(s) where the superob is stored (this is false by default).
  /// If `increment whole record respects where` is false (default is true) then
  /// the increment is applied to all locations in the record, ignoring the
  /// `where` clause.
  oops::OptionalParameter<std::vector<bool>> incrementIfNonMissing{
      "increment if non-missing", this};
  oops::OptionalParameter<std::vector<Variable>> variablesToIncrement{
      "variables to increment", this};
  oops::OptionalParameter<std::vector<int>> incrementValues{"increment values",
                                                            this};
  oops::OptionalParameter<std::vector<bool>> incrementWholeRecord{
      "increment whole record", this};
  oops::OptionalParameter<std::vector<bool>> incrementWholeRecordRespectsWhere{
      "increment whole record respects where", this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_SUPEROBPARAMETERS_H_
