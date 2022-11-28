/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSELECTSTATISTIC_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSELECTSTATISTIC_H_

#include <vector>

#include "ioda/ObsSpace.h"

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling ObsFunctionSelectStatistic ObsFunction
class SelectStatisticParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SelectStatisticParameters, Parameters)

 public:
  /// Input variable from which indices will be selected
  oops::RequiredParameter<std::vector<Variable>> variable{"variable", this};

  /// Conditions used to select locations at which the Select Statistic
  /// obsfunction should be applied.
  /// If not specified, all locations will be selected.
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};

  /// Select... max:
  oops::Parameter<bool> selectMax{"select maximum", false, this};
  /// ...min:
  oops::Parameter<bool> selectMin{"select minimum", false, this};
  /// ...mean:
  oops::Parameter<bool> selectMean{"select mean", false, this};
  /// ...median:
  oops::Parameter<bool> selectMedian{"select median", false, this};

  /// If true, then when the input variable is all missing in a record, force the first obs
  /// of the record to still be selected regardless. Required for OPS compatibility.
  oops::Parameter<bool> forceSelect{"force select", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs a vector (int) of all 0's, except 1 at the locations where the input variable
/// fulfil given criteria.
/// \details ObsFunctionSelectStatistic ObsFunction
///
/// Example 1
///
///  obs function:
///    name: IntObsFunction/SelectStatistic
///    options:
///      variable: MetaData/depthBelowWaterSurface
///      select minimum: true
///      select maximum: true
///
/// will return a vector (int) of all 0's, except 1 at the location where
/// MetaData/depthBelowWaterSurface is minimum, i.e. topmost level, and where it is maximum,
/// i.e. bottom-most. If observations are grouped into records, the output will have 1 at the
/// topmost and bottom-most levels of every record. Select criteria (minimum, maximum, mean,
/// median) are all false by default and can be used in any combination.
///
/// The intention is to use SelectStatistic with the Variable Assignment filter to create a
/// variable that can be used in a Perform Action 'where' condition to accept/reject certain obs,
/// e.g. if the topmost and bottom-most levels of every profile must be kept, then any rejection
/// from previously applied filters can be reversed.
///
/// The output of SelectStatistic may also be used as a priority variable to certain Thinning
/// filters, to specify that observations fulfilling certain conditions must not be thinned. If
/// multiple priority levels are required, this can be done by using LinearCombination ObsFunction
/// to scale and combine outputs from multiple calls of SelectStatistic - if this is needed often,
/// and is bothersome to configure in YAML, please raise an issue to modify this filter. Similarly
/// if inputs other than float-type are required.
///
/// Channels are treated separately - e.g. if selecting minimum, there will be a 1 in each channel
/// in each record for which the input variable value is minimum within the channel and record.
///

class SelectStatistic : public ObsFunctionBase<int> {
 public:
  explicit SelectStatistic(const eckit::LocalConfiguration &);

  ~SelectStatistic();

  ufo::ObsAccessor createObsAccessor(const ioda::ObsSpace &obsdb) const;

  void compute(const ObsFilterData &in,
               ioda::ObsDataVector<int> &out) const;

  const ufo::Variables & requiredVariables() const;

 private:
  SelectStatisticParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSELECTSTATISTIC_H_
