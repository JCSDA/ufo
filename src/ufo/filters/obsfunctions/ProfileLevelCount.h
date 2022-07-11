/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_PROFILELEVELCOUNT_H_
#define UFO_FILTERS_OBSFUNCTIONS_PROFILELEVELCOUNT_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the ProfileLevelCount ObsFunction
class ProfileLevelCountParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ProfileLevelCountParameters, Parameters)

 public:
  /// Conditions used to select locations at which the profile level count
  /// obsfunction should be applied.
  /// If not specified, all locations will be selected.
  oops::Parameter<std::vector<WhereParameters>> where{"where", {}, this};

  /// Operator used to combine the results of successive `where` options at the same location.
  /// The available operators are `and` and `or`.
  oops::Parameter<WhereOperator> whereOperator{"where operator", WhereOperator::AND, this};
};

// -----------------------------------------------------------------------------

/// \brief Count the number of levels in each profile (subject to conditions in the `where` clause).
///
/// This ObsFunction can be used for samples that have been divided into profiles (records).
/// The number of locations in each profile that pass the `where` clause associated
/// with this function is computed. The number is then assigned to all members of the output vector.
/// Note that the `where` clause associated with the governing filter (e.g. Variable Assignment)
/// can be used to control which locations are assigned the count.
class ProfileLevelCount : public ObsFunctionBase<int> {
 public:
  explicit ProfileLevelCount(const eckit::LocalConfiguration &);
  ~ProfileLevelCount();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ProfileLevelCountParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_PROFILELEVELCOUNT_H_
