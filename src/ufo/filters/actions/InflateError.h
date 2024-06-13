/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_INFLATEERROR_H_
#define UFO_FILTERS_ACTIONS_INFLATEERROR_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class InflateErrorParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(InflateErrorParameters, FilterActionParametersBase);

 public:
  oops::OptionalParameter<float> inflationFactor{"inflation factor", this};
  oops::OptionalParameter<Variable> inflationVariable{"inflation variable", this};

  /// This function is overridden to check that either `inflation factor` or `inflation variable`
  /// is specified, but not both.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;
};

// -----------------------------------------------------------------------------
/// \brief Observation error inflation action.
/// \details Inflates Observation error for filter variables by:
/// - constant (if "inflation factor" is specified in yaml)
/// - spatially varying filter data (if "inflation variable" is specified in yaml).
///   If inflation variable is the same size as filter variables, inflation is done
///   variable by variable (e.g. inflation variable 1 is used to inflate filter
///   variable 1; inflation variable 2 is used to inflate filter variable 2, etc).
///   If inflation variable is of size 1, the same inflation variable is used for
///   updating all filter variables.
class InflateError : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef InflateErrorParameters Parameters_;

  explicit InflateError(const Parameters_ &);

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;

  const ufo::Variables & requiredVariables() const override {return allvars_;}

  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;            /// variables required to compute inflation
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_INFLATEERROR_H_
