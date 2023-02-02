/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORDUPLICATECHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORDUPLICATECHECK_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorDuplicateCheck ObsFunction
class ObsErrorFactorDuplicateCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorDuplicateCheckParameters, Parameters)

 public:
  oops::RequiredParameter<bool> use_air_pressure{"use_air_pressure", this};
  oops::RequiredParameter<std::string> variable{"variable", this};
  oops::Parameter<std::string> original_obserr{"original_obserr", "ObsErrorData", this};
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
};

// -----------------------------------------------------------------------------

/// ### example configurations for a FilterBase derived class: ###
///
///     - filter: BlackList
///       filter variables:
///       - name: windEastward
///       action:
///         name: inflate error
///         inflation variable:
///           name: ObsErrorFactorDuplicateCheck@ObsFunction
///           options:
///             use_air_pres: true
///             variable: windEastward
///
class ObsErrorFactorDuplicateCheck : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorFactorDuplicateCheck";}

  explicit ObsErrorFactorDuplicateCheck(const eckit::Configuration &config);
  ~ObsErrorFactorDuplicateCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
  template <typename T> std::vector<T> getGlobalVariable(const ObsFilterData &,
                                                      const std::string &,
                                                      const std::string &) const;
  int getDispl(const ObsFilterData &, const std::vector<int> &) const;
 private:
  ufo::Variables invars_;
  std::unique_ptr<ObsErrorFactorDuplicateCheckParameters> options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORDUPLICATECHECK_H_
