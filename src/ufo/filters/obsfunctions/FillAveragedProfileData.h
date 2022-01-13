/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_FILLAVERAGEDPROFILEDATA_H_
#define UFO_FILTERS_OBSFUNCTIONS_FILLAVERAGEDPROFILEDATA_H_

#include <string>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the FillAveragedProfileData ObsFunction
class FillAveragedProfileDataParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FillAveragedProfileDataParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> variable_to_copy
    {"variable to copy",
      "Name of the variable to be copied from the original to the averaged profiles.",
      this};

  oops::RequiredParameter<std::string> observation_vertical_coordinate
    {"observation vertical coordinate",
     "Name of the observation vertical coordinate.",
     this};

  oops::RequiredParameter<std::string> model_vertical_coordinate
    {"model vertical coordinate",
     "Name of the model vertical coordinate.",
     this};

  oops::Parameter<int> numIntersectionIterations{
    "number of intersection iterations",
    "Number of iterations that are used to find the intersection between "
    "the observed profile and each model level",
     3,
     this,
     {oops::minConstraint(1)}};
};

// -----------------------------------------------------------------------------

/// \brief Copy values from an observed atmospheric profile to the same profile
/// averaged onto model levels.
///
/// When processing atmospheric profile data, the ObsSpace is extended in order to enable
/// reported profiles to be averaged onto model levels before they are sent to the data
/// assimilation. During the ObsSpace extension procedure, selected MetaData variables
/// are initialised to the first entry in the corresponding observed profile.
/// This ObsFunction enables a more accurate transfer of values by following the
/// slanted path of the ascent. The SlantPathLocation algorithm is used to do that;
/// it determines the locations in the original profile that correspond to the intersections
/// of the ascent with each model level. Those locations are used to copy values across to
/// the averaged profile.
template <typename FunctionValue>
class FillAveragedProfileData : public ObsFunctionBase<FunctionValue> {
 public:
  explicit FillAveragedProfileData(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  void fillAverageProfile(const ObsFilterData & in,
                          ioda::ObsDataVector<FunctionValue> & out) const;

  FillAveragedProfileDataParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_FILLAVERAGEDPROFILEDATA_H_
