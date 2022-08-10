/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_PROFILEAVERAGEOBSPRESSURE_H_
#define UFO_FILTERS_OBSFUNCTIONS_PROFILEAVERAGEOBSPRESSURE_H_

#include <string>
#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the ProfileAverageObsPressure ObsFunction
class ProfileAverageObsPressureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ProfileAverageObsPressureParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> model_vertical_coordinate
    {"model vertical coordinate",
     "Name of the model vertical coordinate.",
     this};

  oops::RequiredParameter<std::string> observation_vertical_coordinate
    {"observation vertical coordinate",
     "Name of the observation vertical coordinate.",
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

/// \brief Fill values of pressure in profiles that have been averaged onto model levels.
///
/// When processing atmospheric profile data, the ObsSpace is extended in order to enable
/// reported profiles to be averaged onto model levels before they are sent to the data
/// assimilation. During the ObsSpace extension procedure, all of the values of air pressure
/// in each averaged profile are initialised to the first entry (i.e. the highest pressure)
/// in the corresponding reported profile.
/// The BackgroundErrorVertInterp operator, which relies on values of air pressure, will therefore
/// not compute background errors correctly for most of the profile.
/// This ObsFunction mitigates that problem by copying the GeoVaL of air pressure located at
/// the first entry in the reported profile into the averaged pressure.
/// Due to the relatively low resolution of the background errors, any horizontal drift during
/// the ascent is ignored.
/// (There is a separate observation operator that computes H(x) values for each simulated
/// variable in the averaged profiles, so a similar ObsFunction is not required in those cases.)
class ProfileAverageObsPressure : public ObsFunctionBase<float> {
 public:
  explicit ProfileAverageObsPressure(const eckit::LocalConfiguration &);
  ~ProfileAverageObsPressure();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ProfileAverageObsPressureParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_PROFILEAVERAGEOBSPRESSURE_H_
