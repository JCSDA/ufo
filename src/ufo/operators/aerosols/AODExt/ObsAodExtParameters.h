/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTPARAMETERS_H_
#define UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the AtmVertInterpLay operator.
class ObsAodExtParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAodExtParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<int> Nprofiles
    {"nprofiles",
     "Number of profiles",
     this,
     {oops::minConstraint<int>(2), oops::maxConstraint<int>(3)}};

  oops::RequiredParameter<std::vector<double>> BkgWavelengths
    {"bkg_wavelengths",
     "list of background wavelengths",
     this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODEXT_OBSAODEXTPARAMETERS_H_
