/*
 *
 * Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICEPARAMETERS_H_
#define UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICEPARAMETERS_H_

#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the AodMetOffice operator.
class ObsAodMetOfficeParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAodMetOfficeParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<int> NDustBins{
    "NDustBins", "Number of dust bins", this};
  oops::RequiredParameter<std::vector<double>> AodKExt{
    "AodKExt", "Extinction coefficient for each dust bin", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_AEROSOLS_AODMETOFFICE_OBSAODMETOFFICEPARAMETERS_H_
