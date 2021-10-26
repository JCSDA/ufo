/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_AVGKERNEL_OBSAVGKERNELPARAMETERS_H_
#define UFO_AVGKERNEL_OBSAVGKERNELPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// Configuration options recognized by the AvgKernel operator.
class ObsAvgKernelParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAvgKernelParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter <int> nlayers_kernel
    {"nlayers_kernel",
     "Number of layers in the averaging kernel",
     this};

  oops::RequiredParameter <std::vector<double>> ak
    {"ak",
     "First pressure coordinate coefficient. "
     "P(z) = ak(z) + bk(z) * P(surface).",
     this};

  oops::RequiredParameter <std::vector<double>> bk
    {"bk",
     "Second pressure coordinate coefficient. "
     "P(z) = ak(z) + bk(z) * P(surface).",
     this};

  oops::Parameter <std::string> AvgKernelVar
    {"AvgKernelVar",
     "Name of observation variable used in averaging kernel",
     "averaging_kernel",
     this};

  oops::RequiredParameter <std::vector<std::string>> tracerVariables
    {"tracer variables",
     "Names of model tracer variables",
     this};

  oops::Parameter <bool> troposphericColumn
    {"tropospheric column",
     "Perform calculation in troposphere. "
     "An error will be thrown if both this and 'total column' are true.",
     false,
     this};

  oops::Parameter <bool> totalColumn
    {"total column",
     "Perform calculation in total column. "
     "An error will be thrown if both this and 'tropospheric column' are true.",
     false,
     this};

  oops::Parameter <double> modelUnitsCoeff
    {"model units coeff",
     "Conversion between model units",
     1.0,
     this};

  oops::Parameter <double> hofxUnitsCoeff
    {"hofx units coeff",
     "Conversion between H(x) units",
     1.0,
     this};
};

}  // namespace ufo
#endif  // UFO_AVGKERNEL_OBSAVGKERNELPARAMETERS_H_
