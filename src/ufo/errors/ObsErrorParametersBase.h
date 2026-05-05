/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORPARAMETERSBASE_H_
#define UFO_ERRORS_OBSERRORPARAMETERSBASE_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/errors/ObsErrorReconditioner.h"

namespace ufo {

/// \brief Base obs errors parameters class
class ObsErrorParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsErrorParametersBase, Parameters)
 public:
  /// \brief Name of the covariance model.
  oops::Parameter<std::string> model{"covariance model", "diagonal", this};
  oops::Parameter<double> RMSEtolerance{"Obs Error test tolerance", "RMSE tolerance for"
       "oops::ObsErrorCovariance test", 1.0e-10, this};

  oops::Parameter<ObsErrorReconditionerParameters> reconditioning{"reconditioning",
    ObsErrorReconditionerParameters(), this};
};

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORPARAMETERSBASE_H_
