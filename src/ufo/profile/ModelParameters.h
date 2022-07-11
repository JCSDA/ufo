/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_MODELPARAMETERS_H_
#define UFO_PROFILE_MODELPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ufo {

  /// \brief Options related to GeoVaLs used in the profile QC code.
  class ModelParameters : public oops::Parameters {
     OOPS_CONCRETE_PARAMETERS(ModelParameters, Parameters)

   public:
    /// Number of model theta levels. Hard-coded for unit tests.
    size_t numModelLevels() const {return 70;}

    /// Number of model rho levels. Hard-coded for unit tests.
    size_t numModelLevels_rho() const {return 71;}
  };
}  // namespace ufo

#endif  // UFO_PROFILE_MODELPARAMETERS_H_

