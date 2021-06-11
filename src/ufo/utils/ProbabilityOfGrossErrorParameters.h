/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PROBABILITYOFGROSSERRORPARAMETERS_H_
#define UFO_UTILS_PROBABILITYOFGROSSERRORPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

  /// \brief Options controlling the operation of the calculations involving
  /// probability of gross error.
  class ProbabilityOfGrossErrorParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ProbabilityOfGrossErrorParameters, Parameters)

   public:  // variables
    /// Maximum value of exponent in background QC.
    oops::Parameter<float> PGE_ExpArgMax{"max exponent", 80.0, this};

    /// PGE rejection limit.
    oops::Parameter<float> PGE_PGECrit{"PGE threshold", 0.1, this};

    /// Multiplication factor for observation errors.
    oops::Parameter<float> PGE_ObErrMult{"obs error multiplier", 1.0, this};

    /// Multiplication factor for background errors.
    oops::Parameter<float> PGE_BkgErrMult{"BG error multiplier", 1.0, this};

    /// Critical value for squared difference from background / ErrVar.
    oops::Parameter<float> PGE_SDiffCrit{"obs minus BG threshold", 100.0, this};
  };
}  // namespace ufo

#endif  // UFO_UTILS_PROBABILITYOFGROSSERRORPARAMETERS_H_

