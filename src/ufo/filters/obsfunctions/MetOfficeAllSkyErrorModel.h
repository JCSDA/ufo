/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_METOFFICEALLSKYERRORMODEL_H_
#define UFO_FILTERS_OBSFUNCTIONS_METOFFICEALLSKYERRORMODEL_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling MetOfficeAllSkyErrorModel ObsFunction
class MetOfficeAllSkyErrorModelParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MetOfficeAllSkyErrorModelParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> channelList{"channels", this};

  // these coefficients should have size of channels
  oops::OptionalParameter<std::vector<float>> fixland{"fixland", this};
  oops::OptionalParameter<std::vector<float>> fixsea{"fixsea", this};
  oops::OptionalParameter<std::vector<float>> fixice{"fixice", this};
  oops::OptionalParameter<std::vector<float>> taulinsea{"taulinsea", this};
  oops::OptionalParameter<std::vector<float>> taulinland{"taulinland", this};
  oops::OptionalParameter<std::vector<float>> taulinice{"taulinice", this};
  oops::OptionalParameter<std::vector<float>> tausqsea{"tausqsea", this};
  oops::OptionalParameter<std::vector<float>> tausqland{"tausqland", this};
  oops::OptionalParameter<std::vector<float>> tausqice{"tausqice", this};
  oops::OptionalParameter<std::vector<float>> lwpcoef{"lwpcoef", this};
  oops::OptionalParameter<std::vector<float>> iwpcoef{"iwpcoef", this};

  oops::Parameter<float> ScaleVarErrorAboveEmissError{"ScaleVarErrorAboveEmissError",
    "inflate error due to skin temperature over land", 1.0, this};
  oops::Parameter<bool> UseEmissivityAtlas{"UseEmissivityAtlas",
    "scale error when less than MW emissivity error", false, this};
  oops::Parameter<float> maxlwp{"maxlwp", "default lwp value applied to error model", 2.0, this};
  oops::Parameter<float> maxiwp{"maxiwp", "default iwp value applied to error model", 3.0, this};
};

// -----------------------------------------------------------------------------

///
/// \brief Function calculates individual situation-dependent observation errors .
///

class MetOfficeAllSkyErrorModel : public ObsFunctionBase<float> {
 public:
  explicit MetOfficeAllSkyErrorModel(const eckit::LocalConfiguration &);
  ~MetOfficeAllSkyErrorModel();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  MetOfficeAllSkyErrorModelParameters options_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_METOFFICEALLSKYERRORMODEL_H_
