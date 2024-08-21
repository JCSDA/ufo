/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELLEVELINDEX_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELLEVELINDEX_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the ModelLevelIndex ObsFunction
class ModelLevelIndexParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelLevelIndexParameters, Parameters)

 public:
  oops::RequiredParameter<std::string>
    obsCoordName{"observation vertical coordinate",
      "Name of observation vertical coordinate",
      this};

  oops::Parameter<std::string>
    obsCoordGroup{"observation vertical coordinate group",
      "Name of observation vertical coordinate group",
      "MetaData", this};

  oops::RequiredParameter<std::string>
    modelCoordName{"model vertical coordinate",
      "Name of model vertical coordinate",
      this};
};

// -----------------------------------------------------------------------------

/// \brief Given observed and model vertical coordinates, assign a model level index
/// to each location.
///
/// \details The user specifies observed and model vertical coordinates.
/// At each location, the model level in which the observation lies is determined.
/// This ObsFunction returns an integer vector of model levels at each location.
///
/// If the observation value is missing, or outside either bound of the corresponding
/// model column, the model level is set to the missing integer value.
///
/// Treatment of values at boundaries is governed by the underlying vertical interpolation
/// code. Examples can be found in the accompanying ctest for this ObsFunction.
class ModelLevelIndex : public ObsFunctionBase<int> {
 public:
  explicit ModelLevelIndex(const eckit::LocalConfiguration &);
  ~ModelLevelIndex();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ModelLevelIndexParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELLEVELINDEX_H_
