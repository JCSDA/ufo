/*
 * (C) Copyright 2026 NOAA/NCEP/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MODELLEVELSBELOWSURFACE_H_
#define UFO_FILTERS_OBSFUNCTIONS_MODELLEVELSBELOWSURFACE_H_

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling the ModelLevelsBelowSurface ObsFunction
class ModelLevelsBelowSurfaceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelLevelsBelowSurfaceParameters, Parameters)

 public:
  oops::Parameter<Variable> obsPressure{"observation pressure",
    "Observed pressure (Pa)",
    Variable("MetaData/pressure"), this};
};

// -----------------------------------------------------------------------------

/// \brief Number of model levels an observation is below the model surface.
///
/// \details For obs below the model surface (p_obs > p_sfc), returns
/// ln(p_obs/p_sfc) / ln(p_1/p_2), where p_* is the model pressure at that
/// level. This is the depth of the observation below the model surface in
/// model levels. Returns zero for obs above the surface.
/// The variable `rlow` replicates GSI's functionality, quantifying how
/// many model levels an observation lies below the model surface terrain.
///
class ModelLevelsBelowSurface : public ObsFunctionBase<float> {
 public:
  explicit ModelLevelsBelowSurface(const eckit::LocalConfiguration &);
  ~ModelLevelsBelowSurface();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ModelLevelsBelowSurfaceParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MODELLEVELSBELOWSURFACE_H_
