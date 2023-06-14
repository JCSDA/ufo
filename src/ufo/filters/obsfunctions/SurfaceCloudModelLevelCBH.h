/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SURFACECLOUDMODELLEVELCBH_H_
#define UFO_FILTERS_OBSFUNCTIONS_SURFACECLOUDMODELLEVELCBH_H_

#include <string>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling SurfaceCloudModelLevelCBH ObsFunction
class SurfaceCloudModelLevelCBHParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SurfaceCloudModelLevelCBHParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> cloud_base_height
    {"cloud base height",
     "Name of cloud base height.",
     this};
};

class SurfaceCloudModelLevelCBH : public ObsFunctionBase<float> {
 public:
    explicit SurfaceCloudModelLevelCBH(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());
    void compute(const ObsFilterData &,
                       ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    SurfaceCloudModelLevelCBHParameters options_;
    ufo::Variables invars_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SURFACECLOUDMODELLEVELCBH_H_
