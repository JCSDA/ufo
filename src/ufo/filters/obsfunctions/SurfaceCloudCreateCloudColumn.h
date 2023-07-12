/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SURFACECLOUDCREATECLOUDCOLUMN_H_
#define UFO_FILTERS_OBSFUNCTIONS_SURFACECLOUDCREATECLOUDCOLUMN_H_

#include <string>
#include <vector>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling SurfaceCloudCreateCloudColumn ObsFunction
class SurfaceCloudCreateCloudColumnParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SurfaceCloudCreateCloudColumnParameters, Parameters)

 public:
  /// Cloud fraction at base of clouds
  oops::RequiredParameter<std::string> cloudFractionAtBase
    {"cloud fraction at base",
     "Name of cloud fraction at base.",
     this};

  /// Cloud fraction at base of clouds for metar second layer
  oops::RequiredParameter<std::string> cloudFractionAtBaseMetar2
    {"cloud fraction at base for metar second layer",
     "Name of cloud fraction at base for metar second layer.",
     this};

  /// Cloud fraction at base of clouds for metar third layer
  oops::RequiredParameter<std::string> cloudFractionAtBaseMetar3
    {"cloud fraction at base for metar third layer",
     "Name of cloud fraction at base for metar third layer.",
     this};

  /// Cloud fraction at base of clouds for deep convective cloud (3 layers)
  oops::RequiredParameter<std::string> cloudFractionAtBaseDCC1
    {"cloud fraction at base for dcc first layer",
     "Name of cloud fraction at base for deep convective cloud first layer.",
     this};

  oops::RequiredParameter<std::string> cloudFractionAtBaseDCC2
    {"cloud fraction at base for dcc second layer",
     "Name of cloud fraction at base for deep convective cloud second layer.",
     this};

  oops::RequiredParameter<std::string> cloudFractionAtBaseDCC3
    {"cloud fraction at base for dcc third layer",
     "Name of cloud fraction at base for deep convective cloud third layer.",
     this};

  /// Model level height immediately above measured cloud base height
  oops::RequiredParameter<std::string> modelLevelCloudBaseHeight
    {"model level cloud base height",
     "Name of model level cloud base height.",
     this};

  /// Model level height immediately above measured cloud base height
  /// for metar second layer
  oops::RequiredParameter<std::string> modelLevelCloudBaseHeightMetar2
    {"model level cloud base height for metar second layer",
     "Name of model level cloud base height for metar second layer.",
     this};

  /// Model level height immediately above measured cloud base height
  /// for metar third layer
  oops::RequiredParameter<std::string> modelLevelCloudBaseHeightMetar3
    {"model level cloud base height for metar third layer",
     "Name of model level cloud base height for metar third layer.",
     this};

  /// Model level height immediately above measured cloud base height
  /// for deep convective clouds (3 layers)
  oops::RequiredParameter<std::string> modelLevelCloudBaseHeightDCC1
    {"model level cloud base height for dcc first layer",
     "Name of model level cloud base height for deep convective cloud first layer.",
     this};

  oops::RequiredParameter<std::string> modelLevelCloudBaseHeightDCC2
    {"model level cloud base height for dcc second layer",
     "Name of model level cloud base height for deep convective cloud second layer.",
     this};

  oops::RequiredParameter<std::string> modelLevelCloudBaseHeightDCC3
    {"model level cloud base height for dcc third layer",
     "Name of model level cloud base height for deep convective cloud third layer.",
     this};

  /// Set of 'channels' to use for output of obs error
  /// (really this is levels in the context of SurfaceCloud)
  oops::RequiredParameter<std::string> channels
    {"channels",
     "Name of channels (levels) for error output.",
     this};

  /// Error for clear levels, default value 0.55
  oops::Parameter<float> errorClear
    {"error clear",
     "Name of error for clear levels (default=0.55).",
     0.55f,
     this};

  /// Error for cloudy levels, default value 0.25
  oops::Parameter<float> errorOvercast
    {"error overcast",
     "Name of error for cloudy levels (default=0.25).",
     0.25f,
     this};

  /// Observation error multiplier, default value 1.0
  oops::Parameter<float> obsErrorMultiplier
    {"observation error multiplier",
     "Name of error multiplier for SurfaceCloud (default=1.0),",
     1.0f,
     this};

  /// SurfaceCloud altitude limit, default value 5000m
  oops::Parameter<float> altitudeLimit
    {"surfacecloud altitude limit",
     "Name of SurfaceCloud altitude limit (default=5000.0),",
     5000.0f,
     this};

  /// Number of levels to be assigned depth of 1
  /// - these are furthest from the surface
  oops::Parameter<int> numNDepthLevels1
    {"number of one deep levels",
     "Name of number of single-level-depth levels (default=33),",
     33,
     this};

  /// Number of levels to be assigned depth of 2
  oops::Parameter<int> numNDepthLevels2
    {"number of two deep levels",
     "Name of number of two-level-depth levels (default=23),",
     23,
     this};

  /// Number of levels to be assigned depth of 3
  oops::Parameter<int> numNDepthLevels3
    {"number of three deep levels",
     "Name of number of three-level-depth levels (default=5),",
     5,
     this};

  /// Number of levels to be assigned depth of 4
  /// - these are closest to the surface
  oops::Parameter<int> numNDepthLevels4
    {"number of four deep levels",
     "Name of number of four-level-depth levels (default=8),",
     8,
     this};

  /// Output directory
  oops::Parameter<std::string> outputGroup
    {"output group",
     "Name of output group for obs error (default=DerivedObsError).",
     "DerivedObsError",
     this};
};

class SurfaceCloudCreateCloudColumn : public ObsFunctionBase<float> {
 public:
    explicit SurfaceCloudCreateCloudColumn(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());
    void compute(const ObsFilterData &,
                 ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    SurfaceCloudCreateCloudColumnParameters options_;
    std::vector<int> channels_;
    ufo::Variables invars_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SURFACECLOUDCREATECLOUDCOLUMN_H_

