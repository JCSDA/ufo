/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_GEOCLOUDCREATECLOUDCOLUMN_H_
#define UFO_FILTERS_OBSFUNCTIONS_GEOCLOUDCREATECLOUDCOLUMN_H_

#include <string>
#include <vector>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling GeoCloudCreateCloudColumn ObsFunction
class GeoCloudCreateCloudColumnParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeoCloudCreateCloudColumnParameters, Parameters)

 public:
  /// Cloud top pressure
  oops::RequiredParameter<std::string> cloudTopPressure
    {"cloud top pressure",
     "Name of cloud top pressure.",
     this};

  /// Cloud amount
  oops::RequiredParameter<std::string> cloudAmount
    {"cloud amount",
     "Name of cloud amount.",
     this};

  /// Cloud optical thickness
  oops::RequiredParameter<std::string> cloudOpticalThickness
    {"cloud optical thickness",
     "Name of cloud optical thickness.",
     this};

  /// Spatial coherence test flag
  oops::RequiredParameter<std::string> spatialCoherenceTest
    {"spatial coherence test",
     "Name of spatial coherence test.",
     this};

  /// Set of 'channels' to use for output of obs error
  /// (really this is levels in the context of GeoCloud)
  oops::RequiredParameter<std::string> channels
    {"channels",
     "Name of channels (levels) for error output.",
     this};

  /// Error for clear levels in cloudy scenes, default value 0.6
  oops::Parameter<float> errorClear
    {"error for cloudy scene clear levels",
     "Error for clear levels in cloudy scenes (default=0.6).",
     0.6f,
     this};

  /// Error for cloudy levels, default value 0.4
  oops::Parameter<float> errorOvercast
    {"error for cloudy levels",
     "Error for cloudy levels (default=0.4).",
     0.4f,
     this};

  /// Error for clear levels in clear scenes, default value 0.6
  oops::Parameter<float> errorClearClear
    {"error for clear scenes",
     "Error for clear scenes (default=0.6).",
     0.6f,
     this};

  /// Observation error multiplier, default value 1.0
  oops::Parameter<float> obsErrorMultiplier
    {"observation error multiplier",
     "Error multiplier for GeoCloud (default=1.0).",
     1.0f,
     this};

  /// GeoCloud altitude limit for clear obs, default 10000.0 hPa
  oops::Parameter<float> altitudeLimitClear
    {"clear altitude limit",
     "GeoCloud clear altitude limit (default=10000.0).",
     10000.0f,
     this};

  /// GeoCloud altitude limit for cloudy obs, default 50000.0 hPa
  oops::Parameter<float> altitudeLimitCloudy
    {"cloudy altitude limit",
     "GeoCloud cloudy altitude limit (default=50000.0).",
     50000.0f,
     this};

  /// Tau scaling denominator, default value 1.0
  oops::Parameter<float> tauScaleDenom
    {"tau denominator",
     "Tau denominator (default=1.0).",
     1.0f,
     this};

  /// Tau scaling offset, default value 0.0
  oops::Parameter<float> tauScaleOffset
    {"tau offset",
     "Tau offset (default=0.0).",
     0.0f,
     this};

  /// NDepth array, default 1 deep everywhere except level closest to surface
  oops::Parameter<std::vector<int>> nDepth
    {"ndepth array",
     "NDepth array.",
     {0},
     this};

  /// Maximum number of layers for calculating cloud optical thickness
  oops::Parameter<int> maxZDepthLayers
    {"maximum zdepth layers",
     "Maximum ZDepth layers (default=1).",
     1,
     this};

  /// Use cloud optical thickness to deepen cloud layer
  oops::Parameter<bool> useCOT
    {"use cloud optical thickness",
     "Use cloud optical thickness? (default=false).",
     false,
     this};

  /// Tau lower height limit array
  oops::Parameter<std::vector<float>> tauZDepthHtLow
    {"tau lower height limit array",
     "Tau lower height limit array.",
     {},
     this};

  /// Tau upper height limit array
  oops::Parameter<std::vector<float>> tauZDepthHtUpp
    {"tau upper height limit array",
     "Tau upper height limit array.",
     {},
     this};

  /// Tau gradient array
  oops::Parameter<std::vector<float>> tauZDepthGrad
    {"tau gradient array",
     "Tau gradient array.",
     {},
     this};

  /// Tau constant array
  oops::Parameter<std::vector<float>> tauZDepthConst
    {"tau constant array",
     "Tau constant array.",
     {},
     this};

  /// Tau maximum depth array
  oops::Parameter<std::vector<float>> tauMaxZDepth
    {"tau maximum depth array",
     "Tau maximum depth array.",
     {},
     this};

  /// Output directory
  oops::Parameter<std::string> outputGroup
    {"output group",
     "Name of output group for obs error (default=DerivedObsError).",
     "DerivedObsError",
     this};
};

class GeoCloudCreateCloudColumn : public ObsFunctionBase<float> {
 public:
    explicit GeoCloudCreateCloudColumn(const eckit::LocalConfiguration &
                                               = eckit::LocalConfiguration());
    void compute(const ObsFilterData &,
                 ioda::ObsDataVector<float> &) const;
    const ufo::Variables & requiredVariables() const;

 private:
    GeoCloudCreateCloudColumnParameters options_;
    ufo::Variables invars_;
    std::vector<int> channels_;
};
}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_GEOCLOUDCREATECLOUDCOLUMN_H_


