/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONVELOCITY_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONVELOCITY_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

/// \brief Options controlling Velocity ObsFunction
template <typename FunctionValue>
class VelocityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VelocityParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::Parameter<std::string> channelList{"channels", "", this};
  /// Observation group name. Default is ObsValue.
  oops::Parameter<std::string> group{"type", "ObsValue", this};
  /// Variable name of the u(windEastward) component. Default is windEastward.
  oops::Parameter<std::string> EastwardWindVariable{"eastward wind variable",
                                                    "windEastward", this};
  /// Variable name of the v(windNorthward) component. Default is windNorthward.
  oops::Parameter<std::string> NorthwardWindVariable{"northward wind variable",
                                                     "windNorthward", this};
};

// -----------------------------------------------------------------------------

/// \brief Outputs the wind speed of the u(windEastward) and v(windNorthward)
/// components of wind.
///
/// Example 1
///
///  obs function:
///    name: ObsFunction/Velocity
///    options:
///      type: ObsValue
///
/// will return sqrt( ObsValue/windEastward  * ObsValue/windEastward +
///                   ObsValue/windNorthward * ObsValue/windNorthward )
///
/// Example 2 - multi-channel
///
///  obs function:
///    name: ObsFunction/Velocity
///    options:
///      type: ObsValue
///      channels: 1-4
///
/// will return sqrt( ObsValue/windEastward[channel] * ObsValue/windEastward[channel] +
///                   ObsValue/windNorthward[channel] * ObsValue/windNorthward[channel] )
///

template <typename FunctionValue>
class Velocity : public ObsFunctionBase<FunctionValue> {
 public:
  explicit Velocity(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<FunctionValue> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  VelocityParameters<FunctionValue> options_;
  ufo::Variables invars_;
  std::string group_;
  std::vector<int> channels_;
  std::string eastwardwindvariable_;
  std::string northwardwindvariable_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONVELOCITY_H_
