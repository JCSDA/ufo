/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_AVERAGETEMPERATUREBELOW_H_
#define UFO_FILTERS_OBSFUNCTIONS_AVERAGETEMPERATUREBELOW_H_

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {
///
/// \brief Options: Height below which to find the average, and whether to
///                 use the layer thickness
///
class AverageTemperatureBelowParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AverageTemperatureBelowParameters, Parameters)

 public:
  /// Limit of the height below which to calculate the average (in meters)
  oops::Parameter<float> heightLimit{"height limit", 20000.0f, this};

  /// If true: use the layer thickness in calculating average.
  /// If false: take average of model level temperatures without considering
  ///           separation between model levels.
  oops::Parameter<bool> useThickness{"use layer thickness", false, this};
};

///
/// \brief Function calculates the average model temperature below a given height.
/// Can be calculated as a simple average over model levels, or take into account
/// the thickness of each model layer.
///
/// Warning: This function assumes that geovals are ordered top-to-bottom.
/// If this is not the case in the provided file, then the function will return
/// missing data everywhere!
///
class AverageTemperatureBelow : public ObsFunctionBase<float> {
 public:
  explicit AverageTemperatureBelow(const eckit::LocalConfiguration &);
  ~AverageTemperatureBelow();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  AverageTemperatureBelowParameters options_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_AVERAGETEMPERATUREBELOW_H_
