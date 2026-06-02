/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_REMAPSCANPOSITION_H_
#define UFO_VARIABLETRANSFORMS_CAL_REMAPSCANPOSITION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ufo/variabletransforms/TransformBase.h"

#include "oops/util/parameters/Parameter.h"

namespace ufo {

/// Configuration parameters for the satellite scan position variable transformation
class Cal_RemapScanPositionParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_RemapScanPositionParameters, VariableTransformParametersBase);

 public:
  /// The number of fields of view to remap the current scan position from.  For ATMS this
  /// is 3 and so to preserve backwards compatibility this is the default.
  oops::Parameter<int> numFOV{"number of fields of view", 3, this};
  /// By default, floorRemap is false and the remapping will use the formula
  /// std::ceil(scanpos/numFOV). If floorRemap is True the alternative formula
  /// std::floor((scanpos+1)/numFOV) is used.
  oops::Parameter<bool> floorRemap{"remap to floor", false, this};
};

/*!
* \brief Renumber satellite scan position
*
* \details  Within the Variable Transforms filter, apply the transform "RemapScanPosition"
*  in order to renumber the satellite scan position. At the Met Office ATMS observations are
*  spatially resampled, resulting in 32 fields of view per record sampled from the raw
*  96 FOVs. From the initial observation data the values of MetaData/sensorScanPosition are
*  2, 5, 8, ..., 92, 95 (integers). However, in the calculation of observation bias we
*  require MetaData/sensorScanPosition to be renumbered as 1, 2, 3, ..., 32.
///
*
* See VariableTransformParametersBase for filter setup.
*/
class Cal_RemapScanPosition : public TransformBase {
 public:
  typedef Cal_RemapScanPositionParameters Parameters_;

  Cal_RemapScanPosition(const Parameters_ &options,
                        const ObsFilterData &data,
                        const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                        const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
 private:
  Parameters_ parameters_;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_REMAPSCANPOSITION_H_
