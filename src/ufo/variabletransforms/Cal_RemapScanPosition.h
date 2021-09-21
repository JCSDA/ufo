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

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/*!
* \brief Renumber satellite scan position
*
* \details  Within the Variable Transforms filter, apply the transform "RemapScanPosition"
*  in order to renumber satellite scan position. At the Met Office ATMS observations are
*  spatially resampled, resulting in 32 fields of view per record sampled from the raw
*  96 FOVs. From the initial observation data the values of scan_position@MetaData are 
*  2, 5, 8, ..., 92, 95 (integers). However, in the calculation of observation bias we 
*  require scan_position@MetaData to be renumbered as 1, 2, 3, ..., 32.
///
*
* See VariableTransformsParameters for filter setup.
*/
class Cal_RemapScanPosition : public TransformBase {
 public:
  Cal_RemapScanPosition(const VariableTransformsParameters &options,
                        const ObsFilterData &data,
                        const std::shared_ptr<ioda::ObsDataVector<int>> &flags);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_REMAPSCANPOSITION_H_
