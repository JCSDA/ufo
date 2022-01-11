/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_PSTAR_H_
#define UFO_VARIABLETRANSFORMS_CAL_PSTAR_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/*!
* \brief Calculate pressure at model surface
*
* \details  Within the Variable Transforms filter, apply the transform "PStar" (P*)
*  in order to calculate the pressure at model surface. P* is calculated using one of
*  the three reported pressures (P) along with the corresponding height (Z):
*  - P: Pressure observed at station level, Z: station elevation
*  - P: Pressure at mean sea level, Z = 0m
*  - P: Pressure on standard pressure surface, Z: height of pressure surface.
*  Note that pressures above are listed in order of precedence for which is used to
*  calculate P*, assuming that PstnPrefFlag is true.
*  The calculation of P* is performed in two stages:
*  - Determine the background pressure at the model location.
*  - Calculate the observed pressure at the model surface.
*
*  As well as calculating the observed pressure at the model surface, the flags,
*  observation errors and PGEs for P* are also set.

///
*
* See VariableTransformsParameters for filter setup.
*/
class Cal_PStar : public TransformBase {
 public:
  Cal_PStar(const GenericVariableTransformParameters &options,
            const ObsFilterData &data,
            const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
            const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
  Variables requiredVariables() const override { return gvals_; }

 private:
  Variables gvals_;
};
}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_CAL_PSTAR_H_
