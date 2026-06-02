/*
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_RADARBEAMGEOMETRY_H_
#define UFO_VARIABLETRANSFORMS_CAL_RADARBEAMGEOMETRY_H_

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

class Cal_RadarBeamGeometryParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_RadarBeamGeometryParameters, VariableTransformParametersBase);

 public:
  oops::Parameter<bool> opsCompatibilityMode
    {"OPS compatibility mode",
     "OPS compatibility mode. "
     "If true, use slightly different values of the effective Earth radius at different points "
     "in the calculation. If false, use a consistent value.",
     false,
     this};
};

/*!
* \brief Radar beam geometry calculations.
*
* Given values of radar beam tilt [deg] and azimuth [deg], gate range [m] and
* station elevation [m], compute the following quantities:
*   sin(tilt)
*   cos(azimuth) * cos(tilt)
*   sin(azimuth) * cos(tilt)
*   gate height [m]
*
* These quantities are required in the radar Doppler winds observation operator.
*
* Please note that tilt is what is known in Europe (in the OPERA standard) as elevation angle.
* Tilt defined as the angle above the horizontal plane of the surface, and
* azimuth is measured clockwise with zero degrees representing due North.

* Example:
*
* \code{.yaml}
* obs filters:
* - filter: Variables Transforms
*   Transform: RadarBeamGeometry
* \endcode
*/

class Cal_RadarBeamGeometry : public TransformBase {
 public:
  typedef Cal_RadarBeamGeometryParameters Parameters_;

  Cal_RadarBeamGeometry(const Parameters_ & options,
                        const ObsFilterData & data,
                        const std::shared_ptr<ioda::ObsDataVector<int>> & flags,
                        const std::shared_ptr<ioda::ObsDataVector<float>> & obserr);
  void runTransform(const std::vector<bool> &apply) override;

 private:
  Parameters_ params_;
};
}  // namespace ufo
#endif  // UFO_VARIABLETRANSFORMS_CAL_RADARBEAMGEOMETRY_H_
