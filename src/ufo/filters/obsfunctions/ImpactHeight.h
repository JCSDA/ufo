/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_IMPACTHEIGHT_H_
#define UFO_FILTERS_OBSFUNCTIONS_IMPACTHEIGHT_H_

#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {
///
/// \brief Function calculates the GNSS-RO impact height as the difference
/// between the impact parameter and earth's radius of curvature.
///
class ImpactHeight : public ObsFunctionBase {
 public:
  explicit ImpactHeight(const eckit::LocalConfiguration &);
  ~ImpactHeight();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_IMPACTHEIGHT_H_
