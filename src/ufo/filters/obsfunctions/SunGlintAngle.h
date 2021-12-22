/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SUNGLINTANGLE_H_
#define UFO_FILTERS_OBSFUNCTIONS_SUNGLINTANGLE_H_

#include <string>
#include <vector>

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Calculate Sun glint angles at observation locations.
///
class SunGlintAngle : public ObsFunctionBase<float> {
 public:
  explicit SunGlintAngle(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~SunGlintAngle();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SUNGLINTANGLE_H_
