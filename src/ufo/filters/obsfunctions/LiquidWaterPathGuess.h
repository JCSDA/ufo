/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_LIQUIDWATERPATHGUESS_H_
#define UFO_FILTERS_OBSFUNCTIONS_LIQUIDWATERPATHGUESS_H_

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Calculate liquid water path from first guess at observation locations.
///  method based on summation of cloud properties in Met Office OPS SatRad routines
///
class LiquidWaterPathGuess : public ObsFunctionBase<float> {
 public:
  explicit LiquidWaterPathGuess(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~LiquidWaterPathGuess();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_LIQUIDWATERPATHGUESS_H_
