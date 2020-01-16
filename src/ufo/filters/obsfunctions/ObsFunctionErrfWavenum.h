/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFWAVENUM_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFWAVENUM_H_

#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Error Inflation Factor (EIF) for satellite radiance with wavenumber in the
// range of (2000, 2400] during daytime (sun zenith angle < 89) and containing
// water fraction in the footprint:
// x = wavenumber [1/cm]
// y = surface-to-space transmittance
// z = solar zenith angle [radian]
// EIF = SQRT[ 1 / ( 1 - (x - 2000)) * y * MAX(0, COS(z)) / 4000 ]
// -----------------------------------------------------------------------------

class ObsFunctionErrfWavenum : public ObsFunctionBase {
 public:
  explicit ObsFunctionErrfWavenum(const eckit::LocalConfiguration);
  ~ObsFunctionErrfWavenum();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFWAVENUM_H_
