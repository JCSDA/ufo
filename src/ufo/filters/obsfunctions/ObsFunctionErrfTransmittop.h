/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFTRANSMITTOP_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFTRANSMITTOP_H_

#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Error Inflation Factor for satellite radiance as a function of model
// top-to-space transmittance:
// x = model top-to-space transmittance
// EIF = SQRT ( 1.0 / x )
// -----------------------------------------------------------------------------

class ObsFunctionErrfTransmittop : public ObsFunctionBase {
 public:
  explicit ObsFunctionErrfTransmittop(const eckit::LocalConfiguration);
  ~ObsFunctionErrfTransmittop();

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

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFTRANSMITTOP_H_
