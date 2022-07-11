/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSCATTERING_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSCATTERING_H_

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class ObsFunctionScattering : public ObsFunctionBase<float> {
 public:
  explicit ObsFunctionScattering(const eckit::LocalConfiguration conf
                                       = eckit::LocalConfiguration());
  ~ObsFunctionScattering();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONSCATTERING_H_
