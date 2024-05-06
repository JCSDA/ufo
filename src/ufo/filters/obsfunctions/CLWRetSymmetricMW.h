/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLWRETSYMMETRICMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLWRETSYMMETRICMW_H_

#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options applying to the calculation of symmetric cloud amount
///
typedef CLWRetMWParameters CLWRetSymmetricMWParameters;

///
/// \brief Calculate symmetric (mean) cloud amount from the cloud amount retrieved
/// from the observed and simulated measurements
///
class CLWRetSymmetricMW : public ObsFunctionBase<float> {
 public:
  explicit CLWRetSymmetricMW(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~CLWRetSymmetricMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLWRETSYMMETRICMW_H_
