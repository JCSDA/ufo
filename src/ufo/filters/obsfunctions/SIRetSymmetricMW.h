/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SIRETSYMMETRICMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_SIRETSYMMETRICMW_H_

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/obsfunctions/SIRetMW.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to the calculation of symmetric cloud amount
///
typedef SIRetMWParameters SIRetSymmetricMWParameters;

///
/// \brief Calculate symmetric (mean) cloud amount from the cloud amount retrieved
/// from the observed and simulated measurements
///
class SIRetSymmetricMW : public ObsFunctionBase<float> {
 public:
  explicit SIRetSymmetricMW(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~SIRetSymmetricMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SIRETSYMMETRICMW_H_
