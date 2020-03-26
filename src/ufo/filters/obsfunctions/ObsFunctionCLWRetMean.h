/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRETMEAN_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRETMEAN_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// /brief Calculate symmetric (mean) cloud amount from the cloud amount retrieved
/// from the observed and simulated measurements
///

class ObsFunctionCLWRetMean : public ObsFunctionBase {
 public:
  explicit ObsFunctionCLWRetMean(const eckit::LocalConfiguration conf
                                       = eckit::LocalConfiguration());
  ~ObsFunctionCLWRetMean();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<std::string> group_;
  std::string funcname_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRETMEAN_H_
