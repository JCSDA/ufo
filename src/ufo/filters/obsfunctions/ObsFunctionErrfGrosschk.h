/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFGROSSCHK_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFGROSSCHK_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Otuput is the residual threshold for BackgroundCheck filter as a function
// of transmittance at model top and latutude.
// Filter out data if | obs -h(x) | > residual threshold
// -----------------------------------------------------------------------------

class ObsFunctionErrfGrosschk : public ObsFunctionBase {
 public:
  explicit ObsFunctionErrfGrosschk(const eckit::LocalConfiguration);
  ~ObsFunctionErrfGrosschk();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::string group_;
  std::vector<int> channels_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFGROSSCHK_H_
