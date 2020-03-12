/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFJSFC_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFJSFC_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Output is the error inflation factor as a function of weighted surface temperature
// jacobian and surface emissivity jacobian
// -----------------------------------------------------------------------------


class ObsFunctionErrfJsfc : public ObsFunctionBase {
 public:
  explicit ObsFunctionErrfJsfc(const eckit::LocalConfiguration);
  ~ObsFunctionErrfJsfc();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::string group_;
  std::vector<int> channels_;
 protected:
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONERRFJSFC_H_
