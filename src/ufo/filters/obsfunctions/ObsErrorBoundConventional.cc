/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorBoundConventional.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorBoundConventional> makerSteps_("ObsErrorBoundConventional");
// -----------------------------------------------------------------------------

ObsErrorBoundConventional::ObsErrorBoundConventional(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.reset(new ObsErrorBoundConventionalParameters());
  options_->deserialize(conf);

  const std::string obsvar = options_->obsvar.value();
  invars_ += Variable("ObsErrorData/" + obsvar);
}

// -----------------------------------------------------------------------------

ObsErrorBoundConventional::~ObsErrorBoundConventional() {}

// -----------------------------------------------------------------------------

void ObsErrorBoundConventional::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  size_t nlocs = in.nlocs();

  // Get observation error bounds from options
  const float &obserr_bound_max = options_->obserrBoundMax.value();
  const float &obserr_bound_min = options_->obserrBoundMin.value();
  const float &obserr_bound_factor = options_->obserrBoundFactor.value();
  const float missing = util::missingValue<float>();

  std::vector<float> currentObserr(nlocs);
  const std::string obsvar = options_->obsvar.value();
  in.get(Variable("ObsErrorData/" + obsvar), currentObserr);

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
     out[0][iloc] = obserr_bound_factor * \
         std::clamp(currentObserr[iloc], obserr_bound_min, obserr_bound_max);
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorBoundConventional::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
