/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/CLWRetSymmetricMW.h"

#include <algorithm>
#include <cmath>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<CLWRetSymmetricMW> makerCLWRetSymmetricMW_("CLWRetSymmetricMW");

// -----------------------------------------------------------------------------

CLWRetSymmetricMW::CLWRetSymmetricMW(const eckit::LocalConfiguration & conf)
  : invars_(), conf_(conf) {
  CLWRetMW clwretfunc(conf_);
  ASSERT(clwretfunc.clwVariableGroups().size() == 2);

  invars_ += clwretfunc.requiredVariables();
}

// -----------------------------------------------------------------------------

CLWRetSymmetricMW::~CLWRetSymmetricMW() {}

// -----------------------------------------------------------------------------

void CLWRetSymmetricMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();

  // Get CLW retrievals from function
  CLWRetMW clwretfunc(conf_);
  oops::ObsVariables clwvars(clwretfunc.clwVariableGroups());
  ioda::ObsDataVector<float> clwret(in.obsspace(), clwvars, "ObsFunction", false);
  clwretfunc.compute(in, clwret);

  // Get symmetric CLW amount
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    out[0][iloc] = 0.5 * (clwret[0][iloc] + clwret[1][iloc]);
    if (clwret[0][iloc] >= clwretfunc.getBadValue() || clwret[1][iloc] >= clwretfunc.getBadValue())
        out[0][iloc] = clwretfunc.getBadValue();
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & CLWRetSymmetricMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
