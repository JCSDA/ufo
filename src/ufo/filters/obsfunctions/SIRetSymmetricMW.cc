/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SIRetSymmetricMW.h"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/SIRetMW.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<SIRetSymmetricMW> makerSIRetSymmetricMW_("SIRetSymmetricMW");

// -----------------------------------------------------------------------------

SIRetSymmetricMW::SIRetSymmetricMW(const eckit::LocalConfiguration & conf)
  : invars_(), conf_(conf) {
  SIRetMW siretfunc(conf_);
  ASSERT(siretfunc.siVariableGroups().size() == 2);

  invars_ += siretfunc.requiredVariables();
}

// -----------------------------------------------------------------------------

SIRetSymmetricMW::~SIRetSymmetricMW() {}

// -----------------------------------------------------------------------------

void SIRetSymmetricMW::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();

  // Get SI retrievals from function
  SIRetMW siretfunc(conf_);
  oops::ObsVariables sivars(siretfunc.siVariableGroups());
  ioda::ObsDataVector<float> siret(in.obsspace(), sivars, "ObsFunction", false);
  siretfunc.compute(in, siret);

  // Get symmetric SI amount
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    out[0][iloc] = 0.5 * (siret[0][iloc] + siret[1][iloc]);
    if (siret[0][iloc] >= siretfunc.getBadValue() || siret[1][iloc] >= siretfunc.getBadValue())
        out[0][iloc] = siretfunc.getBadValue();
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SIRetSymmetricMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
