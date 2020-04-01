/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionCLWRetMean.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/obsfunctions/ObsFunctionCLWRet.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<ObsFunctionCLWRetMean> makerObsFuncCLWRetMean_("CLWRetMean");

// -----------------------------------------------------------------------------

ObsFunctionCLWRetMean::ObsFunctionCLWRetMean(const eckit::LocalConfiguration & conf)
  : invars_(), conf_(conf) {
  // Initialize options
  options_.deserialize(conf_);

  ObsFunctionCLWRet clwretfunc(conf_);

  invars_ += clwretfunc.requiredVariables();
}

// -----------------------------------------------------------------------------

ObsFunctionCLWRetMean::~ObsFunctionCLWRetMean() {}

// -----------------------------------------------------------------------------

void ObsFunctionCLWRetMean::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimension
  const size_t nlocs = in.nlocs();

  // Get Mean CLW retrievals from function
  oops::Variables clwvars(options_.varGrp.value());
  ioda::ObsDataVector<float> clwret(in.obsspace(), clwvars, "ObsFunction", false);
  ObsFunctionCLWRet clwretfunc(conf_);
  clwretfunc.compute(in, clwret);

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    out[0][iloc] = 0.5 * (clwret[0][iloc] + clwret[1][iloc]);
    if (clwret[0][iloc] >= clwretfunc.getBadValue() || clwret[1][iloc] >= clwretfunc.getBadValue())
        out[0][iloc] = clwretfunc.getBadValue();
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsFunctionCLWRetMean::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
