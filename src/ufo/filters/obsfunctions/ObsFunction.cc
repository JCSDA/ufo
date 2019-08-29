/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunction.h"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsFunction::ObsFunction(const std::string & id)
  : obsfct_(ObsFunctionFactory::create(id))
{}

// -----------------------------------------------------------------------------

ObsFunction::~ObsFunction() {}

// -----------------------------------------------------------------------------

void ObsFunction::compute(const ObsFilterData & in,
                          ioda::ObsDataVector<float> & out) const {
  obsfct_->compute(in, out);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsFunction::requiredGeoVaLs() const {
  return obsfct_->requiredGeoVaLs();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
