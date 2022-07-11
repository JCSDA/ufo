/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunction.h"

#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/util/DateTime.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

template <typename FunctionValue>
ObsFunction<FunctionValue>::ObsFunction(const Variable & var)
  : obsfct_(ObsFunctionFactory<FunctionValue>::create(var))
{}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
ObsFunction<FunctionValue>::~ObsFunction() {}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
void ObsFunction<FunctionValue>::compute(const ObsFilterData & in,
                             ioda::ObsDataVector<FunctionValue> & out) const {
  obsfct_->compute(in, out);
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
const ufo::Variables & ObsFunction<FunctionValue>::requiredVariables() const {
  return obsfct_->requiredVariables();
}

// -----------------------------------------------------------------------------

// Explicit instantiations for the supported value types
template class ObsFunction<float>;
template class ObsFunction<int>;
template class ObsFunction<std::string>;
template class ObsFunction<util::DateTime>;

// -----------------------------------------------------------------------------

}  // namespace ufo
