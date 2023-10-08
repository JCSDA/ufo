/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/Conditional.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<Conditional<float>> floatMaker("Conditional");
static ObsFunctionMaker<Conditional<int>> intMaker("Conditional");
static ObsFunctionMaker<Conditional<std::string>> stringMaker("Conditional");
static ObsFunctionMaker<Conditional<util::DateTime>> dateTimeMaker("Conditional");

// -----------------------------------------------------------------------------
template <typename FunctionValue>
Conditional<FunctionValue>::Conditional(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.validateAndDeserialize(conf);

  // Populate invars_
  for (const LocalConditionalParameters<FunctionValue> &lcp : options_.cases.value())
    invars_ += getAllWhereVariables(lcp.where);
}

// -----------------------------------------------------------------------------
template <typename FunctionValue>
void Conditional<FunctionValue>::compute(const ObsFilterData & in,
                                         ioda::ObsDataVector<FunctionValue> & out) const {
  // Assign default value to array
  const FunctionValue missing = util::missingValue<FunctionValue>();
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
    std::fill(out[ivar].begin(), out[ivar].end(), options_.defaultvalue.value().value_or(missing));
  }  // ivar

  // Assign values based on the where clauses from the configuration.
  // if firstmatchingcase is true, the first case that is true assigns the value.
  // if firstmatchingcase is false, the last matching case will assign the value.
  std::vector<bool> applied(out.nlocs(), false);
  for (const LocalConditionalParameters<FunctionValue> &lcp : options_.cases.value()) {
    std::vector<bool> apply = processWhere(lcp.where, in, lcp.whereOperator);
    for (size_t iloc = 0; iloc < out.nlocs(); ++iloc) {
      if (apply[iloc] && applied[iloc] == false) {
        for (size_t ivar = 0; ivar < out.nvars(); ++ivar)
          out[ivar][iloc] = lcp.value.value();
        if (options_.firstmatchingcase.value()) applied[iloc] = true;
      }  // if apply
    }  // iloc
  }  // lcp
}  // compute

// -----------------------------------------------------------------------------
template <typename FunctionValue>
const ufo::Variables & Conditional<FunctionValue>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
