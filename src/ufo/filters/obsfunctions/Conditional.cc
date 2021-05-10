/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/obsfunctions/Conditional.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<Conditional> makerConditional_("Conditional");

Conditional::Conditional(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.validateAndDeserialize(conf);

  // Populate invars_
  for (const LocalConditionalParameters &lcp : options_.cases.value())
    invars_ += getAllWhereVariables(lcp.where);
}

// -----------------------------------------------------------------------------
void Conditional::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  // Assign default value to array
  const float missing = util::missingValue(float());
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
    std::fill(out[ivar].begin(), out[ivar].end(), options_.defaultvalue.value().value_or(missing));
  }  // ivar

  // Assign values based on the where clauses from the configuration.
  // if firstmatchingcase is true, the first case that is true assigns the value.
  // if firstmatchingcase is false, the last matching case will assign the value.
  std::vector<bool> applied(out.nlocs(), false);
  for (const LocalConditionalParameters &lcp : options_.cases.value()) {
    std::vector<bool> apply = processWhere(lcp.where, in);
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
const ufo::Variables & Conditional::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
