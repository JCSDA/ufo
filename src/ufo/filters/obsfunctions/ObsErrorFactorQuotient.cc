/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorQuotient.h"

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorQuotient> makerSteps_("ObsErrorFactorQuotient");

// -----------------------------------------------------------------------------

ObsErrorFactorQuotient::ObsErrorFactorQuotient(const eckit::LocalConfiguration config)
  : invars_() {
  // Initialize options
  options_.deserialize(config);

  // Initialize the two desired variables, numerator and denominator
  const Variable &numerator = options_.numerator.value();
  const Variable &denominator = options_.denominator.value();
  ASSERT(numerator.size() == 1);
  ASSERT(denominator.size() == 1);

  invars_ += numerator;
  invars_ += denominator;

  oops::Log::debug() << "ObsErrorFactorQuotient: config (constructor) = "
                     << config << std::endl;
}

// -----------------------------------------------------------------------------

ObsErrorFactorQuotient::~ObsErrorFactorQuotient() {}

// -----------------------------------------------------------------------------

void ObsErrorFactorQuotient::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue<float>();

  // Get the numeratory and denominator names
  const Variable &numerator = options_.numerator.value();
  const Variable &denominator = options_.denominator.value();
  oops::Log::debug() << "  ObsErrorFactorQuotient, numerator name: " << numerator.variable()
                     << "  and group: " << numerator.group() << std::endl
                     << "                        denominator name: " << denominator.variable()
                     << "  and group: " << denominator.group() << std::endl;

  // Populate the arrays.
  ioda::ObsDataVector<float> numer(data.obsspace(), numerator.toOopsObsVariables());
  data.get(numerator, numer);
  ioda::ObsDataVector<float> denom(data.obsspace(), denominator.toOopsObsVariables());
  data.get(denominator, denom);

    // The 1st index of data should have size 1 and 2nd index should be size nlocs.
  int iv = 0;
  if (numer[iv].size() != out[iv].size() || numer[iv].size() != denom[iv].size()) {
    std::ostringstream errString;
    errString << "Something is wrong, numer size not equal out or denom size."
              << " Sizes: " << numer[iv].size() << " and " << out[iv].size() << std::endl;
    oops::Log::error() << errString.str();
    throw eckit::BadValue(errString.str());
  }

  for (size_t jobs = 0; jobs < numer[iv].size(); ++jobs) {
    out[iv][jobs] = missing;
    if (numer[iv][jobs] == missing || denom[iv][jobs] == missing || denom[iv][jobs] == 0) {
      continue;
    }
    out[iv][jobs] = numer[iv][jobs]/denom[iv][jobs];
  }

  if (options_.save) {
    out.save("DerivedValue");
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorQuotient::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
