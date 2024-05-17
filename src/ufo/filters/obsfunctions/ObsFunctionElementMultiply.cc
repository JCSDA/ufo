/*
 * (C) Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionElementMultiply.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<ElementMultiply<float>> ElementMultiplyfloatMaker("ElementMultiply");
static ObsFunctionMaker<ElementMultiply<int>> ElementMultiplyintMaker("ElementMultiply");


// -----------------------------------------------------------------------------

template <typename FunctionValue>
ElementMultiply<FunctionValue>::ElementMultiply(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  for (const Variable & var : options_.variables.value()) {
    invars_ += var;
  }
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
void ElementMultiply<FunctionValue>::compute(const ObsFilterData & in,
                                             ioda::ObsDataVector<FunctionValue> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // number of input variables
  const size_t nv = invars_.size();

  // get exponent coefficients
  std::vector<FunctionValue> exponents(nv, 1);
  if (options_.exponents.value() != boost::none)
    exponents = options_.exponents.value().get();

  // get abort if invalid operation
  bool abortIfInvalid = true;
  if (options_.abortIfInvalid.value() != boost::none)
    abortIfInvalid = options_.abortIfInvalid.value().get();

  // sanity checks
  ASSERT(exponents.size() == nv);

  // check for potential overflow due to large exponents
  std::vector<FunctionValue> abs_exponents;
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    abs_exponents.push_back(std::abs(exponents[ivar]));
  }
  if (*std::max_element(abs_exponents.begin(), abs_exponents.end()) > 10) {
      oops::Log::warning() << "There is at least one large exponent (>10). "
                              "This may result in overflow errors."
                           << std::endl;
  }

  // initialize
  out.zero();

  // compute multipliction and/or exponentiation of input variables
  const FunctionValue missing = util::missingValue<FunctionValue>();
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    ioda::ObsDataVector<FunctionValue> varin(in.obsspace(), invars_[ivar].toOopsObsVariables());
    in.get(invars_[ivar], varin);
    ASSERT(varin.nvars() == out.nvars());
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t ichan = 0; ichan < out.nvars(); ++ichan) {
        if ( varin[ichan][iloc] == missing || out[ichan][iloc] == missing ) {
          out[ichan][iloc] = missing;
        } else if (varin[ichan][iloc] < 0 &&
                  static_cast<int>(std::round(exponents[ivar])) != exponents[ivar]) {
            std::stringstream msg;
            msg << "Coefficient exponent error. "
                   "Trying to raise a negative number to a non-integer exponent. "
                   "Complex numbers are not currently supported. "
                   "Variable is " << invars_[ivar] << " at location " << iloc << ".";
            if (abortIfInvalid) {
              throw eckit::NotImplemented(msg.str(), Here());
            } else {
              msg << " Output will be missing.\n";
              oops::Log::error() << msg.str();
              out[ichan][iloc] = missing;
            }
        } else if (varin[ichan][iloc] == 0 && exponents[ivar] < 0) {
            std::stringstream msg;
            msg << "Div/0 error! "
                   "Variable is " << invars_[ivar] << " at location " << iloc << ".";
            if (abortIfInvalid) {
              throw eckit::BadValue(msg.str(), Here());
            } else {
              msg << " Output will be missing.\n";
              oops::Log::error() << msg.str();
              out[ichan][iloc] = missing;
            }
        } else if (ivar == 0) {
            out[ichan][iloc] =  power(varin[ichan][iloc], exponents[ivar]);
        } else {
            out[ichan][iloc] *=  power(varin[ichan][iloc], exponents[ivar]);
        }
      }  // ichan
    }  // nlocs
  }  // nvars
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
const ufo::Variables & ElementMultiply<FunctionValue>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
FunctionValue ElementMultiply<FunctionValue>::power(FunctionValue value, FunctionValue exponent)
  const {
  if (exponent == FunctionValue(1.0)) {
    return value;
  } else {
    return std::pow(value, exponent);
  }
}
// -----------------------------------------------------------------------------

}  // namespace ufo
