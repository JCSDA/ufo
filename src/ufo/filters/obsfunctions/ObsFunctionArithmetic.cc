/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionArithmetic.h"

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

// Note that this class is used for both the Arithmetic and the Linear Combination
// obs functions hence both makers are specified here.
static ObsFunctionMaker<Arithmetic<float>> ArithmeticfloatMaker("Arithmetic");
static ObsFunctionMaker<Arithmetic<int>> ArithmeticintMaker("Arithmetic");

static ObsFunctionMaker<Arithmetic<float>> LinearCombinationfloatMaker("LinearCombination");
static ObsFunctionMaker<Arithmetic<int>> LinearCombinationintMaker("LinearCombination");

// -----------------------------------------------------------------------------

template <typename FunctionValue>
Arithmetic<FunctionValue>::Arithmetic(const eckit::LocalConfiguration & conf)
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
void Arithmetic<FunctionValue>::compute(const ObsFilterData & in,
                                               ioda::ObsDataVector<FunctionValue> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // number of input variables
  const size_t nv = invars_.size();

  // get exponent coefficients
  std::vector<FunctionValue> coefs(nv, 1);
  if (options_.coefs.value() != boost::none)
    coefs = options_.coefs.value().get();

  // get exponent coefficients
  std::vector<FunctionValue> exponents(nv, 1);
  if (options_.exponents.value() != boost::none)
    exponents = options_.exponents.value().get();

  // get total exponent
  FunctionValue total_exponent = static_cast<FunctionValue>(1);
  if (options_.total_exponent.value() != boost::none)
    total_exponent = options_.total_exponent.value().get();

  // get total multiplicative coefficient
  FunctionValue total_coeff = static_cast<FunctionValue>(1);
    if (options_.total_coeff.value() != boost::none)
    total_coeff = options_.total_coeff.value().get();

  // set intercept
  FunctionValue intercept = static_cast<FunctionValue>(0);
  if (options_.intercept.value() != boost::none)
    intercept = options_.intercept.value().get();

  // use channels not obs
  std::vector<int> channels;
  if (options_.useChannelNumber) {
    channels = invars_[0].channels();
    ASSERT(channels.size() > 0);
  }

  // sanity checks
  ASSERT(coefs.size() == nv);
  ASSERT(exponents.size() == nv);
  std::vector<FunctionValue> abs_exponents;
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    abs_exponents.push_back(std::abs(exponents[ivar]));
  }
  if (*std::max_element(abs_exponents.begin(), abs_exponents.end()) > 10
          || std::abs(total_exponent) > 25) {
      oops::Log::warning() << "There is at least one large exponent (>25). "
                              "This may result in overflow errors."
                           << std::endl;
  }

  // initialize
  out.zero();

  // compute linear combination of input variables
  const FunctionValue missing = util::missingValue(missing);
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    ioda::ObsDataVector<FunctionValue> varin(in.obsspace(), invars_[ivar].toOopsVariables());
    in.get(invars_[ivar], varin);
    ASSERT(varin.nvars() == out.nvars());
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t ichan = 0; ichan < out.nvars(); ++ichan) {
        if ( varin[ichan][iloc] == missing || out[ichan][iloc] == missing ) {
          out[ichan][iloc] = missing;
        } else if (varin[ichan][iloc] < 0 && static_cast<int>(exponents[ivar]) != exponents[ivar]) {
            out[ichan][iloc] = missing;
            oops::Log::warning() << "coefficient exponent "
                                    "Trying to raise a negative number to a non-integer exponent. "
                                    "Output for " << invars_[ivar] << " at location " << iloc <<
                                    " set to missing." << std::endl;
        } else {
          if (options_.useChannelNumber) {
            out[ichan][iloc] += coefs[ivar] * power(channels[ichan], exponents[ivar]);
          } else {
            out[ichan][iloc] += coefs[ivar] * power(varin[ichan][iloc], exponents[ivar]);
          }
          if (ivar == nv - 1) {
            if (out[ichan][iloc] < 0 && static_cast<int>(total_exponent) != total_exponent
                && out[ichan][iloc] != missing) {
                out[ichan][iloc] = missing;
                oops::Log::warning() << "total coefficient exponent "
                                        "Trying to raise a negative number to a non-integer "
                                        "exponent. Output for " << invars_[ivar] <<
                                        " at location " << iloc << " set to missing." << std::endl;
            } else {
                out[ichan][iloc] = total_coeff*power(out[ichan][iloc], total_exponent)
                                 + intercept;
            }
          }
        }
      }  // ichan
    }  // nlocs
  }  // nvars
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
const ufo::Variables & Arithmetic<FunctionValue>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
FunctionValue Arithmetic<FunctionValue>::power(FunctionValue value, FunctionValue exponent) const {
  if (exponent == FunctionValue(1.0)) {
    return value;
  } else {
    return std::pow(value, exponent);
  }
}
// -----------------------------------------------------------------------------

}  // namespace ufo
