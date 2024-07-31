/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionArithmetic.h"

#include <string>
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

  // get exponents
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

  // get total log base
  std::string total_log_base;
  if (options_.total_log_base.value() != boost::none)
    total_log_base = options_.total_log_base.value().get();

  // get log bases
  std::vector<std::string> log_bases(nv, "");
  if (options_.log_bases.value() != boost::none)
    log_bases = options_.log_bases.value().get();

  // take absolute value
  std::vector<bool> absolute_value(nv, false);
  if (options_.absolute_value.value() != boost::none)
    absolute_value = options_.absolute_value.value().get();

  // truncate to nearest integer multiple
  std::vector<int> truncate(nv, 0);
  if (options_.truncate.value() != boost::none)
    truncate = options_.truncate.value().get();

  // use channels not obs
  std::vector<int> channels;
  if (options_.useChannelNumber) {
    channels = invars_[0].channels();
    ASSERT(channels.size() > 0);
  }

  // sanity checks
  ASSERT(coefs.size() == nv);
  ASSERT(exponents.size() == nv);
  ASSERT(log_bases.size() == nv);
  ASSERT(absolute_value.size() == nv);
  ASSERT(truncate.size() == nv);
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
  const FunctionValue missing = util::missingValue<FunctionValue>();
  const int missing_int = util::missingValue<int>();
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    ioda::ObsDataVector<FunctionValue> varin(in.obsspace(), invars_[ivar].toOopsObsVariables());
    in.get(invars_[ivar], varin);
    ASSERT(varin.nvars() == out.nvars());
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t ichan = 0; ichan < out.nvars(); ++ichan) {
        if (options_.useChannelNumber) {
          varin[ichan][iloc] = (channels[ichan] == missing_int) ?
                      missing : static_cast<float>(channels[ichan]);
        }
        if ( varin[ichan][iloc] == missing || out[ichan][iloc] == missing ) {
          out[ichan][iloc] = missing;
        } else if (varin[ichan][iloc] < 0 && static_cast<int>(exponents[ivar]) != exponents[ivar]) {
            out[ichan][iloc] = missing;
            oops::Log::warning() << "coefficient exponent "
                                    "Trying to raise a negative number to a non-integer exponent. "
                                    "Output for " << invars_[ivar] << " at location " << iloc <<
                                    " set to missing." << std::endl;
        } else {
          FunctionValue value = varin[ichan][iloc];
          if (absolute_value[ivar]) {
            value = std::abs(value);
          }
          if (truncate[ivar] > 0) {
            value = std::trunc(value / truncate[ivar]) * truncate[ivar];
          }
          out[ichan][iloc] += coefs[ivar] * logpower(
            value,
            exponents[ivar],
            log_bases[ivar]);
          if (ivar == nv - 1) {
            if (out[ichan][iloc] < 0 && static_cast<int>(total_exponent) != total_exponent
                && out[ichan][iloc] != missing) {
                out[ichan][iloc] = missing;
                oops::Log::warning() << "total coefficient exponent "
                                        "Trying to raise a negative number to a non-integer "
                                        "exponent. Output for " << invars_[ivar] <<
                                        " at location " << iloc << " set to missing." << std::endl;
            } else {
                out[ichan][iloc] = total_coeff*logpower(
                  out[ichan][iloc],
                  total_exponent,
                  total_log_base) + intercept;
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
FunctionValue Arithmetic<FunctionValue>::logpower(
    FunctionValue value,
    FunctionValue exponent,
    std::string log_base) const {
  if (log_base.empty()) {
    return power(value, exponent);
  } else {
    // Multiply by exponent for numerical stability
    return exponent * logbasen(value, log_base);
  }
}

template <typename FunctionValue>
FunctionValue Arithmetic<FunctionValue>::logbasen(FunctionValue value, std::string log_base) const {
  if (value <= 0.0) {
    std::stringstream msg;
    msg << "Invalid log value '" << value << "' for log base '" + log_base + "'.";
    throw eckit::BadValue(msg.str(), Here());
  }
  if (log_base == "e") {
    return std::log(value);
  }
  if (std::isdigit(static_cast<unsigned char>(log_base[0]))) {
    FunctionValue base = stof(log_base);
    if (base == 2.0) {
      return std::log2(value);
    } else if (base == 10.0) {
      return std::log10(value);
    } else if (base != 1 && base != 0) {
      return std::log(value) / std::log(base);
    }
  }
  std::stringstream msg;
  msg << "Invalid log base '" + log_base + "' for value '" << value << + "'.";
  throw eckit::BadValue(msg.str(), Here());
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
