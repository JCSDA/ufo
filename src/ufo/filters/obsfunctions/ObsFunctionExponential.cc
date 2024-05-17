/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionExponential.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

static ObsFunctionMaker<Exponential>
                       makerExponential_("Exponential");

// -----------------------------------------------------------------------------

Exponential::Exponential(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  for (const Variable & var : options_.variables.value()) {
    invars_ += var;
  }
}

// -----------------------------------------------------------------------------

Exponential::~Exponential() {}

// -----------------------------------------------------------------------------

void Exponential::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "Compute Exponential ObsFunction." << std::endl;

  // dimension
  const size_t nlocs = in.nlocs();

  // number of input variables
  const size_t nv = invars_.size();

  // get coefficients for exponential function
  const float coeffA = options_.coeffA.value();
  const float coeffB = options_.coeffB.value();
  const float coeffC = options_.coeffC.value();
  const float coeffD = options_.coeffD.value();

  // sanity checks / initialize
  out.zero();
  size_t nchans = out.nvars();

  std::unique_ptr<ioda::Accumulator<size_t>> count_exparg_toobig_accumulator =
    (in.obsspace()).distribution()->createAccumulator<size_t>();

  // compute exponential function of input variable
  ASSERT(nv == 1);  // currently only works with one (possibly multi-channel) x variable
  const float missing = util::missingValue<float>();
  for (size_t ivar = 0; ivar < nv; ++ivar) {
    ioda::ObsDataVector<float> varin(in.obsspace(), invars_[ivar].toOopsObsVariables());
    in.get(invars_[ivar], varin);
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        if ( varin[ichan][iloc] == missing || out[ichan][iloc] == missing ) {
          out[ichan][iloc] = coeffD;
        } else {
          if (coeffB * varin[ichan][iloc] > 40.0) {  // avoid float overflow in arg of exp()
            out[ichan][iloc] = coeffD;
            count_exparg_toobig_accumulator->addTerm(iloc+ichan*nchans, 1);
          } else {
            out[ichan][iloc] = coeffA * exp(coeffB * varin[ichan][iloc]) + coeffC;
          }
        }
      }  // ichan
    }  // nlocs
  }  // nvars
  const std::size_t count_exparg_toobig = count_exparg_toobig_accumulator->computeResult();
  if (count_exparg_toobig > 0) {
    oops::Log::warning() << "Exponential ObsFunction WARNING: " << count_exparg_toobig
      << " instance(s) of exponential argument too large." << std::endl;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & Exponential::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
