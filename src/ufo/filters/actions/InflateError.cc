/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/actions/InflateError.h"

#include <algorithm>
#include <set>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<InflateError> makerInflateErr_("inflate error");

// -----------------------------------------------------------------------------

InflateError::InflateError(const eckit::Configuration & conf)
  : allvars_(), strfactor_(conf.getString("inflation")),
    conf_(conf.getSubConfiguration("options")) {
  if (!isFloat(strfactor_)) {
    allvars_ += Variable(strfactor_, conf_);
  }
}

// -----------------------------------------------------------------------------

void InflateError::apply(const Variables & vars,
                         const std::vector<std::vector<bool>> & flagged,
                         const ObsFilterData & data,
                         ioda::ObsDataVector<int> &,
                         ioda::ObsDataVector<float> & obserr) const {
  oops::Log::debug() << " input obserr: " << obserr << std::endl;
  // Check float was read:
  if (isFloat(strfactor_)) {
    std::vector<float> factors(data.nlocs());
    float factor;
    readFloat(strfactor_, factor);
    oops::Log::debug() << "processing a float: " << factor << std::endl;
    std::fill(factors.begin(), factors.end(), factor);
    ASSERT(factors.size() == obserr.nlocs());
    for (size_t jv = 0; jv < vars.nvars(); ++jv) {
      size_t iv = obserr.varnames().find(vars.variable(jv).variable());
      for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
        if (flagged[iv][jobs]) obserr[iv][jobs] *= factors[jobs];
      }
    }
  // Check string was read:
  } else {
    Variable factorvar(strfactor_, conf_);
    oops::Log::debug() << "processing data: " << strfactor_ << std::endl;
    size_t nfiltervars = vars.size();
    size_t nlocs = data.nlocs();
    // Check channels
    if (conf_.has("channels")) {
      const std::string chlist = conf_.getString("channels");
      std::set<int> channelset = oops::parseIntSet(chlist);
      std::vector<int> channels;
      std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels));
      std::vector<std::string> fvars(channels.size());
      for (size_t ich = 0; ich < channels.size(); ++ich) {
        fvars[ich] = factorvar.variable()+"_"+std::to_string(channels[ich]);
      }
      oops::Variables ovars(fvars);
      ioda::ObsDataVector<float> factors(data.obsspace(), ovars, "ObsFunction", false);
      data.get(factorvar, factors);
      ASSERT(factors.nvars() == vars.nvars());
      for (size_t jv = 0; jv < vars.nvars(); ++jv) {
        size_t iv = obserr.varnames().find(vars.variable(jv).variable());
        for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
          if (flagged[iv][jobs]) obserr[iv][jobs] *= factors[jv][jobs];
        }
      }
    } else {
      std::vector<float> factors(data.nlocs());
      data.get(factorvar, factors);
      ASSERT(factors.size() == obserr.nlocs());
      for (size_t jv = 0; jv < vars.nvars(); ++jv) {
        size_t iv = obserr.varnames().find(vars.variable(jv).variable());
        for (size_t jobs = 0; jobs < obserr.nlocs(); ++jobs) {
          if (flagged[iv][jobs]) obserr[iv][jobs] *= factors[jobs];
        }
      }
    }
  }
  oops::Log::debug() << " output obserr: " << obserr << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
