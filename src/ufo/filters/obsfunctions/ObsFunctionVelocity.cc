/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionVelocity.h"

#include <math.h>
#include <algorithm>
#include <set>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<Velocity<float>> floatMaker("Velocity");

// -----------------------------------------------------------------------------

template <typename FunctionValue>
Velocity<FunctionValue>::Velocity(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get parameters from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  group_ = options_.group.value();
  eastwardwindvariable_ = options_.EastwardWindVariable.value();
  northwardwindvariable_ = options_.NorthwardWindVariable.value();
  // Include list of required data from ObsSpace
  invars_ += Variable(group_ + "/" + eastwardwindvariable_, channels_);
  invars_ += Variable(group_ + "/" + northwardwindvariable_, channels_);
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
void Velocity<FunctionValue>::compute(const ObsFilterData & in,
                                      ioda::ObsDataVector<FunctionValue> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // number of input variables
  const size_t nv = invars_.size();

  // number of channels
  const size_t nchans = out.nvars();

  // sanity check we have 2 wind components
  ASSERT(nv == 2);

  // compute wind speed
  const FunctionValue missing = util::missingValue<FunctionValue>();
  ioda::ObsDataVector<FunctionValue> u(in.obsspace(), invars_[0].toOopsObsVariables());
  ioda::ObsDataVector<FunctionValue> v(in.obsspace(), invars_[1].toOopsObsVariables());
  in.get(invars_[0], u);
  in.get(invars_[1], v);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      if (u[ichan][iloc] != missing && v[ichan][iloc] != missing) {
        out[ichan][iloc] = sqrt(u[ichan][iloc]*u[ichan][iloc] + v[ichan][iloc]*v[ichan][iloc]);
      } else {
        out[ichan][iloc] = missing;
      }
    }  // nchans
  }  // nlocs
}

// -----------------------------------------------------------------------------

template <typename FunctionValue>
const ufo::Variables & Velocity<FunctionValue>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
